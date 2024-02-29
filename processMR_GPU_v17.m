%% processMR_GPU  -- process multi-reflection forces for non-flat lightsails
%
% Applies iterative approximate raytracing appraoch, from each triangle centroid along reflection vector, to find next
% triangle hit by reflected ray, up to MAX_REFL
%
% commented-out lines assist with debugging.
%
% Development notes:
%  Our simulations already calculate the discretized optical thrust/force (and absorption heat loading) resulting from the 
%  interaction of the incident beam with the N triangles comprising our mesh, but now, we want to also model the 
%  additional forces and heating resulting from the light reflected from the lightsail which impacts another region of 
%  the lightsail.  So, the initial starting point for this calculation is a set of up to N rays.  (If a triangle is
%  broken or reverse-illumated, it doesn't produce a reflected array, so the number of starting rays could be less than 
%  N, so we'll say we have M starting rays instead.  If our lightsail is working, though, M should be very close to N.)
%
%% Attempt 1:  Iterative 1-to-N intersection calculations with depth-first handling of reflected rays.  
% I started by writing a vectorized 1-to-N intersection algorithm following several online examples, based on the method
% of Möller and Trumbore ('97).  So we trace 1 ray (origin[3], direction[3]), to N triangles (vert0[3xN], vert1[3xN], 
% vert3[3xN]), producing an output of intersect[1xN, bool] and t[1xN]. 't' is raytricing "time," or the distance
% (as a multiple of the direction vector length) from the origin to the intersection point with the triangle plane, and
% is valid only when 'intersect' is true.  We accomplish this by replicating (repmat) the origin and direction vectors 
% into arrays matching the length the triangle points, then performing all calculation steps by via vectorized operations.  
% 
% There are a few early escapes/short-circuits to escape the full calculations, which I included, but do not expect 
% (and have never observed) to be relevant in any of our simulations.  But they are cheap to include, so why not...  
%
% We then iterate over all the M source rays, which is slow, but not substantially slow compared to the vectorized inner
% operations.  Since we are considering only one source ray at each iteration, it made a lot more sense to continue 
% processing that ray's continued reflection off of other triangles, until the ray escapes, or after MAX_REFL 
% reflections have occurred.  (as oppossed to trying to batch-process each generation of reflected rays via external
% iteration)
%
% During each iterative step, we compute the optical response of the ray's interactions with each impacted triangle, 
% adding the resulting optical force/thrust and absorbed power to that triangle's state vectors via two outputs:  
% t_oth_mr[3xN] (total optical force/thrust from reflected beams)
% t_abs_mr[1xN] (total absorbed power from reflected beams)
% Note that neither contain the effects of the initial beam interaction, which is already processed in the existing
% simulation code
%
% Ultimately this approach was pretty easy to devlop, but it was SLOW due to the added 0(N) x O(MAX_REFL) non-vectorized
% iterative calculations.  Simulating reasonable meshes now takes 10-1000 longer.  Fuck.
%

%% Attempt 2:  Fully vectorized along dimension N.  
% I wanted to try fully vectorizing the problem of M-to-N ray calculations, to get rid of the O(N) outer loop in the V1 
% approach.  But to do this, I had to process vector arrays of length NxM (i.e., O(N^2)).  The simplest approach, which was
% still quite confusing to code, relied on expanding the MxN permutation problem via repmat() and repelem() of the 
% starting arrays, so that each vectorized step of the algorithm would then work on MxN (N^2) length arrays.  I further
% realized that I only had to generate the expanded arrays of triangle vectors/edges once, i.e, repmat()ing each by a 
% factor of M at the first generation of rays.  Thereafter, I could downslect a subsect of these initial arrays to match
% the number of remaining arrays at each reflection iteration.  
%
% Other than the challenge of simply keeping track of the multitude of dimensions, masks, and index arrays needed for 
% these calculations, the one noteworthy challenge encountered here was that, with an entire generation 
% of M-to-N ray calculations contained in a single linear vector array, we would inevitably wind up with some triangles 
% receiving a multitude of reflected rays.  For example, spheres and paraboloids focus many reflected rays from outer/lower
% triangles to a handful of triangles near the top/center of ther sail.  Using simple vectorized-index assignment
% commands to accumulate the resulting absorption and optical forces would effectively ignore the contributions of all 
% but one incident ray at each triangle.  So instead we had to figure out how to use unique() and accumarray() to overcome
% this behavior.  Luckily that didn't seem to slow things down.  So I wound up thinking I had devised a clever and faster
% way to solve the M-to-N ray intersection challenge via full vectorization.  To benchmark it, I had the script 
% repeatedly perform the raytracing operations on a static starting set of rays and triangles, through which I became
% smulgly pleased to see a ~10-100x speed improvement vs. the M-iterated, 1-to-N, depth-first approach.  
%
% But it turned out my benchmarking approach was not representative of the actual simulation environment.  (I guess the 
% MATLAB interpreter must be doing more opimization than I thought it did?)  When I put the M-to_N vectorized code to
% work in the main simulator, it turned out to be several times *slower* than the iterated 1-to-N version!  I guess the
% advantage of full vectorization is entirely lost due to the overhead of fully expanding the input permutation vectors 
% in memory, combined with the staggering memory access delays when working with such huge arrays.  
%
% This was a frustrating conclusion, as the fully vectorized code was dramatacally harder to write than I had
% anticipated.  Perhaps it might be zippy-fast on very small meshes, but that's not the case I needed to improve.  

%% Attempt 3: GPU arrays
% Now that I've been pondering it, I'm pretty sure that compiled code (MEX) would be a great approach for
% this problem, since it would let us eliminate the senseless copying / expansion of the input permutations needed for
% MxN-length vectorization.  All we'd need to fit in cache memory would be the N-length mesh arrays and the M_length ray
% arrays, which we can easily do, and the compiled code could easily work through the permutative arithmetic steps with 
% a handful of M or N-length buffers and a handful of registers.  Writing this concept out here now makes the idea 
% sound even better.
%
% So naturally I decided skip MEX and to try something different, because I'm too lazy to re-learn how to code and 
% debug in MEX, and have delusional aspirations of sharing the code with other collaborators, where the MEX aspect might
% pose a barrier to portability.  
%
% So instead I decided to try out MATLAB's more recent support for GPU-based acceleration.  For the mesh sizes I'm
% contemplating running, I think MEX would actually be faster.  But for extremly large meshes, GPU acceleration
% will eventaully win out due to parallelization and memory bandwith.  More imporantly, the rest of the simulation 
% loop in principle also appears to be well suited for GPU acceleration, and I'm hoping this first step of offloading 
% the raytracing calcs to the GPU might help us towards that ultimate goal in the future.  Finally, whereas I've endured
% the brutal frustration of re-coding and re-optimizing my code for C/MEX in the past, I've never done GPU stuff, so 
% I had further motivation to try something new this time. 
%
%  attempt 3.1:  Just do the ray-triangle intersection, and the subsequent ray reflection on the GPU.  Leave calculation/
%                accumulation of t_oth_mr and t_abs_mr to the CPU, and also compute t_reflp between GPU cycles on the
%                CPU.  
%
%  Result on moderate sizes array:   numtris = 3456     numnodes = 1801       numedges = 5256
%
% db_bm_timer_gpu
% db_bm_timer_gpu =
%     0.8081
% db_bm_timer_vec
% db_bm_timer_vec =
%    77.8848
% db_bm_timer_itr
% db_bm_timer_itr =
%   120.4719
%
% Fucking FINALLY this is fast enough to run our simulations!

%% Failed attempt 4:  Make the whole simulation run in GPU arrays...  
%  this was a pain in the ass.  I converted the entire main simulation loop to run using GPUArrays, and it was slower!  
%  the code has been reverted back to using the vectorized (non-gpu) matlab functions.



%% optional debug setup and figure generation:
% close all hidden;  % close prior figure

% simulate_v17;  % for debug, add a "return" in the simulate code just after the light has been turned on, after trinorms
%                 % and specular optics have been calculated.  That way we can test and debug the MR-RT code on a fresh
%                 % set of valid data.
% t_ev10 = t_ev1;
% t_ev20 = t_ev2;
% t_refln0 = t_refln;
% t_vert00 = t_vert0;

showplots = 0;  % turn on for debugging plots

plotcolormode = 0;  % 1 == num hits per tri, other == power absorbed per tri                                     %
%if showplots, close all, end                                                                                    %
if showplots                                                                                                     %
    if ~exist('bjet'),  load('bjet.mat'); end                                                                    %
    if plotcolormode == 1                                                                                        %
        cdatalabel = 'Number of hits';                                                                           %
    else                                                                                                         %
        cdatalabel='Multibounce absorbed power (W)';                                                             %
    end                                                                                                          %
    howmanyplots = 0;                                                                                            %
end                                                                                                              %



%% benchmarker
% ?







% if 1
    %%  GPU version
    %
    % first set up debug plot window if needed:
    if showplots                                                                                                                 %
        prettycolors = { 'g' 'r' 'c' 'm' 'y' };                                                                                  %
        figure                                                                                                                   % 
        hold on;                                                                                                                 %
        %    quiver3( t_cx, t_cy, zreliefmag*t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), zreliefmag*t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead', 'off');
        %  if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end  %plot mechanical force
        %  vectors (can be confusing).

        quiver3(t_cx, t_cy, t_cz, t_refln(1,:), t_refln(2,:), t_refln(3,:), 'g');                                                %
        hold on;                                                                                                                 %
        xlabel('x (mm)');                                                                                                        %
        ylabel('y (mm)');                                                                                                        %
        zlabel('z (mm)');                                                                                                        %
        axis equal;                                                                                                              %
    end                                                                                                                          %
    
    
    
    
    %tri_broken2orig_map = 1:numtris;  % fixme, move to main sim code, retain as state vector
    
    
%     %FIXME FIXME temporary debug, let's break some triangles:
     t_notbroken(:) = true;
     anybroken = 0;
%    t_notbroken(1:10:end) = 0;
%    anybroken = 1;
    
    
    tic;
    
    ray_pwrs = t_reflp(:, t_notbroken);
    
    % now initialize GPU arrays
    if anybroken
        bigvecsize = [sum(t_notbroken) 3];
        gpu_ev1 =      gpuArray(  t_ev1(:,t_notbroken)');
        gpu_ev2 =      gpuArray(  t_ev2(:,t_notbroken)');
        gpu_v0 =       gpuArray(t_vert0(:,t_notbroken)');
        gpu_t_norms =  gpuArray(t_norms(:,t_notbroken) ); %not used in arrayfun
        gpu_ray_dirs = gpuArray(t_refln(:,t_notbroken) );
        tri_broken2orig_map = thisseemsdumb(t_notbroken);
        gpu_oxyzs = gpuArray([ t_cx(t_notbroken);  t_cy(t_notbroken);  t_cz(t_notbroken)]);  

    else
        bigvecsize = fliplr(size(t_ev1));
        gpu_ev1 =      gpuArray(  t_ev1'); % const
        gpu_ev2 =      gpuArray(  t_ev2'); % const
        gpu_v0 =       gpuArray(t_vert0'); % const
        gpu_t_norms =  gpuArray(t_norms ); %not used in arrayfun, const
        gpu_ray_dirs = gpuArray(t_refln ); % varies
        gpu_oxyzs = gpuArray([ t_cx;  t_cy;  t_cz]);

    end
    
    t_mr_nh = zeros(size(t_na));
    t_oth_mr = zeros(size(t_oth));
    t_abs_mr = zeros(size(t_na));
    mr_nhits = zeros(1,MAX_REFL);
    
    %
    % gpu_ev1 = gpuArray(t_ev1');
    % gpu_ev2 = gpuArray(t_ev2');
    % gpu_v0 = gpuArray(t_vert0');
    % gpu_tnorms = gpuArray(t_norms);
    
    % gpu_ev1 = gpuArray(t_ev1' + rand(bigvecsize));
    % gpu_ev2 = gpuArray(t_ev2' + rand(bigvecsize));
    % gpu_v0 = gpuArray(t_vert0' + rand(bigvecsize));
    
    
    % this way is no faster:
    % gpu_ev1 = gpuArray(t_ev1);
    % gpu_ev2 = gpuArray(t_ev2);
    % gpu_v0 = gpuArray(t_vert0);
    
    
    %% benchmarking code setup:
    %nbm = 1;
    %for nfucks=1:nbm  % benchmarking loop, generates new slightly tweaked values
        
        %ntest  = length(t_na)-1000;  % subtract some initial rays to ensure non-square MxN, to help detect dimension confusion during development
   %     ntest  = length(t_na);
        
%         gpu_ox = gpuArray(t_cx(1:ntest)+(nfucks>1).*rand(1,ntest)); %fixme delete
%         gpu_oy = gpuArray(t_cy(1:ntest)+(nfucks>1).*rand(1,ntest)); %fixme delete
%         gpu_oz = gpuArray(t_cz(1:ntest)+(nfucks>1).*rand(1,ntest)); %fixme delete
%         gpu_oxyzs = [ gpu_ox; gpu_oy; gpu_oz]; % FIXME this is temporary, the correct start value is given above, just adding some noise for benchmarking
%         
%         gpu_dx = gpuArray(t_refln(1,(1:ntest))+(nfucks>1).*rand(1,ntest)); %fixme delete
%         gpu_dy = gpuArray(t_refln(2,(1:ntest))+(nfucks>1).*rand(1,ntest)); %fixme delete
%         gpu_dz = gpuArray(t_refln(3,(1:ntest))+(nfucks>1).*rand(1,ntest)); %fixme delete
%         gpu_ray_dirs = [gpu_dx; gpu_dy; gpu_dz];  % FIXME this is temporary, for debugging
        %% end of benchmarking code setup
        
        
        ri = 1;
        %% start multi ray loop here
        while 1
            
            % This is where I'd cull rays below a certain threshold, but I don't think it would be worth the effort.
            
            gpu_rn = gpuArray(1:size(gpu_oxyzs,2));  % integers counting from 1 to M, the number of rays
            
            if ri == 1
                [gpu_intersect, gpu_t_all] = arrayfun( @RTI_MK_GPU, ...
                    ...gpu_ox(nray), gpu_oy(nray), gpu_oz(nray), ...
                    gpu_oxyzs(1,:),    gpu_oxyzs(2,:),    gpu_oxyzs(3,:), ...
                    ...ray_dx(nray),ray_dy(nray),ray_dy(nray), ...
                    gpu_ray_dirs(1,:), gpu_ray_dirs(2,:), gpu_ray_dirs(3,:), ...
                    gpu_v0(:,1),       gpu_v0(:,2),       gpu_v0(:,3), ...
                    gpu_ev1(:,1),      gpu_ev1(:,2),      gpu_ev1(:,3), ...
                    gpu_ev2(:,1),      gpu_ev2(:,2),      gpu_ev2(:,3) );
            end
            % alternate version, didnt run any faster...
            % [gpu_intersect, gpu_t] = arrayfun( @RTI_MK_GPU, ...
            %     ...gpu_ox(nray), gpu_oy(nray), gpu_oz(nray), ...
            %     gpu_ox,         gpu_oy,         gpu_oz, ...
            %     ...ray_dx(nray),ray_dy(nray),ray_dy(nray), ...
            %     gpu_dx,         gpu_dy,         gpu_dz, ...
            %     gpu_v0(1,:)',   gpu_v0(2,:)',   gpu_v0(3,:)', ...
            %     gpu_ev1(1,:)',  gpu_ev1(2,:)',  gpu_ev1(3,:)', ...
            %     gpu_ev2(1,:)',  gpu_ev2(2,:)',  gpu_ev2(3,:)' );
            
            
            % intersect:  Num rows = num tris-nb;     num cols = num rays.   True if col index (ray) intersects row index (tri-nb)
            
            [gpu_m, gpu_i] = max(gpu_intersect);  %length of m and i:  num starting rays.
            % m is a bool mask--whether or not that ray hit.  Its index is the starting ray
            % index. its sum is the number of hit rays.
            % i the triangle index (nb) of the detected hit.  i(k) is inavlid when m(k) == 0.
            %  we will use this to index the triangles hit by the ray
             
            
            
            gpu_hit_src_rn = gpu_rn(gpu_m);  % length is number of hit rays, values span index range of starting rays, skips indices of missed rays
            gpu_hit_dest_tn = gpu_i(gpu_m);   % length is number of hit rays, values span index range of non-broken tris for this simstep
            gpu_ar_costheta = dot(gpu_t_norms(:,gpu_hit_dest_tn), gpu_ray_dirs(:,gpu_m)); % ok this is one of our outputs.  Already trimmed to length of hit rays.
            % TODO:  Further trim entire output set to filter out negative cosine angles????
            gpu_ar_theta = acos(gpu_ar_costheta);
            
            
            
            
            
            gpu_hit_array_idxs =  sub2ind( size(gpu_t_all), gpu_hit_dest_tn, gpu_hit_src_rn);
            gpu_hit_t = gpu_t_all(gpu_hit_array_idxs); % length is number of hit rays
            %gpu_all_array_idxs =  sub2ind( size(gpu_t_all), gpu_i, gpu_rn);  % this is the same as the 'hit' version but includes non-hit intersections if present, for debug only, deleteme
            %gpu_all_t = gpu_t_all(gpu_all_array_idxs);
            
            
            
            
            %% OK, done with GPU calcs, let's pull data back to CPU
            %intersect = gather(gpu_intersect);
            % t = gather(gpu_t_all);
            %  M = gather(gpu_m);
            %  I = gather(gpu_i);
            %  T = gather(gpu_t);
            %  ResultSrcIdxs = gather(gpu_hit_src_rn);
            %  ResultArrIdxs = gather(gpu_hit_dest_tn);
            %  ResultTs = gather(gpu_hit_t);
            % ResultAllTs = gather(gpu_all_t);
            
            % first let's get a 1D array, that'll be quick.
            Result_hit_tri_idxs_nb = gather(gpu_hit_dest_tn);
            ray_dest_idxs = tri_broken2orig_map(Result_hit_tri_idxs_nb);
            
            % check for early escape option, if all rays failed:
            if isempty(ray_dest_idxs) %if that's empty, we're done!
                break;
            end
            mr_nhits(ri) = length(ray_dest_idxs);  % record number of hit rays
            
            % still here?  I guess we should bring back the rest of the data
            Result_remain_dirs = gather(gpu_ray_dirs(:,gpu_m));
            Result_costheta = gather(gpu_ar_costheta);
            Result_theta = gather(gpu_ar_theta);
            %Result_theta = acos(Result_costheta);
            Result_mask = gather(gpu_m);

            
            %% Bring these back from the GPU, but only if only needed for plot/debug:
            if showplots
                Result_src_xyzs = gather(gpu_oxyzs(:,gpu_m));
                gpu_dest_xyzs = gpu_oxyzs(:,gpu_m) + gpu_hit_t.* gpu_ray_dirs(:,gpu_m);
                Result_dest_xyzs = gather(gpu_dest_xyzs);
            end
            
            %% Now we can go back and update the gpu arrays for the next iteration, if needed
            if (ri >= MAX_REFL )
                %Don't need to update GPU arrays, this is the last iteration
            else % I hope these commands might run in parallel on the GPU while we process the output on the CPU, that's why I'm initiating them here instead of at end-of-loop
                %if ~showplots, 
                    gpu_dest_xyzs = gpu_oxyzs(:,gpu_m) + gpu_hit_t.* gpu_ray_dirs(:,gpu_m); 
                %end
                % [t_cx(ray_dest_idxs) + hit_ts.* ray_dirs(1) ;   ...
                %  t_cy(ray_dest_idxs) + hit_ts.* ray_dirs(2) ;   ...
                %  t_cz(ray_dest_idxs) + hit_ts.* ray_dirs(3) ];
                gpu_refl_dirs = gpu_ray_dirs(:,gpu_m) - 2.* gpu_ar_costheta .* gpu_t_norms(:,gpu_hit_dest_tn);  %we could move this to below the escape if below to save a bit of time
                gpu_oxyzs = gpu_dest_xyzs;
                gpu_ray_dirs = gpu_refl_dirs;
                
                % trigger the next big M-to-Raytracing calculation... again I'm hoping this will start on the GPU while
                % we perform our CPU output calcs below...
                [gpu_intersect, gpu_t_all] = arrayfun( @RTI_MK_GPU, ...
                    ...gpu_ox(nray), gpu_oy(nray), gpu_oz(nray), ...
                    gpu_oxyzs(1,:),    gpu_oxyzs(2,:),    gpu_oxyzs(3,:), ...
                    ...ray_dx(nray),ray_dy(nray),ray_dy(nray), ...
                    gpu_ray_dirs(1,:), gpu_ray_dirs(2,:), gpu_ray_dirs(3,:), ...
                    gpu_v0(:,1),       gpu_v0(:,2),       gpu_v0(:,3), ...
                    gpu_ev1(:,1),      gpu_ev1(:,2),      gpu_ev1(:,3), ...
                    gpu_ev2(:,1),      gpu_ev2(:,2),      gpu_ev2(:,3) );

            end
            
            %% Back at the ranch, we can do the rest of our CPU calcs
            % downselect ray power array to the successful rays only
            ray_pwrs = ray_pwrs(Result_mask);

            %tri_map_indxs = tri_broken2orig_map( ... % this is already Result_hit_tri_idxs
            
            % perform unique() stuff so we can map the rays back to triangles in an accumulator
            [ray_dest_idxs_unique, ~, idxUnique] = unique(ray_dest_idxs);
            t_mr_nh(ray_dest_idxs_unique) = t_mr_nh(ray_dest_idxs_unique) + accumarray(idxUnique, 1)';  % keep track of number of rays reaching each triangle
            
            % do LUL using GPU calculated values for costheta.  TODO:  compare GPU vs CPU calculation of acos?
            
            %ok now calculate reflection coefficients for the re-reflected rays
            % default material:
                ar_abs = t_Iabs(ray_dest_idxs) .* ray_pwrs ; %watt
                ar_rpwr = ray_pwrs .* t_Irefl(ray_dest_idxs) ; % watt
                
                %ray_pwrs = ray_pwrs .* Irefl;


            if useSpecularTMM  % now calculate optical response and assign forces/power to triangles
                ar_LULi = max(min(floor(Result_theta ./ LUL.angleStep)+1,LUL.num),1);
                mymask = t_mat( ray_dest_idxs ) > 0;
                ar_abs( mymask ) = LUL.A(  sub2ind(LUL.sz, t_mat(mymask), ar_LULi) )  .*  ray_pwrs;
                ar_rpwr( mymask ) = LUL.R( sub2ind(LUL.sz, t_mat(mymask), ar_LULi) )  .*  ray_pwrs;
                
                %ar_abs = ar_specA  .* ar_specR ;;
                %ar_rpwr = ray_pwrs .* ar_specR ;
                
                %ray_pwrs = ray_pwrs .* ar_specR;
            else
                % %ar_abs = t_Iabs(Result_mask) .* ray_pwrs ; %watt
                % ar_abs = Iabs .* ray_pwrs ; %watt
                % %ar_rpwr = ray_pwrs .* t_Irefl(Result_mask) ; % watt
                % ar_rpwr = ray_pwrs .* Irefl ; % watt
                % 
                % ray_pwrs = ray_pwrs .* Irefl;
            end
            ray_pwrs = ar_rpwr;
            
            if ~isempty(idxUnique)
                t_abs_mr(ray_dest_idxs_unique) = t_abs_mr(ray_dest_idxs_unique) + accumarray(idxUnique, ar_abs)'; % Watt
                
                % calculate optical force ("optical thrust") for each triangle
                t_oth_mr(1,ray_dest_idxs_unique) = t_oth_mr(1,ray_dest_idxs_unique) + ...
                    accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(1,ray_dest_idxs) .* Result_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* Result_remain_dirs(1,:)    )';    % and this is the absorption force
                t_oth_mr(2,ray_dest_idxs_unique) = t_oth_mr(2,ray_dest_idxs_unique) + ...
                    accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(2,ray_dest_idxs) .* Result_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* Result_remain_dirs(2,:)    )';    % and this is the absorption force
                t_oth_mr(3,ray_dest_idxs_unique) = t_oth_mr(3,ray_dest_idxs_unique) + ...
                    accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(3,ray_dest_idxs) .* Result_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* Result_remain_dirs(3,:)    )';    % and this is the absorption force
            else
                t_abs_mr(ray_dest_idxs)   = t_abs_mr(ray_dest_idxs) + ar_abs; % Wattt_oth_mr(1,ray_dest_idxs_unique) + ...
                t_oth_mr(1,ray_dest_idxs) = t_oth_mr(1,ray_dest_idxs) + ...
                    2 .* ar_rpwr ./ c0 .* t_norms(1,ray_dest_idxs) .* Result_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* Result_remain_dirs(1,:) ;    % and this is the absorption force
                t_oth_mr(2,ray_dest_idxs) = t_oth_mr(2,ray_dest_idxs) + ...
                    2 .* ar_rpwr ./ c0 .* t_norms(2,ray_dest_idxs) .* Result_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* Result_remain_dirs(2,:) ;    % and this is the absorption force
                t_oth_mr(3,ray_dest_idxs) = t_oth_mr(3,ray_dest_idxs) + ...
                    2 .* ar_rpwr ./ c0 .* t_norms(3,ray_dest_idxs) .* Result_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* Result_remain_dirs(3,:) ;    % and this is the absorption force
            end
            
            if showplots                                                                                                                             %
                mycolor = prettycolors{ mod(ri-1, 5)+1 };                                                                                            %
                    for npl = 1:length(ray_dest_idxs)                                                                                                %
                        plot3( [Result_src_xyzs(1,npl)  Result_dest_xyzs(1,npl)], ...                                                                %
                               [Result_src_xyzs(2,npl)  Result_dest_xyzs(2,npl)], ...                                                                %
                               [Result_src_xyzs(3,npl)  Result_dest_xyzs(3,npl)], mycolor , 'LineWidth', ri, 'Marker', '.');                         %
                        howmanyplots = howmanyplots + 1;                                                                                             %
                    end                                                                                                                              %
                disp(['Step ' num2str(ri) ' numplots = ' num2str(howmanyplots)]);                                                                    %
            end                                                                                                                                      %
             
           
            
            ri = ri+1;
            
            
            if ri > MAX_REFL
                break;
            end
            
        end
        
 %   end
%     toc
% disp(' '); disp(' ');    
% disp(['GPU code: ' num2str(toc) ' s']);
% disp(num2str(mr_nhits))
% disp(['Force: ' num2str([ sum(t_oth_mr(1,:)) sum(t_oth_mr(2,:)) sum(t_oth_mr(3,:)) ]) ]);
% disp(['Absorb: ' num2str(sum(t_abs_mr))]);
% %disp(['Number of lines rendered: ' num2str(howmanyplots) ] );
% disp(' ');


% db_bm_timer_gpu_val = toc;
 

if showplots
    if plotcolormode == 1
        cdata=t_mr_nh;
    else
        cdata=t_abs_mr;
    end
    trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
        set(gca,'Colormap',bjet);
    hcolorbar = colorbar;
    %hcolorbar.Label.String = '#ray hits per triangle';
    %hcolorbar.Label.String = 'Multibounce absorbed power (W)';
    hcolorbar.Label.String = cdatalabel;
    view([0 0]);
    axis tight;
    set(gca,'YLim',[0 2*radiusmm]);
    maxoth = max(vecnorm(t_oth));
    maxoth_mr = max(vecnorm(t_oth_mr));
    mr_ratio = maxoth_mr / maxoth;
    quiver3(t_cx, t_cy, t_cz, t_oth(1,:), t_oth(2,:), t_oth(3,:), 'm');  %incident thrust
    quiver3(t_cx, t_cy, t_cz, t_oth_mr(1,:), t_oth_mr(2,:), t_oth_mr(3,:), mr_ratio, 'b'); %multi-reflection thrust
    title('GPU');
end


%     % end of GPU version
% 
% end
% 
% %return; % return here to skip other versions
% 
% 




















% 
% if 1
% if showplots
% prettycolors = { 'r' 'g' 'c' 'm' 'y' };
% figure
% hold on;
% %    quiver3( t_cx, t_cy, zreliefmag*t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), zreliefmag*t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead', 'off');
% %  if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end  %plot mechanical force
%   %  vectors (can be confusing).
% 
% quiver3(t_cx, t_cy, t_cz, t_refln(1,:), t_refln(2,:), t_refln(3,:), 'g');
% hold on;
% xlabel('x (mm)');
% ylabel('y (mm)');
% zlabel('z (mm)');
% axis equal;
% end       
%          howmanyplots = 0;       
% %n_xyz = [n_x' n_y' n_z'];  % node position vector
% 
% %t_cxyz = [t_cx' t_cy' t_cz'];  %ray's position vector
% 
% % ray's direction:    % t_refln'
% numtests=1;
% 
% startrealtic = tic;
% 
%     
% eps        = 1e-8;
% t_indxs = 1:numtris;
% numtrisgood = sum(t_notbroken);
% numrays = numtrisgood;
% src_tris  = t_indxs(t_notbroken);
% mr_nhits = zeros(1,MAX_REFL);
% %t_mr_na = zeros(1, length(t_na));
% t_mr_nh = zeros(size(t_na));  %how many rays reach each tri?
% t_oth_mr = zeros(size(t_oth));
% t_abs_mr = zeros(size(t_na));
% 
% tic
% mytimer2 = 0;
% mytimer = 0;
% for ntest=1:numtests
% 
% t_vert0f = [n_x(t_na(t_notbroken)); n_y(t_na(t_notbroken)); n_z(t_na(t_notbroken))];
% t_ev1f = t_ev1(:,t_notbroken);
% t_ev2f = t_ev2(:,t_notbroken);
% 
% t_cxyz = [t_cx(t_notbroken); t_cy(t_notbroken); t_cz(t_notbroken)];
% 
% ray_origins = t_cxyz;
% ray_dirs = t_refln(:, t_notbroken);
% ray_pwrs = t_reflp(:, t_notbroken);
% 
% t_vert0_exp = repmat(t_vert0f, 1, numrays);
% t_ev1_exp = repmat(t_ev1f, 1, numrays );
% t_ev2_exp = repmat(t_ev2f, 1, numrays );
% 
% numcombs = numrays*numtrisgood;
% 
% %intersect = false(1, numcombs); % by default there are no intersections
% t = inf+zeros(1, numcombs); u=t; v=t;
% tvec = zeros(3,numcombs);
% pvec = zeros(3,numcombs);
% det = zeros(1,numcombs);
% 
% 
%     
% ray_origins = t_cxyz;
% ray_dirs = t_refln(:, t_notbroken);
% ray_pwrs = t_reflp(:, t_notbroken);
%     
% startrttic = tic;
% ri = 1;
% 
% 
% % start loop somehow
% while 1
%     % ok, replicate the arrays
%     numcombs = numrays*numtrisgood;
%     ray_origins_exp = repelem(ray_origins, 1, numtrisgood);
%     ray_dirs_exp = repelem(ray_dirs, 1, numtrisgood);
%     
%     %  disp(['ORI ' num2str(size(ray_origins_exp))]);
%     %  disp(['DIR ' num2str(size(ray_dirs_exp))]);
%     %  disp(['vert0 ' num2str(size(t_vert0_exp))]);
%     %  disp(['edg1 ' num2str(size(t_ev1_exp))]);
%     %  disp(['edg2 ' num2str(size(t_ev2_exp))]);
%     startsubtic = tic;
%     
%     %intersect(1:numcombs) = false; % by default there are no intersections
%     %t(1:numcombs) = inf; 
%     u(1:numcombs)=inf;
%     v(1:numcombs)=inf;
%     tvec(:,1:numcombs)  = ray_origins_exp(:,1:numcombs) - t_vert0_exp(:,1:numcombs);          % vector from vert0 to ray origin
%     pvec(:,1:numcombs)  = cross(ray_dirs_exp(:,1:numcombs), t_ev2_exp(:,1:numcombs) );  % begin calculating determinant - also used to calculate U parameter
%     det(1:numcombs)     = sum( t_ev1_exp(:,1:numcombs) .* pvec(:,1:numcombs) );   % determinant of the matrix M = dot(edge1,pvec)
% 
%     angleOK = (abs(det(1:numcombs))>eps); % if determinant is near zero then ray lies in the plane of the triangle
%   
%     if all(~angleOK)  % if all parallel than no intersections
%         %disp('tacos'); disp(' '); 
%         %return; 
%         break;
%     else 
%         det(~angleOK) = nan;              % change to avoid division by zero
%         u(1:numcombs)    = sum(tvec(:,1:numcombs).*pvec(:,1:numcombs))./det(1:numcombs);    % 1st barycentric coordinate
%         v(1:numcombs) = nan+zeros(1,numcombs); 
%         t(1:numcombs) = v(1:numcombs) ;
% %  disp(['V ' num2str(size(v))]);
%     ok = angleOK & u(1:numcombs)>=eps & u(1:numcombs)<=1.0-eps; % mask, check for first barycentric coord in tri
% %  disp(['OK ' num2str(size(ok))]);
%   % if all line/plane intersections are outside the triangle than no intersections
%     if ~any(ok), intersect = ok;  
%     
%     else
%       qvec = cross(tvec(:,ok), t_ev1_exp(:,ok)); % prepare to test V parameter
%       v(:,ok) = sum(ray_dirs_exp(:,ok).*qvec) ./ det(:,ok); % 2nd barycentric coordinate
%       t(:,ok) = sum(t_ev2_exp(:,ok).*qvec)./det(:,ok);
% 
% 
%       % test if line/plane intersection is within the triangle
%       ok = (ok & v(1:numcombs)>=eps & u(1:numcombs)+v(1:numcombs)<=1.0-eps);
% 
%       intersect = (ok & t(1:numcombs)>=eps); 
%     end
%     end
% 
%   
%     
%     %[intersect, hit_t] = TriangleRayIntersectionMultiMK( ray_origins_exp, ray_dirs_exp, t_vert0_exp(:,1:(numcombs)), t_ev1_exp(:,1:(numcombs)), t_ev2_exp(:,1:(numcombs)));
%     mytimer = mytimer + toc(startsubtic); 
% 
%     intersectArr = reshape(intersect(1:numcombs),[ numtrisgood numrays]);
%     %tArr = reshape(t,[ numtrisgood numrays]);
%     
%     
%     
%     [hit_mask_ar, hitidxs_ar] = max(intersectArr);
%     hit_ts = t( sub2ind([ numtrisgood  numrays ],  hitidxs_ar, 1:length(hitidxs_ar)   ) );
%     %hit_ts = tArr( sub2ind([ numtrisgood numrays ],   1:length(hitidxs_ar), hitidxs_ar  ) )
%     
%   
%     
%     hit_mask_ar = (hit_mask_ar > 0);  % do I need this?  I think so, max() could reteturn values of 2 or more ...
%     ray_dest_idxs = hitidxs_ar(hit_mask_ar);
%     
%    % mytimer2 = mytimer2 + toc;
%     % OK do actual optics physics stuff here..
%     %myhittvals = hit_t(intersect);
%     %rr_srct = ti;
%     %rr_dest = myhitts(1);
%     %rr_dir = t_refln(:,ti);
%     %rr_pwr = t_reflp(ti);
%     %ri = 1;
%     %while 1
%     %if t_broken(rr_dest), break, end
%     ar_costheta = dot(t_norms(:,ray_dest_idxs), ray_dirs(:,hit_mask_ar));
%     hit_mask_ar( ar_costheta < 0 ) = 0;  % also filter out back-illumianted tris from the source hit mask
%     ar_costheta = ar_costheta( ar_costheta >= 0);  % remove back-illuminated from the costheta array
%     ray_dest_idxs = hitidxs_ar(hit_mask_ar);  %remove back illuminated from the destination index array
%     hit_ts = hit_ts(hit_mask_ar);
%     numrays = length(ray_dest_idxs);
%     mr_nhits(ri) = numrays; %keep track of number or rays for this iteration
%     
%     %temp_mat1 = zeros(numrays,numrays);
%     %temp_mat2 = zeros(numrays,numrays);
%     %temp_mat3 = zeros(numrays,numrays);
%     %temp_mat4 = zeros(numrays,numrays);
%     
%     % Because some triangles receive multiple ray hits, there may be multiple repeated index values in ray_dest_idxs
%     % So we have to either construct a big array to expand our calculations, or try figuring out unique() and
%     % accumarray() to get the desired behavior.  Not sure which is faster.
%     
%     if numrays == 0
%         break;
%     end
%     
%     [ray_dest_idxs_unique, ~, idxUnique] = unique(ray_dest_idxs);
%     t_mr_nh(ray_dest_idxs_unique) = t_mr_nh(ray_dest_idxs_unique) + accumarray(idxUnique, 1)';  % keep track of number of rays reaching each triangle
%     
%     ray_pwrs = ray_pwrs(hit_mask_ar);  % select only successfull rays
%     ray_dirs = ray_dirs(:,hit_mask_ar);  % select only successful rays
%     
%     if useSpecularTMM  % now calculate optical response and assign forces/power to triangles
%         ar_theta = acos(ar_costheta);
%         ar_LULi = max(min(floor(ar_theta ./ LULAngleStep)+1,LUL.num),1);
%         
%         ar_specA = LUL.A(ar_LULi);
%         ar_specR = LUL.R(ar_LULi);
%         
%         ar_abs = ar_specA .* ray_pwrs .* t_notbroken(ray_dest_idxs);
%         ar_rpwr = ray_pwrs .* ar_specR .* t_notbroken(ray_dest_idxs);
%     else
%         ar_abs = Iabs .* ray_pwrs .* t_notbroken(ray_dest_idxs); %watt
%         ar_rpwr = ray_pwrs .* Irefl .* t_notbroken(ray_dest_idxs); % watt
%     end
%     if ~isempty(idxUnique)
%         t_abs_mr(ray_dest_idxs_unique) = t_abs_mr(ray_dest_idxs_unique) + accumarray(idxUnique, ar_abs)'; % Watt
% 
%         % calculate optical force ("optical thrust") for each triangle
%         t_oth_mr(1,ray_dest_idxs_unique) = t_oth_mr(1,ray_dest_idxs_unique) + ...
%             accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(1,ray_dest_idxs) .* ar_costheta  + ... %  this is the reflection force
%             ar_abs  ./ c0 .* ray_dirs(1,:)    )';    % and this is the absorption force
%         t_oth_mr(2,ray_dest_idxs_unique) = t_oth_mr(2,ray_dest_idxs_unique) + ...
%             accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(2,ray_dest_idxs) .* ar_costheta  + ... %  this is the reflection force
%             ar_abs  ./ c0 .* ray_dirs(2,:)    )';    % and this is the absorption force
%         t_oth_mr(3,ray_dest_idxs_unique) = t_oth_mr(3,ray_dest_idxs_unique) + ...
%             accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(3,ray_dest_idxs) .* ar_costheta  + ... %  this is the reflection force
%             ar_abs  ./ c0 .* ray_dirs(3,:)    )';    % and this is the absorption force
%     else
%         t_abs_mr(ray_dest_idxs)   = t_abs_mr(ray_dest_idxs) + ar_abs; % Wattt_oth_mr(1,ray_dest_idxs_unique) + ...
%         t_oth_mr(1,ray_dest_idxs) = t_oth_mr(1,ray_dest_idxs) + ...
%                2 .* ar_rpwr ./ c0 .* t_norms(1,ray_dest_idxs) .* ar_costheta  + ... %  this is the reflection force
%                     ar_abs  ./ c0 .* ray_dirs(1,:) ;    % and this is the absorption force
%         t_oth_mr(2,ray_dest_idxs) = t_oth_mr(2,ray_dest_idxs) + ...
%                2 .* ar_rpwr ./ c0 .* t_norms(2,ray_dest_idxs) .* ar_costheta  + ... %  this is the reflection force
%                     ar_abs  ./ c0 .* ray_dirs(2,:) ;    % and this is the absorption force
%         t_oth_mr(3,ray_dest_idxs) = t_oth_mr(3,ray_dest_idxs) + ...
%                2 .* ar_rpwr ./ c0 .* t_norms(3,ray_dest_idxs) .* ar_costheta  + ... %  this is the reflection force
%                     ar_abs  ./ c0 .* ray_dirs(3,:) ;    % and this is the absorption force
%     end
%     
%     ri = ri+1;
%     
%     
%     if numrays == 0
%         break;
%     end
%    
%     if showplots
%     np2 = 0;
%     mycolor = prettycolors{ mod(ri-1, 5)+1 }; 
%     for npl = 1:numrays
%         if  hit_mask_ar(npl)
%             np2 = np2+1;
%             plot3( ray_origins(1, npl) + hit_ts(np2).* ray_dirs(1,np2) .* [0 1], ...
%                    ray_origins(2, npl) + hit_ts(np2).* ray_dirs(2,np2) .* [0 1], ...
%                    ray_origins(3, npl) + hit_ts(np2).* ray_dirs(3,np2) .* [0 1], mycolor);
%                howmanyplots = howmanyplots + 1;
%         end
%     end
%     end
%     %disp(['Debug howmanyplots = ' num2str(howmanyplots)]);
%     if ri > MAX_REFL
%         break;
%     end
%   
%     
%     
%     ray_origins = ray_origins(:,hit_mask_ar) + hit_ts.* ray_dirs ;
%               % [t_cx(ray_dest_idxs) + hit_ts.* ray_dirs(1) ;   ...
%               %  t_cy(ray_dest_idxs) + hit_ts.* ray_dirs(2) ;   ...
%               %  t_cz(ray_dest_idxs) + hit_ts.* ray_dirs(3) ];
%     
%     % calculate new reflected direction
%     %ray_origins = t_cxyz(:,ray_dest_idxs);
%     ray_pwrs = ar_rpwr;
%     ray_dirs = ray_dirs - 2 .* ar_costheta .* t_norms(:,ray_dest_idxs);
% 
% end
% 
% % toc(startrttic)
% % disp(['Function call time: ' num2str(mytimer) ]);
% % disp(['Reshape timel time: ' num2str(mytimer2) ]);
% % disp(num2str(mr_nhits))
% % disp(num2str([ sum(t_oth_mr(1,:)) sum(t_oth_mr(2,:)) sum(t_oth_mr(3,:)) ]));
% % disp(num2str(sum(t_abs_mr)));
% 
% 
% if showplots
%     if plotcolormode == 1
%         cdata=t_mr_nh;
%     else
%         cdata=t_abs_mr;
%     end
%     trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
%     set(gca,'Colormap',bjet);
%     hcolorbar = colorbar;
%     %hcolorbar.Label.String = '#ray hits per triangle';
%     %hcolorbar.Label.String = 'Multibounce absorbed power (W)';
%     hcolorbar.Label.String = cdatalabel;
%     view([0 0]);
%     axis tight;
%     set(gca,'YLim',[0 2*radiusmm]);
%     maxoth = max(vecnorm(t_oth));
%     maxoth_mr = max(vecnorm(t_oth_mr));
%     mr_ratio = maxoth_mr / maxoth;
%     quiver3(t_cx, t_cy, t_cz, t_oth(1,:), t_oth(2,:), t_oth(3,:), 'm');
%     quiver3(t_cx, t_cy, t_cz, t_oth_mr(1,:), t_oth_mr(2,:), t_oth_mr(3,:), mr_ratio, 'b');
%     title('Vector');
% end
% 
% 
% mytimer2 = mytimer2+toc(startrttic);
% end
% db_bm_timer_vec_val = toc(startrealtic);
% 
% % disp(['Vector code: ' num2str(toc(startrealtic)) 's   ' num2str(mytimer2)  's   ' num2str(mytimer) 's']);
% % disp(num2str(mr_nhits))
% %  disp(num2str([ sum(t_oth_mr(1,:)) sum(t_oth_mr(2,:)) sum(t_oth_mr(3,:)) ]));
% %  disp(num2str(sum(t_abs_mr)));
% % % disp(['Number of lines rendered: ' num2str(howmanyplots) ] );
% % disp(' ');
% 
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %hold on;
% if 1
%     
% if showplots
%     prettycolors = { 'r' 'g' 'c' 'm' 'y' };
%     figure
%     hold on;
%     %    quiver3( t_cx, t_cy, zreliefmag*t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), zreliefmag*t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead', 'off');
%     %  if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end  %plot mechanical force
%     %  vectors (can be confusing).
%     
%     quiver3(t_cx, t_cy, t_cz, t_refln(1,:), t_refln(2,:), t_refln(3,:), 'g');
%     hold on;
%     xlabel('x (mm)');
%     ylabel('y (mm)');
%     zlabel('z (mm)');
%     axis equal;
%     howmanyplots = 0;
% end
% 
% tic;
% startrttic = tic;
% 
% %for ntest=1:numtests
% 
% mr_nhits = zeros(1,MAX_REFL);
% t_mr_na = zeros(size(t_na));
% t_mr_nh = zeros(size(t_na));
% t_oth_mr = zeros(size(t_oth));
% t_abs_mr = zeros(size(t_na));
% 
% 
%  hitts = [];
% 
% 
%  for ti = 1:length(t_na)
%      if t_broken(ti)
%          continue
%      end
%      %intersect = TriangleRayIntersection(orig, dir, vert1, vert2, vert3);
%      ar_xyz = [t_cx(ti); t_cy(ti); t_cz(ti)];
%      [intersect, t] = TriangleRayIntersectionMK( ar_xyz, t_refln(:,ti), t_vert0, t_ev1, t_ev2);
%      %  intersect = TriangleRayIntersection( t_cxyz, t_refln', t_naxyz, t_nbxyz, t_ncxyz );
%      myhitts = find(intersect);
%      if myhitts
%          myhittvals = t(intersect);
%          rr_srct = ti;
%          rr_dest = myhitts(1);
%          rr_dir = t_refln(:,ti);
%          rr_pwr = t_reflp(ti);
%          ri = 1;
%          while 1
%              if t_broken(rr_dest), break, end
%              ar_costheta = dot(t_norms(:,rr_dest), rr_dir);
%              if ar_costheta < 0, warning('Negative angle!');  break, end  %Update:  We should not exclude negative
%              %angles here, apparently???  I'm still confused.
%              mr_nhits(ri) = mr_nhits(ri) + 1;
%              t_mr_nh(ti) = t_mr_nh(ti) +1;
%              t_mr_na(rr_dest) = t_mr_na(rr_dest) +1;
%              hitts(end+1) = myhittvals(1); % for debug only
%             
%              if showplots
%                  mycolor = prettycolors{ mod(ri, 5)+1 };
%                  plot3( ar_xyz(1) + myhittvals(1).* rr_dir(1) .* [0 1], ...
%                      ar_xyz(2) + myhittvals(1).* rr_dir(2) .* [0 1], ...
%                      ar_xyz(3) + myhittvals(1).* rr_dir(3) .* [0 1], mycolor);
%                  howmanyplots = howmanyplots + 1;
%              end
%              if useSpecularTMM
%                  ar_theta = acos(ar_costheta);
%                  %if min(ar_theta) < 0
%                  %    warning('Negative theta in useSpecularTMM code');
%                  %   t_theta(t_theta < 0) = 0;
%                  %end
%                  %if max(t_theta) > pi/2
%                  %    warning('Reverse illuminated triangles in specularTMM code');
%                  %    t_theta(t_theta > pi/2) = pi/2 - t_theta(t_theta > pi/2);
%                  %end
%                  ar_LULi = max(min(floor(ar_theta ./ LULAngleStep)+1,LUL.num),1);
%                  
%                  ar_specA = LUL.A(ar_LULi);
%                  ar_specR = LUL.R(ar_LULi);
%                  
%                  ar_abs = ar_specA .* rr_pwr;
%                  ar_rpwr = rr_pwr * ar_specR;
%              else
%                  ar_abs = Iabs .* rr_pwr; %watt
%                  ar_rpwr = rr_pwr * Irefl; % watt
%              end
%              
%              t_abs_mr(rr_dest) = t_abs_mr(rr_dest) + ar_abs; % Watt
%              
%              % calculate optical force ("optical thrust") for each triangle
%              t_oth_mr(:,rr_dest) = t_oth_mr(:,rr_dest) + ...
%                  2 * ar_rpwr ./ c0 .* t_norms(:,rr_dest) .* ar_costheta  + ... %  this is the reflection force
%                  ar_abs       ./ c0 .* rr_dir ;    % and this is the absorption force
%              
%              if ri >= MAX_REFL  % end of loop termination
%                  break
%              end
%              
%              % calculate intersect point:
%              ar_xyz =   ar_xyz + myhittvals(1).* rr_dir;
%              %     [t_cx(rr_srct) + myhittvals(1).* rr_dir(1) ;   ...
%              %      t_cy(rr_srct) + myhittvals(1).* rr_dir(2) ;   ...
%              %      t_cz(rr_srct) + myhittvals(1).* rr_dir(3) ];
%              
%              % calculate new reflected light vector
%              rr_dir = rr_dir - 2 .* ar_costheta .* t_norms(:,rr_dest);
%              
%              % calculate new reflected ray intersections
%              %ar_xyz = [t_cx(myhitts(1)); t_cy(myhitts(1)); t_cz(myhitts(1))];
%              
%              [intersect, t] = TriangleRayIntersectionMK(  ar_xyz ,...
%                  rr_dir, t_vert0, t_ev1, t_ev2 );   % this version traces from centroids
%              
%              %[intersect, t] = TriangleRayIntersectionMK( ar_xyz, rr_dir', t_vert0, t_edg1, t_edg2, 'border','exclusive' );   % this version traces from point of reflection.
%              
%              
%              myhitts = find(intersect);
%              if myhitts
%                  myhittvals = t(intersect);
%                  rr_srct = rr_dest;
%                  rr_dest = myhitts(1);
%                  %rr_dir = t(refln(:,ti));
%                  rr_pwr = ar_rpwr;
%                  ri = ri+1;
%              else
%                  break;
%              end
%          end
%      end
%      
%      %fprintf('%g  ', myhitts); fprintf('     ');  fprintf('%g  ',myhittvals);  fprintf('\n');
%  end
% 
% % disp('Iteratative version');
% %toc
% db_bm_timer_itr_val = toc;
% % fprintf('Hits:  '); fprintf('%d  ', mr_nhits); fprintf('\n');
% % %disp(num2str(mr_nhits))
% % disp(num2str([ sum(t_oth_mr(1,:)) sum(t_oth_mr(2,:)) sum(t_oth_mr(3,:)) ]));
% % disp(num2str(sum(t_abs_mr)));
% % disp(['Number of lines rendered: ' num2str(howmanyplots) ] );
% 
% if showplots
% if plotcolormode == 1
%     cdata=t_mr_na;
% else
%     cdata=t_abs_mr;
% end
% trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
% set(gca,'Colormap',bjet);
% hcolorbar = colorbar;
% %hcolorbar.Label.String = '#ray hits per triangle';
% %hcolorbar.Label.String = 'Multibounce absorbed power (W)';
% hcolorbar.Label.String = cdatalabel;
% view([0 0]);
% axis tight;
% set(gca,'YLim',[0 2*radiusmm]);
% maxoth = max(vecnorm(t_oth));
% maxoth_mr = max(vecnorm(t_oth_mr));
% mr_ratio = maxoth_mr / maxoth;
% quiver3(t_cx, t_cy, t_cz, t_oth(1,:), t_oth(2,:), t_oth(3,:), 'm');
% quiver3(t_cx, t_cy, t_cz, t_oth_mr(1,:), t_oth_mr(2,:), t_oth_mr(3,:), mr_ratio, 'b');
% title('For loop');
% end
% 
% end
% % toc(startrttic)
% % disp(num2str(mr_nhits))
% %  disp(num2str([ sum(t_oth_mr(1,:)) sum(t_oth_mr(2,:)) sum(t_oth_mr(3,:)) ]));
% %  disp(num2str(sum(t_abs_mr)));
% % disp(' ');
% 
% %end
