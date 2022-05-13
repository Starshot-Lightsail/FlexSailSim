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
%
%  attempt 3.2 (v18):  I am converting the main simulation loop to all-GPU calculations, so we don't need to transfer info
%                 and from the GPU any more.  Changing to use all g_ source vectors, and deleting the "vec" and "iter"
%                 versions for good.

%% optional debug setup and figure generation:
% close all hidden;  % close prior figure

% simulate_v17;  % for debug, add a "return" in the simulate code just after the light has been turned on, after trinorms
%                 % and specular optics have been calculated.  That way we can test and debug the MR-RT code on a fresh
%                 % set of valid data.
% t_ev10 = t_ev1;
% t_ev20 = t_ev2;
% t_refln0 = t_refln;
% t_vert00 = t_vert0;

% showplots = 0;  % turn on for debugging plots  NOTE:  no longer works, requires repair

% plotcolormode = 0;  % 1 == num hits per tri, other == power absorbed per tri
% %if showplots, close all, end
% if showplots
%     if ~exist('bjet'),  load('bjet.mat'); end
%     if plotcolormode == 1
%         cdatalabel = 'Number of hits';
%     else
%         cdatalabel='Multibounce absorbed power (W)';
%     end
%     howmanyplots = 0;
% end



%% benchmarker
% ?







% if 1
    %%  GPU version
    %
    % first set up debug plot window if needed:
%     if showplots
%         prettycolors = { 'r' 'g' 'c' 'm' 'y' };
%         figure
%         hold on;
%         %    quiver3( t_cx, t_cy, zreliefmag*t_cz, t_texn(1,:).*(t_tex>0), t_texn(2,:).*(t_tex>0), zreliefmag*t_texn(3,:).*(t_tex>0), 'b', 'LineWidth',0.5, 'ShowArrowHead', 'off');
%         %  if colormode==0   quiver3(n_x,n_y,n_z,n_fx,n_fy,n_fz,'k'); end  %plot mechanical force
%         %  vectors (can be confusing).
%         
%         quiver3(t_cx, t_cy, t_cz, t_refln(1,:), t_refln(2,:), t_refln(3,:), 'g');
%         hold on;
%         xlabel('x (mm)');
%         ylabel('y (mm)');
%         zlabel('z (mm)');
%         axis equal;
%     end
    
    
    
    
    %tri_broken2orig_map = 1:numtris;  % fixme, move to main sim code, retain as state vector
    
    
%     %FIXME FIXME temporary debug, let's break some triangles:
%      g_t_notbroken(:) = true;
%      anybroken = 0;
%    g_t_notbroken(1:10:end) = 0;
%    anybroken = 1;
    
    
%    tic;
    
    ray_pwrs = g_t_reflp(:, g_t_notbroken);  % ray_pwrs is now a gpu array
    
    % now initialize GPU arrays
    if anybroken
        bigvecsize = [sum(g_t_notbroken) 3];
        gpu_ev1 =        g_t_ev1(:,g_t_notbroken)';
        gpu_ev2 =        g_t_ev2(:,g_t_notbroken)';
        gpu_v0 =         g_t_vert0(:,g_t_notbroken)';
        gpu_t_norms =    g_t_norms(:,g_t_notbroken); %not used in arrayfun
        gpu_ray_dirs =   g_t_refln(:,g_t_notbroken);
        tri_broken2orig_map = g_thisseemsdumb(g_t_notbroken);
        gpu_oxyzs = [ g_t_cx(g_t_notbroken);  g_t_cy(g_t_notbroken);  g_t_cz(g_t_notbroken)];  

    else
        bigvecsize = fliplr(size(t_ev1));
        gpu_ev1 =      g_t_ev1'; % const
        gpu_ev2 =      g_t_ev2'; % const
        gpu_v0 =       g_t_vert0'; % const
        gpu_t_norms =  g_t_norms; %not used in arrayfun, const
        gpu_ray_dirs = g_t_refln; % varies
        gpu_oxyzs = [ g_t_cx;  g_t_cy;  g_t_cz];

    end
    
    g_t_mr_nh = gpuArray(zeros(size(t_na)));
    g_t_oth_mr = gpuArray(zeros(size(t_oth)));
    g_t_abs_mr = gpuArray(zeros(size(t_na)));
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
            %Result_hit_tri_idxs_nb = gather(gpu_hit_dest_tn);
            ray_dest_idxs = tri_broken2orig_map(gpu_hit_dest_tn);
            
            % check for early escape option, if all rays failed:
            if isempty(ray_dest_idxs) %if that's empty, we're done!
                break;
            end
            mr_nhits(ri) = length(ray_dest_idxs);  % record number of hit rays
            
            % still here?  I guess we should bring back the rest of the data
            gpu_remain_dirs = gpu_ray_dirs(:,gpu_m);
            %Result_costheta = gather(gpu_ar_costheta);
            %Result_theta = gather(gpu_ar_theta);
            %Result_theta = acos(Result_costheta);
            %Result_mask = gather(gpu_m);

            
%             %% Bring these back from the GPU, but only if only needed for plot/debug:
%             if showplots
%                 Result_src_xyzs = gather(gpu_oxyzs(:,gpu_m));
%                 gpu_dest_xyzs = gpu_oxyzs(:,gpu_m) + gpu_hit_t.* gpu_ray_dirs(:,gpu_m);
%                 Result_dest_xyzs = gather(gpu_dest_xyzs);
%             end
            
            %% Now we can go back and update the gpu arrays for the next iteration, if needed
            if (ri >= MAX_REFL )
                %Don't need to update GPU arrays, this is the last iteration
            else % I hope these commands might run in parallel on the GPU while we process the output on the CPU, that's why I'm initiating them here instead of at end-of-loop
                % V18 edit, I'm still kicking off this in hopes of parallel processing speedup
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
            
            %% Back at the ranch, we can do the rest of our CPU calcs   Edit:  We are now all on GPU
            % downselect ray power array to the successful rays only
            ray_pwrs = ray_pwrs(gpu_m);

            %tri_map_indxs = tri_broken2orig_map( ... % this is already Result_hit_tri_idxs
            
            % perform unique() stuff so we can map the rays back to triangles in an accumulator
            [ray_dest_idxs_unique, ~, idxUnique] = unique(ray_dest_idxs);  % are these gpuArrays now?
            g_t_mr_nh(ray_dest_idxs_unique) = g_t_mr_nh(ray_dest_idxs_unique) + accumarray(idxUnique, 1)';  % keep track of number of rays reaching each triangle
            
            % do LUL using GPU calculated values for costheta.  TODO:  compare GPU vs CPU calculation of acos?
            
            if useSpecularTMM  % now calculate optical response and assign forces/power to triangles
                ar_LULi = max(min(floor(gpu_ar_theta ./ LULAngleStep)+1,LUL.num),1);  % ar_LULi should now be a gpuArray
                
                ar_specA = g_LUL.A(ar_LULi);
                ar_specR = g_LUL.R(ar_LULi);
                
                ar_abs = ar_specA .* ray_pwrs ;
                ar_rpwr = ray_pwrs .* ar_specR ;
                
                ray_pwrs = ray_pwrs .* ar_specR;
            else
                ar_abs = Iabs .* ray_pwrs ; %watt
                ar_rpwr = ray_pwrs .* Irefl ; % watt
                
                ray_pwrs = ray_pwrs .* Irefl;
            end
            
            if ~isempty(idxUnique)
                g_t_abs_mr(ray_dest_idxs_unique) = g_t_abs_mr(ray_dest_idxs_unique) + accumarray(idxUnique, ar_abs)'; % Watt
                
                % calculate optical force ("optical thrust") for each triangle
                g_t_oth_mr(1,ray_dest_idxs_unique) = g_t_oth_mr(1,ray_dest_idxs_unique) + ...
                    accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(1,ray_dest_idxs) .* gpu_ar_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* gpu_remain_dirs(1,:)    )';    % and this is the absorption force
                g_t_oth_mr(2,ray_dest_idxs_unique) = g_t_oth_mr(2,ray_dest_idxs_unique) + ...
                    accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(2,ray_dest_idxs) .* gpu_ar_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* gpu_remain_dirs(2,:)    )';    % and this is the absorption force
                g_t_oth_mr(3,ray_dest_idxs_unique) = g_t_oth_mr(3,ray_dest_idxs_unique) + ...
                    accumarray(idxUnique, 2 * ar_rpwr ./ c0 .* t_norms(3,ray_dest_idxs) .* gpu_ar_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* gpu_remain_dirs(3,:)    )';    % and this is the absorption force
            else
                g_t_abs_mr(ray_dest_idxs)   = g_t_abs_mr(ray_dest_idxs) + ar_abs; % Wattt_oth_mr(1,ray_dest_idxs_unique) + ...
                g_t_oth_mr(1,ray_dest_idxs) = g_t_oth_mr(1,ray_dest_idxs) + ...
                    2 .* ar_rpwr ./ c0 .* t_norms(1,ray_dest_idxs) .* gpu_ar_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* gpu_remain_dirs(1,:) ;    % and this is the absorption force
                g_t_oth_mr(2,ray_dest_idxs) = g_t_oth_mr(2,ray_dest_idxs) + ...
                    2 .* ar_rpwr ./ c0 .* t_norms(2,ray_dest_idxs) .* gpu_ar_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* gpu_remain_dirs(2,:) ;    % and this is the absorption force
                g_t_oth_mr(3,ray_dest_idxs) = g_t_oth_mr(3,ray_dest_idxs) + ...
                    2 .* ar_rpwr ./ c0 .* t_norms(3,ray_dest_idxs) .* gpu_ar_costheta  + ... %  this is the reflection force
                    ar_abs  ./ c0 .* gpu_remain_dirs(3,:) ;    % and this is the absorption force
            end
     
            %Note:  I broke plot functionality by going to full GPU calcs.  The variables needed are already
            %overwritten!
%             if showplots
%                 mycolor = prettycolors{ mod(ri, 5)+1 }; 
%                 for npl = 1:length(ray_dest_idxs)
%                         plot3( [Result_src_xyzs(1,npl)  Result_dest_xyzs(1,npl)], ...
%                                [Result_src_xyzs(2,npl)  Result_dest_xyzs(2,npl)], ...
%                                [Result_src_xyzs(3,npl)  Result_dest_xyzs(3,npl)], mycolor );
%                         howmanyplots = howmanyplots + 1;
%                 end
%             end
            %disp(['Debug howmanyplots = ' num2str(howmanyplots)]);
            
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
 
 
% 
% if showplots
% if plotcolormode == 1
%     cdata=t_mr_nh;
% else
%     cdata=t_abs_mr;
% end
% trisurf(TRI(goodTs',:),n_x,n_y,n_z,  cdata(goodTs),'FaceLighting','Gouraud','AmbientStrength',0.5,'SpecularStrength',0.3,'DiffuseStrength',0.5);
%     set(gca,'Colormap',bjet);
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
% title('GPU');
% end
% 
% 
%     % end of GPU version
% 
% end
% 
% %return; % return here to skip other versions
% 
% 



