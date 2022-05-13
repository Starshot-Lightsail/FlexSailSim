            
testn_of = (zeros(3,numnodes)); % Matrix for optical forces on each node
testn_of2 = testn_of;

numtries = 10000;
tic;
for ntest=1:numtries
testn_of(1,:) = t_oth(1,:) * M_t2n ./ 3;
testn_of(2,:) = t_oth(2,:) * M_t2n ./ 3;
testn_of(3,:) = t_oth(3,:) * M_t2n ./ 3;
end
time1 = toc;
tic;
for ntest=1:numtries
testn_of2 = t_oth * M_t2n ./ 3;
end
time2 = toc;

time2/time1

% Result:  the single-line version is 4.2x slower for M_t2n size of [1536, 817]
% Result 2:  1.9x slower for [6144        3169]
% Result 3:  3-4x faster for [216   127]
%
% Conclusion:  breakeven point is around [600   331]