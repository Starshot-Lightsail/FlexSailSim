% plot benchmark plot

% lol BM means poop

figure
hold on
plot(BM_NT, BM_ITR)
plot(BM_NT, BM_VEC)
plot(BM_NT, BM_GPU)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('RT compute time (s)');
xlabel('Number of triangles in mesh');
title('Benchmark for raytracing, 0.8-height sphere');
legend({ 'Iterative', 'Vectorized', 'GPU' } );

figure
hold on
plot(BM_NT2, BM_ITR2)
plot(BM_NT2, BM_VEC2)
plot(BM_NT2, BM_GPU2)
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('RT compute time (s)');
xlabel('Number of triangles in mesh');
title('Benchmark for raytracing, 0.62-height sphere');
legend({ 'Iterative', 'Vectorized', 'GPU' } );

plot(BM_NT, BM_ITR, '-')
plot(BM_NT, BM_VEC, '-')
plot(BM_NT, BM_GPU, '-')


figure
hold on
hl1 = plot(BM_NT2, BM_ITR2);
hl2 = plot(BM_NT2, BM_VEC2);
hl3 = plot(BM_NT2, BM_GPU2);
set(gca,'YScale','log')
set(gca,'XScale','log')
ylabel('RT compute time (s)');
xlabel('Number of triangles in mesh');
title('Benchmark for raytracing, shalow vs. deep sphere');


plot(BM_NT, BM_ITR, '--', 'Color', get(hl1,'Color'))
plot(BM_NT, BM_VEC, '--', 'Color', get(hl2,'Color'))
plot(BM_NT, BM_GPU, '--', 'Color', get(hl3,'Color'))

legend({ 'Iterative', 'Vectorized', 'GPU' } );