% gpu timing tester thing

A = gpuArray(rand(2000));
B = gpuArray(rand(2000));

D = gpuArray(rand(2000));
E = gpuArray(rand(2000));

tic
C = A*B;
B = A*C;
A = C*A;

toc
tic
d = gather(E(1));
toc

%f = @() A * B;
%t = gputimeit(f)
