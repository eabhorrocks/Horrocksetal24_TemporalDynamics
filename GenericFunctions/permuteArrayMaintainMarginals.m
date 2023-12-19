function outputArray = permuteArrayMaintainMarginals(A,N)

% A is input array
% N is number of rearrangements to perform

[n, m] = size(A);
%N = 10000;

for k = 1:N
    r = randperm(n,2);
    r1 = r(1);
    r2 = r(2);
    c = randperm(m,2);
    c1 = c(1);
    c2 = c(2);
    d = randi([-min(A(r1,c1),A(r2,c2)),min(A(r1,c2),A(r2,c1))]);
    A([r1,r2],[c1,c2]) = A([r1,r2],[c1,c2]) + [d,-d;-d,d];
end

outputArray = A;
end