function [A] = setUpperTri2NaN(A,k)

if nargin<2
    k=-1;
end

ii=ones(size(A));
idx=tril(rot90(ii),-1);
A(~idx)=nan;

end