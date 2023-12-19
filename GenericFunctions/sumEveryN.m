function output = sumEveryN(input, N, dim)

if N<2
    error('N must be > 1')
end

s = size(input);

% only works for 2d 
id = find(~ismember([1 2], dim));

for ii = 1:s(id)
    in_temp = input(:,ii);
    output(:,ii) = reshape(sum(reshape(in_temp,N,[]),'omitnan'),[],1);
end

end