function rgb = hex2rgb(hexstring)

for icol = 1:numel(hexstring)
    str = hexstring{icol};
    rgb{icol} = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

end