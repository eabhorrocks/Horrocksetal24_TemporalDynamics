%% code for saving images as jpeg with the actual size
% Downloaded from here: http://www.mathworks.com/matlabcentral/answers/102382-how-do-i-specify-the-output-sizes-of-jpeg-png-and-tiff-images-when-using-the-print-function-in-mat
function screen2jpeg(filename)
%SCREEN2JPEG Generate a JPEG file of the current figure with
%   dimensions consistent with the figure's screen dimensions.
%
%   SCREEN2JPEG('filename') saves the current figure to the
%   JPEG file "filename".
%
%    Sean P. McCarthy
%    Copyright (c) 1984-98 by MathWorks, Inc. All Rights Reserved
if nargin < 1
     error('Not enough input arguments!')
end
oldscreenunits = get(gcf,'Units');
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');
set(gcf,'Units','pixels');
scrpos = get(gcf,'Position');
newpos = scrpos/100;
set(gcf,'PaperUnits','inches',...
     'PaperPosition',newpos)
print('-djpeg', filename, '-r100');
drawnow
set(gcf,'Units',oldscreenunits,...
     'PaperUnits',oldpaperunits,...
     'PaperPosition',oldpaperpos)
 
end