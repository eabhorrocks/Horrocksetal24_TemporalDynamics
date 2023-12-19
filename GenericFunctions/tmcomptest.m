function tmcomptest(X,alpha)
%TMCOMPTEST Tukey's HSD multiple comparisons test among proportions. 
% A statistical analysis used for count data is to compare the proportions 
% of a binary outcome of interest for two or more groups. 
% This m-file it is an implemented alternative to the Marascuillo's test and 
% identifies significant pairwise contrasts that this last one does not. 
% If the previously Chi-square test results in a p-value smaller to an alpha-
% value, the results are deemed significant, the null hypothesis that all
% proportions are equal is rejected, and it is concluded that there is a 
% significant difference between at least two of the proportions. 
% We can look at the pairwise comparisons obtained out of the Tukey's honest
% significant difference test empoyed. You will observe the comparisons 
% identified as significant. We can compare up to 20 proportions.
%
% Syntax: tmcomptest(X,alpha) 
%      
% Input:
%      X - data matrix size rx2. Rows=number of groups. Columnn 1=number of
%          successes or interested events. Column 2=number of trials.
%  alpha - significance level (default = 0.05)
%
% Output:
%        - Complete Tukey's HSD multiple comparisons test among proportions
%          table
%
% Example: From Zar (1999), after a significant Chi-square test of 4 group 
% proportions (P = 0.03),  we are interested to analyse which pair or pairs
% of groups have significantly different proportions with a significance of
% 0.05. Data table is as follows:
%
%                                              Group
%                                  -----------------------------
%                                    1       2       3       4  
%                                  -----------------------------
%  Successes (events of interest)   32      43      16       9
%                          Fails    55      65      64      16
%                                  -----------------------------
%                         Trials    87     108      80      25
%                                  -----------------------------
%                  Proportion (p) 0.3678  0.3981  0.2000  0.3600
%                                  -----------------------------
%    
% Data vector is:
%  X=[32 87;43 108;16 80;9 25];
%
% Calling on Matlab the function: 
%            mcomptest(X)
%
% Answer is:
%
% Group proportions of the interested events are:
% -------------------------------------------
% Group          Proportions
% -------------------------------------------
%   1               0.3678
%   2               0.3981
%   3               0.2000
%   4               0.3600
% -------------------------------------------
%
% Tukey's multiple comparisions test among proportions, k = 4
% -------------------------------------------------------------------------
%  Comparision      Difference*        SE        q        qc      Decision
% -------------------------------------------------------------------------
%     4  3           10.3574         4.603     2.250     3.633       NS
%     4  2            1.9932         4.458     0.447     3.633       NS
%     4  1            0.2396         4.559     0.053     3.633       NS
%     3  2           12.3507         2.980     4.145     3.633        S
%     3  1           10.5970         3.128     3.387     3.633       NS
%     2  1            1.7536         2.911     0.602     3.633       NS
% -------------------------------------------------------------------------
% With a given significance level of: 0.05
% The multiple comparisions can be significant (S) or not significant (NS).
% * After asin transformation of proportions.
%
% Created by A. Trujillo-Ortiz, R. Hernandez-Walls, K. Barba-Rojo and 
%            A. Castro-Perez
%            Facultad de Ciencias Marinas
%            Universidad Autonoma de Baja California
%            Apdo. Postal 453
%            Ensenada, Baja California
%            Mexico.
%            atrujo@uabc.mx
%
% Copyright. July 4, 2007.
%
% To cite this file, this would be an appropriate format:
% Trujillo-Ortiz, A., R. Hernandez-Walls, K. Barba-Rojo and A. Castro-Perez.
%   (2007). tmpomptest:Tukey's HSD multiple comparisons test among proportions.
%   A MATLAB file. [WWW document]. URL http://www.mathworks.com/matlabcentral/
%   fileexchange/loadFile.do?objectId=15499
%
% References:
%  Zar, J. H. (1999), Biostatistical Analysis (2nd ed.). 
%           NJ: Prentice-Hall, Englewood Cliffs. p.564. 
%
switch nargin
    case{2}
        if isempty(X) == false && isempty(alpha) == false
            if (alpha <= 0 || alpha >= 1)
                fprintf('Warning: Significance level error; must be 0 < alpha < 1 \n');
                return;
            end
        end
    case{1}
        alpha = 0.05;
    otherwise 
        error('Requires at least one input argument.');
end
s = X(:,1);  %number of successes
t = X(:,2);  %number of trials
is = fix(s)==s;
it = fix(t)==t;
if (sum(is) ~= length(s)) | (sum(it) ~= length(t)),
    error('Data is not in integer format. Please, check it.');
    return;
end
f = t-s;
if any(f<0)
    error('Number of events must be less than the number of trials. Check it.');
end
k = length(X);
if k < 2,
    error('Requires at least two groups.');
    return;
end
if k > 20,
    error('Cannot compare more than 20 proportions.');
    return,
end
%arcsin transform of each proportion as follows
p = [];
for i = 1:k
    x = 0.5*(asin(sqrt(s(i)/(t(i)+1)))+asin(sqrt((s(i)+1)/(t(i)+1))))*180/pi;
    p = [p x];
end
e = s;
CO=[];dpx=[];se=[];
for i = k:-1:2
    for s = i-1:-1:1
        if s ~= i
            Co = [i s];
            CO = [CO;Co];
            difp = abs(p(i)-p(s));
            dpx = [dpx;difp];
            sep = sqrt((410.35/(t(i)+0.5))+(410.35/(t(s)+0.5)));
            se = [se;sep];
        end
    end
end
q = dpx./se;
      
%Tukey's critical values up to 20 groups for infinity degrees of freedom for
%alpha-value=0.01, 0.05 and 0.10
if alpha == 0.01;
    c = [3.643 4.120 4.403 4.603 4.757 4.882 4.987 5.078 5.157 5.227 5.290 ... 
        5.348 5.400 5.448 5.493 5.535 5.574 5.611 5.645];
elseif alpha == 0.05;
    c = [2.772 3.314 3.633 3.858 4.030 4.170 4.286 4.387 4.474 4.552 4.622 ... 
        4.685 4.743 4.796 4.845 4.891 4.934 4.974 5.012];
else alpha == 0.10;
    c = [2.326 2.902 3.240 3.478 3.661 3.808 3.931 4.037 4.129 4.211 4.285 ... 
        4.351 4.412 4.468 4.519 4.568 4.612 4.654 4.694];
end
        
v = k-1;
qc = c(v);
c = size(CO);
c = c(1);
Ds = [];
for i = 1:c
    if (q(i) >= qc);
        ds = ' S';
    else (q(i) < qc);
        ds = 'NS';
    end;
    Ds = [Ds;ds];
end;
warning ('off');
qc = [ones(size(q))]*qc;
%qc = [ones(size(q),1)]*qc;

warning ('on');
ss = 1:k
ps = e./t
disp(' '); 
disp('Group proportions of the interested events are:');
disp('-------------------------------------------');
disp('Group          Proportions');
disp('-------------------------------------------');
fprintf('  %d    %17.4f\n',[ss',ps].');
disp('-------------------------------------------');
disp(' ');
fprintf('Tukey''s multiple comparisions test among proportions, k = %.i\n', k);
disp('-------------------------------------------------------------------------');
disp(' Comparision      Difference*        SE        q        qc      Decision');
disp('-------------------------------------------------------------------------');
for ml=1:k*((k-1)/2)
    fprintf('    %d  %d       %11.4f        %6.3f     %4.3f     %4.3f       %s\n', ... 
        CO(ml,1),CO(ml,2),dpx(ml),se(ml),q(ml),qc(ml,:),Ds(ml,:));
end;
disp('-------------------------------------------------------------------------');
fprintf('With a given significance level of:% 3.2f\n', alpha);
disp('The multiple comparisions can be significant (S) or not significant (NS).');
disp('* After asin transformation of proportions.');
return,