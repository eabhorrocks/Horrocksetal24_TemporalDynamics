function [bestParams, char, bestR2] = fitGaussianTemplates_meanPSTH(response,customOffset, plotFlag)

%% template matching gaussian fits
gaussFun =  @(params,xdata) params(1) + params(2).*exp(-(((xdata-params(3)).^2)/(2*(params(4).^2))));
opts = optimset('Display','off');
%% get x and y train vals

ytrainvals = response;
xtrainvals = 1:numel(response);

% some guesses for initial params
[maxVal, maxidx] = max(response);
[minVal, minidx] = min(response);

%% fit different template models to the response feature

%%% low pass (decay) %%%
% positive amplitude
% [base, amp, mu, sigma];
lb_low1 = [-200, 0.1, 0-customOffset*10, 0.1/10];
ub_low1 = [200, 200, customOffset, 500/10];
x0_low1 = [0, 1, -customOffset*2, 100/10];

% negative amplitude
%lb_low2 = [0, -200, 5.5, 0.7];
lb_low2 = [-200, -200, numel(response)-customOffset, 0.1/10];
ub_low2 = [200, -0.1, numel(response)+customOffset*10 500/10];
x0_low2 = [0, -1, numel(response)+customOffset*2, 100/10];

[param_out(1,:), resnorm(1)] = lsqcurvefit(gaussFun,x0_low1,xtrainvals,ytrainvals,lb_low1,ub_low1,opts);
[param_out(2,:), resnorm(2)] = lsqcurvefit(gaussFun,x0_low2,xtrainvals,ytrainvals,lb_low2,ub_low2,opts);


%%% highpass (rise) %%%
% positive amplitude
% [base, amp, mu, sigma];
%lb_high1 = [0, 0.1, 5.5, 0.7];
lb_high1 = [-200, 0.1, numel(response)-customOffset, 0.1/10];
ub_high1 = [200, 200, numel(response)+customOffset*10, 500/10];
x0_high1 = [0, 1, numel(response)+customOffset*2, 100/10];

% negative amplitude
lb_high2 = [-200, -200, -customOffset*10, 0.1/10];
%ub_high2 = [200, -0.1, 1.5, 10];
ub_high2 = [200, -0.1, customOffset, 500/10];
x0_high2 = [0, 1, -customOffset*2, 100/10];

[param_out(3,:), resnorm(3)] = lsqcurvefit(gaussFun,x0_high1,xtrainvals,ytrainvals,lb_high1,ub_high1,opts);
[param_out(4,:), resnorm(4)] = lsqcurvefit(gaussFun,x0_high2,xtrainvals,ytrainvals,lb_high2,ub_high2,opts);



%%% bandpass (peak) %%%
lb_band = [-200, 0.1, customOffset, 0.1/10];
ub_band = [200, 200, numel(response)-customOffset, 120/10];
x0_band = [0, 1, maxidx, 50/10];

[param_out(5,:), resnorm(5)] = lsqcurvefit(gaussFun,x0_band,xtrainvals,ytrainvals,lb_band,ub_band,opts);


%%% inverse (trough) %%%
%lb_inv = [0, -200, 2, 0.7];
%ub_inv = [200, -0.1, 5 10];
lb_inv = [-200, -200, customOffset, 0.1/10];
ub_inv = [200, -0.1, numel(response)-customOffset,120/10];
x0_inv = [0, 1, minidx, 50/10];

[param_out(6,:), resnorm(6)] = lsqcurvefit(gaussFun,x0_inv,xtrainvals,ytrainvals,lb_inv,ub_inv,opts);


%%% flat %%%
lb_flat = [-200, -0.1, -inf, 50/10];
ub_flat = [200, 0.1, inf, 500/10];
x0_flat = [0, 0, numel(response)/2, 100/10];

[param_out(7,:), resnorm(7)] = lsqcurvefit(gaussFun,x0_flat,xtrainvals,ytrainvals,lb_flat,ub_flat,opts);


%% determine which function meets required criteria and is best-fitting

% get quality of fits in order
[~, fitIdx] = mink(resnorm, 7);
idx_idx = 1; % default to the best fitting function
done = false;
prominence_thresh = 0.3;

while ~done
    
    bestIdx = fitIdx(idx_idx);
    bestParams = param_out(bestIdx,:);
    
    if bestIdx == 1 || bestIdx == 2
        char = 1; % decay
        
    elseif bestIdx == 3 || bestIdx == 4
        char = 2; % rise
        
    elseif bestIdx == 5
        char = 3; % peak
        [~, ~, ~, prom] = findpeaks(feval(gaussFun,[param_out(bestIdx,:)], 1:0.1:numel(response)));
        if prom < prominence_thresh % check prominence of peak meets threshold
            idx_idx = idx_idx+1;
            continue % if not try the next best function
        end
        
    elseif bestIdx == 6
        char = 4; % trough
        [~, ~, ~, prom] = findpeaks(-1*feval(gaussFun,[param_out(bestIdx,:)], 1:0.1:numel(response)));
        if prom < prominence_thresh % check prominence of dip meets threshold
            idx_idx = idx_idx+1;
            continue % if not try the next best function
        end
        
    elseif bestIdx == 7
        char = 5; % flat
    end
    
    if range( feval(gaussFun,[param_out(bestIdx,:)], 1:0.1:numel(response))) < 0.2 % if no large moduation, set as flat.
        char = 5;
    end
    
    done = true; % finish classification
end

% calculate R2 for the best fit
SS_tot = sum((ytrainvals - mean(ytrainvals)).^2);
bestR2 = 1-(resnorm(bestIdx)/SS_tot);


% optional plotting
if plotFlag
    hold off
    plot(response)
    hold on
    plot(1:0.1:numel(response), feval(gaussFun,[param_out(bestIdx,:)], 1:0.1:numel(response)), 'r-')
    title([num2str(char),', ', num2str(bestR2)])
    ylim([0 1])
    pause
    
end