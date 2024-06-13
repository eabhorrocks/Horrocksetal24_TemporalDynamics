%% plot example latent factors
speedcols = inferno(6);


for isesh = 1:5;
figure
allTraj = [];
for icond = 1:12
    allTraj = vertcat(allTraj, s(isesh).session.cond(icond).meanTrajectory(:,:));
end

axisLimits = [min(allTraj,[],1); max(allTraj,[],1)];



for ifactor = 1:3



     [~, idx] = max(abs(axisLimits(:,ifactor)));
         a = sign(axisLimits(idx,ifactor));
thisAxisLimit = [min(allTraj(:,ifactor),[],1)-0.3; max(allTraj(:,ifactor),[],1)+0.3];
    if a<0, thisAxisLimit = [thisAxisLimit(2)*a, thisAxisLimit(1)*a]; end
         

subplot(3,2,ifactor), hold on
title(num2str(s(isesh).session.s.propSharedVariance(ifactor)))
    for ispeed = 1:6
        plot(a*s(isesh).session.cond(ispeed).meanTrajectory(:,ifactor),'Color', speedcols(ispeed,:))
    end
ax = gca; ax.XTick = 0:20:200; ax.XTickLabel = -0.2:0.2:1.8; ax.YLim = thisAxisLimit;
    defaultAxesProperties(gca, false)
    subplot(3,2,ifactor+3), hold on
    for ispeed = 1:6
        plot(a*s(isesh).session.cond(ispeed+6).meanTrajectory(:,ifactor),'Color', speedcols(ispeed,:))
    end
ax = gca; ax.XTick = 0:20:200; ax.XTickLabel = -0.2:0.2:1.8; ax.YLim = thisAxisLimit;
    defaultAxesProperties(gca, false)
     
end


end