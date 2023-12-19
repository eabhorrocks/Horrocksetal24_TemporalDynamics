N = 360;
R = linspace(0,1000,N)./1000; % (distance in km)
Az = linspace(0,360,N); % in degrees
%[~,~,windSpeed] = peaks(N); % radial wind speefor i = 1:8
binSize = 360/8;
windSpeed = zeros(360);
windSpeed(:,360-22:360) = test(1); windSpeed(:,1:22) = test(1);
for i=2:8
windSpeed(:,23:23+(i-1)*binSize) = test(i);
end


figure
[~,c]= polarPcolor(R,Az,windSpeed);
ylabel(c,' radial wind speed (m/s)');
set(gcf,'color','w')
c.Limits = [-0.15 0.15]