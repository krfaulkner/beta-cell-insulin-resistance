% use calculated steady state values instead of simulated

% PARAMETERS
% Topp parameters (from 2000 J. Theor. Biol. paper)
% R0=864;
% EG0=1.44; 
% SI=0.72;
% sigma=43.2;
% alpha=2*(10^4);
% k=432;
% d0=0.06;
% r1=0.84*(10^(-3));
% r2=0.24*(10^(-5));

R0=10000;
EG0=0.5; 
SI=0.1;
% b=0.1;
sigma=600;
n=1;
alpha=5*(10^4);
k=400;
d0=0.06;
r1=0.84*(10^(-3));
r2=0.24*(10^(-5));
beta=200;

% unit conversions
aG=0.0555; %mM glu/(mg/dL)
aI=0.00144; %nM ins/(mu U/mL)
R0=R0*aG;
SI=SI/aI;
sigma=sigma*aI;
alpha=alpha*(aG^2);
r1=r1/aG;
r2=r2/(aG^2);

epsList=sort([0:0.5:4 3.3605]);


% R0=430;
% EG0=1.44; 
% SI=0.9;
% sigma=43.2;
% alpha=500;
% k=70;
% d0=0.047;
% r1=0.9*(10^(-3));
% r2=0.187*(10^(-5));

figure(1)
clf
subplot(2,1,1)    
ylabel('Glucose (mM)')
xlabel('Time (hr)')
ylim([0 40])
xlim([-0.3 2])
title("Glucose Tolerance Tests for S_{\beta}="+num2str(min(epsList))+", "+num2str(min(epsList(2:(length(epsList)-1))))+",..., "+num2str(max(epsList)))
subplot(2,1,2)
ylabel('Insulin (nM)')
xlabel('Time (hr)')
ylim([0 2.5])
xlim([-0.3 2])
% subplot(3,1,3)
% ylabel('Beta cell mass (mg)')
% xlabel('Time (hr)')


figure(2)
clf 
xlabel('Glucose (mM)')
ylabel('Insulin (nM)')
% xlim([0 40])
% ylim([0 1.25])
title("Phase Plane - Glucose Tolerance Tests for S_{\beta}="+num2str(min(epsList))+", "+num2str(min(epsList(2:(length(epsList)-1))))+",..., "+num2str(max(epsList)))

N=length(epsList);
p=double.empty(N,0);
q=double.empty(N,0);
r=double.empty(N,0);
s=double.empty(N,0);

% x0 = [250;10;400];
% x0 = [160;4;12];
kx0 = [350;1.5;beta];
tspan=(0:0.01:2)./24;
tneg=(-0.3:0.01:-0.01)./24;
polyB=[-r2 +r1 -d0];
Gstar1=max(roots(polyB));
Istar1=R0./(70.*Gstar1)-EG0./70;
bstar=(k/sigma).*Istar1.*(1+(3.4.*Istar1).^n).*(alpha+(Gstar1.^2))./(Gstar1.^2);


for m=1:N
    epsilon=epsList(m);
     
    polyG=[-bstar*sigma*(SI^2)-k*EG0*SI+k*epsilon+k*epsilon*(EG0^2), k*SI*R0-2*k*epsilon*EG0*R0, k*epsilon*(R0^2)+k*epsilon*alpha*(EG0^2)-k*alpha*SI*EG0, k*alpha*SI*R0-2*k*epsilon*alpha*EG0*R0, k*epsilon*alpha*(R0^2)];
    Groots=roots(polyG);
    Groots2=[];
    for j=1:length(Groots)
        if imag(Groots(j))==0 && Groots(j)>0
            Groots2=[Groots2 Groots(j)];
        end
    end
    if isempty(Groots2)
        disp("No plausible roots for G found")
        disp(Groots)
        keyboard
    elseif length(Groots2)>1
        disp("Multiple roots for G found")
        disp(Groots)
        keyboard
    elseif length(Groots2)==1
        Gstar=Groots2(1);
    end
    Istar=(R0-EG0*Gstar)/(SI*Gstar);
    kx0 = [Gstar+20;Istar;bstar];

    kTopp = @(t,x) [R0-(EG0+SI*x(2))*x(1); (1/(1+(epsilon*x(2))).^n)*sigma*x(3)*(x(1)^2)/(alpha+(x(1)^2))-k*x(2); (-d0+r1*x(1)-r2*(x(1)^2))*x(3)];
%     kTopp = @(t,x) [R0-(EG0+SI*x(2))*x(1); (1/(1+epsilon*x(2)))*sigma*x(3)*(x(1)^2)/(alpha+(x(1)^2))-k*x(2); 0];
% %     try (epsI)^n and see how that changes amplitudes of I
%     [kts,kx]=ode45(kTopp,tspan,kx0);
%     kx0 = [kx(length(kx(:,1)),1)+20;kx(length(kx(:,1)),2);200];    
    [kts,kx]=ode45(kTopp,tspan,kx0);
%     kx0 = [kx(length(kx(:,1)),1)+20;kx(length(kx(:,1)),2);200];    
%     [kts,kx]=ode45(kTopp,tspan,kx0);
%     kx0 = [kx(length(kx(:,1)),1)+20;kx(length(kx(:,1)),2);200];    
%     [kts,kx]=ode45(kTopp,tspan,kx0);
    

    figure(1)
    subplot(2,1,1)
    hold on
    % plot(ts, x(:,1),'o','DisplayName', 'Glucose-Topp')
    color1=[0.99*(epsilon-min(epsList))/(max(epsList)-min(epsList)),0.9*(epsilon-min(epsList))/(max(epsList)-min(epsList)),1-(epsilon-min(epsList))/(max(epsList)-min(epsList))];
    color2=[0.6-0.35*(epsilon-min(epsList))/(max(epsList)-min(epsList)),0.6-0.3*(epsilon-min(epsList))/(max(epsList)-min(epsList)),0.3-0.2*(epsilon-min(epsList))/(max(epsList)-min(epsList))];
    p(m)=plot(24*[tneg'; kts], [Gstar*ones(length(tneg),1); kx(:,1)],'Linewidth', 1+2*int16(epsilon==3.3605),'Color', color2,'DisplayName', "S_{\beta}="+num2str(epsilon));
%     plot(24*tneg, Gstar*ones(length(tneg),1),'o','MarkerSize', 2,'Color', color2)
%     legend('Location', 'northOutside')
    subplot(2,1,2)
    hold on
    % plot(ts, x(:,2),'o','DisplayName', 'Insulin-Topp')
    q(m)=plot(24*[tneg'; kts], [Istar*ones(length(tneg),1); kx(:,2)],'Linewidth', 1+2*int16(epsilon==3.3605),'Color', color2,'DisplayName', "S_{\beta}="+num2str(epsilon));
% [1-(max(epsList)-epsilon)/(max(epsList)-min(epsList)),(max(epsList)-epsilon)/(max(epsList)-min(epsList)),1]
%     legend('Location', 'northOutside')
    % subplot(2,2,3)
    % hold on
    % % plot(ts, x(:,3),'o','DisplayName', 'Beta Cell-Topp')
    % plot(24*kts, kx(:,3),'Linewidth', 2,'DisplayName', 'Beta cell')
    % ylabel('Beta cell mass ')
    % xlabel('Time (hr)')
    % legend('Location', 'northOutside')
%     subplot(3,1,3)
%     hold on
%     % plot(ts, x(:,2),'o','DisplayName', 'Insulin-Topp')
%     s(m)=plot(24*[tneg'; kts], [bstar*ones(5,1); kx(:,3)],'Linewidth', 2,'Color', color2,'DisplayName', "\epsilon="+num2str(epsilon));

    
    figure(2)
    hold on
    r(m)=plot(kx(:,1),kx(:,2),'Linewidth', 2,'Color', color2,'DisplayName', "S_{\beta}="+num2str(epsilon));
   
end

figure(1)
subplot(2,1,1)
legend(p([1 N]))
legend('Location', 'NorthEast')
subplot(2,1,2)
legend(q([1 N]))
legend('Location', 'NorthEast')
% subplot(3,1,3)
% legend(s([1 N]))
% legend('Location', 'NorthEast')

figure(2)
legend(r([1 N]))
legend('Location', 'NorthEast')

set(findall(0,'-property','FontSize'),'FontSize',16)

% figure
% hold on
% N=length(x(:,1));

% for n=1:N-1
%     plot(x(n:n+1,1),x(n:n+1,2), 'Color',[exp(-((n-N)^2)/(6*N)),exp(-((n-N/2)^2)/(6*N)),exp(-((n-1)^2)/(6*N))])
% end

% rainbowplot(x(:,1)',x(:,2)')


function rainbowplot(x, y)
%------------------------------------------------------------
%  RAINBOWPLOT Colorful linear 2-D plot
%  This function plots a line colored like a rainbow. 
%  The line is defined by x versus y pairs, which is the same
%  as a regular 2-D plot.
%
%   Usage Examples,
%
%   x = 1:100; y = randn(1,100);  
%   rainbowplot(x, y);
%
%   Kun Liu 
%   Version 1.00
%   June, 2006
%------------------------------------------------------------
if size(x, 1)~=1 || size(y, 1)~=1 
    error('x and y must be one dimensional vector...');    
end
if size(x, 2) ~= size(y, 2)
    error('x and y must have the same number of elements...');
end
length = size(x, 2);
d = length:-1:1;
p = patch([x nan],[y nan], [d nan], 'EdgeColor', 'interp');
end
