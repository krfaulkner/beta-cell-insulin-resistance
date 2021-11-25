

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
% SI=0.1;
b=0.1;
sigma=140000;
% epsilon=3;
alpha=2*(10^4);
k=200;
d0=6;
r1=0.84*(10^(-1));
r2=0.24*(10^(-3));

% unit conversions
aG=0.0555; %mM glu/(mg/dL)
aI=0.00144; %nM ins/(mu U/mL)
R0=R0*aG;
% SI=SI/aI;
sigma=sigma*aI;
alpha=alpha*(aG^2);
r1=r1/aG;
r2=r2/(aG^2);


SIList=100:-10:10;
epsList2=[5:-0.05:4.05,1./(0.25:.05:60), -0.05:-0.05:-0.15];
SIList2=1./(0.007:0.001:0.1);
epsList=[0 3.3605];

meanAUCgf=1.5;
semAUCgf=0.1;
meanAUCgm=1.9;
semAUCgm=0.2;

Fcol=[0.5,0.15,0.2];
Mcol=[0.1,0.2,0.6];


M=length(SIList2);
N=length(epsList2);
AUCins=double.empty(M,0);
AUCglu=double.empty(M,0);


tspan=(0:0.005:2)./24;

polyB=[-r2 +r1 -d0];
Gstar1=max(roots(polyB));
Istar1=R0./(70.*Gstar1)-EG0./70;
bstar=(k/sigma).*Istar1.*(1+3.*Istar1).*(alpha+(Gstar1.^2))./(Gstar1.^2);

for m=1:M
    SI=SIList2(m);
    for n=1:N
        epsilon=epsList2(n);
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
            Gstar=max(Groots);
        elseif length(Groots2)==1
            Gstar=Groots2(1);
        end
        Istar=(R0-EG0*Gstar)/(SI*Gstar);
        kx0 = [Gstar+20;Istar;bstar];
        kTopp = @(t,x) [R0-(EG0+SI*x(2))*x(1); (1/(1+epsilon*x(2)))*sigma*x(3)*(x(1)^2)/(alpha+(x(1)^2))-k*x(2); 0]; 
        [kts,kx]=ode45(kTopp,tspan,kx0);
        AUCins(m,n)=trapz(24*kts,kx(:,2));
        AUCglu(m,n)=trapz(24*kts,kx(:,1));
        color2=[0.8-0.4*(SI-min(SIList))/(max(SIList)-min(SIList)),0.4-0.1*(SI-min(SIList))/(max(SIList)-min(SIList)),0.4-0.1*(SI-min(SIList))/(max(SIList)-min(SIList))];
        if abs(3.4-epsilon)<0.2 && abs(70-SI)<2
            epsStar=epsilon;
            SiStar=SI;
            tstar=kts;
            xstar=kx;
            AUCstar=[AUCglu(m,n) AUCins(m,n)];

        end
    end

end

makePlots
