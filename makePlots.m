% for plotting solutions to this system

limsx=[-0.1 4];
limsy=[10 120];

meanEpsilonFe = 3.3605;
semEpsFe = 0.1893;
meanEpsilonMa = -0.0476;
semEpsMa = 0.1243;
    
figure(1)
clf
s=surf(epsList2(1:1200), SIList2, AUCglu(:,1:1200),'FaceAlpha',1);
colormap('parula')
s.EdgeColor = 'none';
xlim([0 4]) 
ylim([10 150]) 
xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
zlabel("AUC glucose (mM*min)")
view(140,30)

figure(2)
clf
s=surf(epsList2(1:1200), SIList2, AUCins(:,1:1200),'FaceAlpha',1);
colormap('parula')
s.EdgeColor = 'none';
xlim([0 4]) 
ylim([10 150]) 
xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
zlabel("AUC insulin (nM*min)")
view(140,30)



figure(3)
clf
weeks=[4,9,21,39];
malewt=[24.6,25,18.6,17.8];
malewt_up=[43.4, 31.3, 26.2,24.8];
malewt_down=[18.5,20.9,14.3,13.8];
femalewt=[100.3,101.5,90.5,65.5];
femalewt_up=[139,125,135.5,98.5];
femalewt_down= [78.6,85.8,68,48.5];
malemut=[22.1, 23.7, 18.9, 18.2];
malemut_up=[55.6, 32.5, 26.3, 26.4];
malemut_down=[13.5, 18.6, 14.6, 13.8];
femalemut=[38.8, 38.6, 33.8, 26.6];
femalemut_up=[54.2, 55.6, 54.4, 41.2];
femalemut_down= [26.5, 29.4, 24.4, 19.6];

malemut_hfd=18.1;
malemut_hfd_up=23.5;
malemut_hfd_down=14.60;
femalemut_hfd=25.5;
femalemut_hfd_up=33.25;
femalemut_hfd_down=20.5;

hold on
Fcolw=0.7.*Fcol;
Mcolw=0.7.*Mcol;
Fcolm=1.5.*Fcol;
Mcolm=1.5.*Mcol;
plot(weeks,femalewt,':','Color',Fcolw,'LineWidth', 2)
plot(weeks,malewt,':','Color',Mcolw,'LineWidth', 2)
plot(weeks,femalemut,':','Color',Fcolm,'LineWidth', 2)
plot(weeks,malemut,':','Color',Mcolm,'LineWidth', 2)

for i=1:length(weeks)
    q(1)=plot(weeks(i)*ones(1,2), [femalewt_up(i), femalewt_down(i)],'Color',Fcolw,'LineWidth', 2, 'DisplayName', 'Female wild type');
    plot(weeks(i),femalewt(i),'o','MarkerEdgeColor',Fcolw,'MarkerFaceColor',Fcolw,'MarkerSize', 6)
    plot(weeks(i)*ones(1,3)+linspace(-1,1,3),(femalewt_up(i))*ones(1,3),'LineWidth', 2,'Color',Fcolw,'MarkerSize', 10)
    plot(weeks(i)*ones(1,3)+linspace(-1,1,3),(femalewt_down(i))*ones(1,3),'LineWidth', 2,'Color',Fcolw,'MarkerSize', 10)
    q(2)=plot(weeks(i)*ones(1,2), [malewt_up(i), malewt_down(i)],'Color',Mcolw,'LineWidth', 2, 'DisplayName', 'Male wild type');
    plot(weeks(i),malewt(i),'o','MarkerEdgeColor',Mcolw,'MarkerFaceColor',Mcolw,'MarkerSize', 6)
    plot(weeks(i)*ones(1,3)+linspace(-1,1,3),(malewt_up(i))*ones(1,3),'LineWidth', 2,'Color',Mcolw,'MarkerSize', 10)
    plot(weeks(i)*ones(1,3)+linspace(-1,1,3),(malewt_down(i))*ones(1,3),'LineWidth', 2,'Color',Mcolw,'MarkerSize', 10)
    q(3)=plot(weeks(i)*ones(1,2), [femalemut_up(i), femalemut_down(i)],'Color',Fcolm,'LineWidth', 2, 'DisplayName', 'Female knockout');
    plot(weeks(i),femalemut(i),'o','MarkerEdgeColor',Fcolm,'MarkerFaceColor',Fcolm,'MarkerSize', 6)
    plot(weeks(i)*ones(1,3)+linspace(-1,1,3),(femalemut_up(i))*ones(1,3),'LineWidth', 2,'Color',Fcolm,'MarkerSize', 10)
    plot(weeks(i)*ones(1,3)+linspace(-1,1,3),(femalemut_down(i))*ones(1,3),'LineWidth', 2,'Color',Fcolm,'MarkerSize', 10)
    q(4)=plot(weeks(i)*ones(1,2), [malemut_up(i), malemut_down(i)],'Color',Mcolm,'LineWidth', 2, 'DisplayName', 'Male knockout');
    plot(weeks(i),malemut(i),'o','MarkerEdgeColor',Mcolm,'MarkerFaceColor',Mcolm,'MarkerSize', 6)
    plot(weeks(i)*ones(1,3)+linspace(-1,1,3),(malemut_up(i))*ones(1,3),'LineWidth', 2,'Color',Mcolm,'MarkerSize', 10)
    plot(weeks(i)*ones(1,3)+linspace(-1,1,3),(malemut_down(i))*ones(1,3),'LineWidth', 2,'Color',Mcolm,'MarkerSize', 10)
end
xlim([0 45])
xticks([4 9 21 39])
xlabel("Age (weeks)")
ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
legend(q([1:4]))
legend('Location', 'NorthEast')


% --- LFD Plots --- %

figure(4)
clf
[a0,b0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd4-SEM_AUC_fe_wt_lfd4], 'LineColor', Fcol,'LineStyle',':');
[a2,b2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd4+SEM_AUC_fe_wt_lfd4], 'LineColor', Fcol,'LineStyle',':');
[c0,d0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd4-SEM_AUC_ma_wt_lfd4], 'LineColor', Mcol,'LineStyle',':');
[c2,d2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd4+SEM_AUC_ma_wt_lfd4], 'LineColor', Mcol,'LineStyle',':');
contour(epsList2, SIList2, AUCglu, 'LineWidth', 2, 'ShowText', 'on')
colormap([linspace(0,0.8,250)' linspace(0,0.8,250)' linspace(0,0.8,250)'])
hold on
[a1,b1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd4], 'LineColor', Fcol,'LineStyle',':', 'ShowText', 'on');
f1=fill([a0(1,2:length(a0(1,:))-100),fliplr(a2(1,2:length(a2(1,:))))],[a0(2,2:length(a0(1,:))-100),fliplr(a2(2,2:length(a2(1,:))))],Fcol,'EdgeColor','none');
f1.FaceAlpha=0.25;
[c1,d1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd4], 'LineColor', Mcol,'LineStyle',':', 'ShowText', 'on');
f2=fill([c0(1,2:length(c0(1,:))-50),fliplr(c2(1,2:length(c2(1,:))))],[c0(2,2:length(c0(1,:))-50),fliplr(c2(2,2:length(c2(1,:))))],Mcol,'EdgeColor','none');
f2.FaceAlpha=0.25;
% plot(meanEpsilonFe*ones(1,length(SIList2)), SIList2,'Color',Fcol,'LineWidth', 2,'LineStyle',':')
% f3=fill([(meanEpsilonFe-semEpsFe)*ones(1,length(SIList2)),(meanEpsilonFe+semEpsFe)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Fcol,'EdgeColor','none');
% f3.FaceAlpha=0.1;
% plot(meanEpsilonMa*ones(1,length(SIList2)), SIList2,'Color',Mcol,'LineWidth', 2,'LineStyle',':')
% f4=fill([(meanEpsilonMa-semEpsMa)*ones(1,length(SIList2)),(meanEpsilonMa+semEpsMa)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Mcol,'EdgeColor','none');
% f4.FaceAlpha=0.1;
plot([0, 0], [femalemut_up(1), femalemut_down(1)],'Color',Fcol,'LineWidth', 2)
plot(0,femalemut(1),'o','MarkerEdgeColor',Fcol,'MarkerFaceColor',Fcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(femalemut_up(1))*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(femalemut_down(1))*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot([0, 0], [malemut_up(1), malemut_down(1)],'Color',Mcol,'LineWidth', 2)
plot(0,malemut(1),'o','MarkerEdgeColor',Mcol,'MarkerFaceColor',Mcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(malemut_up(1))*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(malemut_down(1))*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
xlim(limsx) 
ylim(limsy) 
xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
title("LFD, 4 weeks: AUC glucose (mM*min)")


figure(5)
clf
[a0,b0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd9-SEM_AUC_fe_wt_lfd9], 'LineColor', Fcol,'LineStyle',':');
[a2,b2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd9+SEM_AUC_fe_wt_lfd9], 'LineColor', Fcol,'LineStyle',':');
[c0,d0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd9-SEM_AUC_ma_wt_lfd9], 'LineColor', Mcol,'LineStyle',':');
[c2,d2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd9+SEM_AUC_ma_wt_lfd9], 'LineColor', Mcol,'LineStyle',':');
contour(epsList2, SIList2, AUCglu, 'LineWidth', 2, 'ShowText', 'on')
colormap([linspace(0,0.8,250)' linspace(0,0.8,250)' linspace(0,0.8,250)'])
hold on
[a1,b1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd9], 'LineColor', Fcol,'LineStyle',':', 'ShowText', 'on');
f1=fill([a0(1,2:length(a0(1,:))-50),fliplr(a2(1,2:length(a2(1,:))-50))],[a0(2,2:length(a0(1,:))-50),fliplr(a2(2,2:length(a2(1,:))-50))],Fcol,'EdgeColor','none');
f1.FaceAlpha=0.25;
[c1,d1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd9], 'LineColor', Mcol,'LineStyle',':', 'ShowText', 'on');
f2=fill([c0(1,2:length(c0(1,:))-50),fliplr(c2(1,2:length(c2(1,:))))],[c0(2,2:length(c0(1,:))-50),fliplr(c2(2,2:length(c2(1,:))))],Mcol,'EdgeColor','none');
f2.FaceAlpha=0.25;
plot(meanEpsilonFe*ones(1,length(SIList2)), SIList2,'Color',Fcol,'LineWidth', 2,'LineStyle',':')
f3=fill([(meanEpsilonFe-semEpsFe)*ones(1,length(SIList2)),(meanEpsilonFe+semEpsFe)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Fcol,'EdgeColor','none');
f3.FaceAlpha=0.25;
plot(meanEpsilonMa*ones(1,length(SIList2)), SIList2,'Color',Mcol,'LineWidth', 2,'LineStyle',':')
f4=fill([(meanEpsilonMa-semEpsMa)*ones(1,length(SIList2)),(meanEpsilonMa+semEpsMa)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Mcol,'EdgeColor','none');
f4.FaceAlpha=0.25;
plot([0, 0], [femalemut_up(2), femalemut_down(2)],'Color',Fcol,'LineWidth', 2)
plot(0,femalemut(2),'o','MarkerEdgeColor',Fcol,'MarkerFaceColor',Fcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(femalemut_up(2))*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(femalemut_down(2))*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot([0, 0], [malemut_up(2), malemut_down(2)],'Color',Mcol,'LineWidth', 2)
plot(0,malemut(2),'o','MarkerEdgeColor',Mcol,'MarkerFaceColor',Mcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(malemut_up(2))*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(malemut_down(2))*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
xlim(limsx) 
ylim(limsy) 
xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
title("LFD, 9 weeks: AUC glucose (mM*min)")

figure(6)
clf
[a0,b0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd21-SEM_AUC_fe_wt_lfd21], 'LineColor', Fcol,'LineStyle',':');
[a2,b2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd21+SEM_AUC_fe_wt_lfd21], 'LineColor', Fcol,'LineStyle',':');
[c0,d0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd21-SEM_AUC_ma_wt_lfd21], 'LineColor', Mcol,'LineStyle',':');
[c2,d2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd21+SEM_AUC_ma_wt_lfd21], 'LineColor', Mcol,'LineStyle',':');
contour(epsList2, SIList2, AUCglu, 'LineWidth', 2, 'ShowText', 'on')
colormap([linspace(0,0.8,250)' linspace(0,0.8,250)' linspace(0,0.8,250)'])
hold on
[a1,b1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd21], 'LineColor', Fcol,'LineStyle',':', 'ShowText', 'on');
f1=fill([a0(1,2:length(a0(1,:))-50),fliplr(a2(1,2:length(a2(1,:))))],[a0(2,2:length(a0(1,:))-50),fliplr(a2(2,2:length(a2(1,:))))],Fcol,'EdgeColor','none');
f1.FaceAlpha=0.25;
[c1,d1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd21], 'LineColor', Mcol,'LineStyle',':', 'ShowText', 'on');
f2=fill([c0(1,2:length(c0(1,:))),fliplr(c2(1,2:length(c2(1,:))))],[c0(2,2:length(c0(1,:))),fliplr(c2(2,2:length(c2(1,:))))],Mcol,'EdgeColor','none');
f2.FaceAlpha=0.25;
% plot(meanEpsilonFe*ones(1,length(SIList2)), SIList2,'Color',Fcol,'LineWidth', 2,'LineStyle',':')
% f3=fill([(meanEpsilonFe-semEpsFe)*ones(1,length(SIList2)),(meanEpsilonFe+semEpsFe)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Fcol,'EdgeColor','none');
% f3.FaceAlpha=0.1;
% plot(meanEpsilonMa*ones(1,length(SIList2)), SIList2,'Color',Mcol,'LineWidth', 2,'LineStyle',':')
% f4=fill([(meanEpsilonMa-semEpsMa)*ones(1,length(SIList2)),(meanEpsilonMa+semEpsMa)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Mcol,'EdgeColor','none');
% f4.FaceAlpha=0.1;
plot([0, 0], [femalemut_up(3), femalemut_down(3)],'Color',Fcol,'LineWidth', 2)
plot(0,femalemut(3),'o','MarkerEdgeColor',Fcol,'MarkerFaceColor',Fcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(femalemut_up(3))*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(femalemut_down(3))*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot([0, 0], [malemut_up(3), malemut_down(3)],'Color',Mcol,'LineWidth', 2)
plot(0,malemut(3),'o','MarkerEdgeColor',Mcol,'MarkerFaceColor',Mcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(malemut_up(3))*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(malemut_down(3))*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
xlim(limsx) 
ylim(limsy) 
xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
title("LFD, 21 weeks: AUC glucose (mM*min)")


figure(7)
clf
[a0,b0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd39-SEM_AUC_fe_wt_lfd39], 'LineColor', Fcol,'LineStyle',':');
[a2,b2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd39+SEM_AUC_fe_wt_lfd39], 'LineColor', Fcol,'LineStyle',':');
[c0,d0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd39-SEM_AUC_ma_wt_lfd39], 'LineColor', Mcol,'LineStyle',':');
[c2,d2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd39+SEM_AUC_ma_wt_lfd39], 'LineColor', Mcol,'LineStyle',':');
contour(epsList2, SIList2, AUCglu, 'LineWidth', 2, 'ShowText', 'on')
colormap([linspace(0,0.8,250)' linspace(0,0.8,250)' linspace(0,0.8,250)'])
hold on
[a1,b1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_lfd39], 'LineColor', Fcol,'LineStyle',':', 'ShowText', 'on');
f1=fill([a0(1,2:length(a0(1,:))-50),fliplr(a2(1,2:length(a2(1,:))))],[a0(2,2:length(a0(1,:))-50),fliplr(a2(2,2:length(a2(1,:))))],Fcol,'EdgeColor','none');
f1.FaceAlpha=0.25;
[c1,d1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_lfd39], 'LineColor', Mcol,'LineStyle',':', 'ShowText', 'on');
f2=fill([c0(1,2:length(c0(1,:))-50),fliplr(c2(1,2:length(c2(1,:))))],[c0(2,2:length(c0(1,:))-50),fliplr(c2(2,2:length(c2(1,:))))],Mcol,'EdgeColor','none');
f2.FaceAlpha=0.25;
% plot(meanEpsilonFe*ones(1,length(SIList2)), SIList2,'Color',Fcol,'LineWidth', 2,'LineStyle',':')
% f3=fill([(meanEpsilonFe-semEpsFe)*ones(1,length(SIList2)),(meanEpsilonFe+semEpsFe)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Fcol,'EdgeColor','none');
% f3.FaceAlpha=0.1;
% plot(meanEpsilonMa*ones(1,length(SIList2)), SIList2,'Color',Mcol,'LineWidth', 2,'LineStyle',':')
% f4=fill([(meanEpsilonMa-semEpsMa)*ones(1,length(SIList2)),(meanEpsilonMa+semEpsMa)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Mcol,'EdgeColor','none');
% f4.FaceAlpha=0.1;
plot([0, 0], [femalemut_up(4), femalemut_down(4)],'Color',Fcol,'LineWidth', 2)
plot(0,femalemut(4),'o','MarkerEdgeColor',Fcol,'MarkerFaceColor',Fcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(femalemut_up(4))*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(femalemut_down(4))*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot([0, 0], [malemut_up(4), malemut_down(4)],'Color',Mcol,'LineWidth', 2)
plot(0,malemut(4),'o','MarkerEdgeColor',Mcol,'MarkerFaceColor',Mcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(malemut_up(4))*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(malemut_down(4))*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
xlim(limsx) 
ylim(limsy) 
xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
title("LFD, 39 weeks: AUC glucose (mM*min)")


% --- HFD Plots --- %

figure(10)
clf
[a0,b0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd9-SEM_AUC_fe_wt_hfd9], 'LineColor', Fcol,'LineStyle',':');
[a2,b2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd9+SEM_AUC_fe_wt_hfd9], 'LineColor', Fcol,'LineStyle',':');
[c0,d0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd9-SEM_AUC_ma_wt_hfd9], 'LineColor', Mcol,'LineStyle',':');
[c2,d2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd9+SEM_AUC_ma_wt_hfd9], 'LineColor', Mcol,'LineStyle',':');
contour(epsList2, SIList2, AUCglu, 'LineWidth', 2, 'ShowText', 'on')
colormap([linspace(0,0.8,250)' linspace(0,0.8,250)' linspace(0,0.8,250)'])
hold on
[a1,b1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd9], 'LineColor', Fcol,'LineStyle',':', 'ShowText', 'on');
f1=fill([a0(1,2:length(a0(1,:))-50),fliplr(a2(1,2:length(a2(1,:))))],[a0(2,2:length(a0(1,:))-50),fliplr(a2(2,2:length(a2(1,:))))],Fcol,'EdgeColor','none');
f1.FaceAlpha=0.25;
[c1,d1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd9], 'LineColor', Mcol,'LineStyle',':', 'ShowText', 'on');
f2=fill([c0(1,2:length(c0(1,:))),fliplr(c2(1,2:length(c2(1,:))))],[c0(2,2:length(c0(1,:))),fliplr(c2(2,2:length(c2(1,:))))],Mcol,'EdgeColor','none');
f2.FaceAlpha=0.25;
plot(meanEpsilonFe*ones(1,length(SIList2)), SIList2,'Color',Fcol,'LineWidth', 2,'LineStyle',':')
f3=fill([(meanEpsilonFe-semEpsFe)*ones(1,length(SIList2)),(meanEpsilonFe+semEpsFe)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Fcol,'EdgeColor','none');
f3.FaceAlpha=0.1;
plot(meanEpsilonMa*ones(1,length(SIList2)), SIList2,'Color',Mcol,'LineWidth', 2,'LineStyle',':')
f4=fill([(meanEpsilonMa-semEpsMa)*ones(1,length(SIList2)),(meanEpsilonMa+semEpsMa)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Mcol,'EdgeColor','none');
f4.FaceAlpha=0.1;
plot([0, 0], [femalemut_hfd_up, femalemut_hfd_down],'Color',Fcol,'LineWidth', 2)
plot(0,femalemut_hfd,'o','MarkerEdgeColor',Fcol,'MarkerFaceColor',Fcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(femalemut_hfd_up)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(femalemut_hfd_down)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
plot([0, 0], [malemut_hfd_up, malemut_hfd_down],'Color',Mcol,'LineWidth', 2)
plot(0,malemut_hfd,'o','MarkerEdgeColor',Mcol,'MarkerFaceColor',Mcol,'MarkerSize', 6)
plot(linspace(-0.05,0.05,3),(malemut_hfd_up)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
plot(linspace(-0.05,0.05,3),(malemut_hfd_down)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
xlim(limsx) 
ylim(limsy) 
yticks([20:20:120])
xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
title("HFD, 9 weeks: AUC glucose (mM*min)")
% 
% figure(2)
% clf
% [a0,b0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd21-SEM_AUC_fe_wt_hfd21], 'LineColor', Fcol,'LineStyle',':');
% [a2,b2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd21+SEM_AUC_fe_wt_hfd21], 'LineColor', Fcol,'LineStyle',':');
% [c0,d0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd21-SEM_AUC_ma_wt_hfd21], 'LineColor', Mcol,'LineStyle',':');
% [c2,d2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd21+SEM_AUC_ma_wt_hfd21], 'LineColor', Mcol,'LineStyle',':');
% contour(epsList2, SIList2, AUCglu, 'LineWidth', 2, 'ShowText', 'on')
% colormap([linspace(0,0.8,250)' linspace(0,0.8,250)' linspace(0,0.8,250)'])
% hold on
% [a1,b1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd21], 'LineColor', Fcol,'LineStyle',':', 'ShowText', 'on');
% f1=fill([a0(1,2:length(a0(1,:))-50),fliplr(a2(1,2:length(a2(1,:))))],[a0(2,2:length(a0(1,:))-50),fliplr(a2(2,2:length(a2(1,:))))],Fcol,'EdgeColor','none');
% f1.FaceAlpha=0.25;
% [c1,d1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd21], 'LineColor', Mcol,'LineStyle',':', 'ShowText', 'on');
% f2=fill([c0(1,2:length(c0(1,:))),fliplr(c2(1,2:length(c2(1,:))))],[c0(2,2:length(c0(1,:))),fliplr(c2(2,2:length(c2(1,:))))],Mcol,'EdgeColor','none');
% f2.FaceAlpha=0.25;
% plot(meanEpsilonFe*ones(1,length(SIList2)), SIList2,'Color',Fcol,'LineWidth', 2,'LineStyle',':')
% f3=fill([(meanEpsilonFe-semEpsFe)*ones(1,length(SIList2)),(meanEpsilonFe+semEpsFe)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Fcol,'EdgeColor','none');
% f3.FaceAlpha=0.1;
% plot(meanEpsilonMa*ones(1,length(SIList2)), SIList2,'Color',Mcol,'LineWidth', 2,'LineStyle',':')
% f4=fill([(meanEpsilonMa-semEpsMa)*ones(1,length(SIList2)),(meanEpsilonMa+semEpsMa)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Mcol,'EdgeColor','none');
% f4.FaceAlpha=0.1;
% plot(zeros(1,3), mean_AUC_fe_ko_hfd21*ones(1,3)+SEM_AUC_fe_ko_hfd21*[-1 0 1],'Color',Fcol,'LineWidth', 2)
% plot(0,mean_AUC_fe_ko_hfd21,'o','MarkerEdgeColor',Fcol,'MarkerFaceColor',Fcol,'MarkerSize', 6)
% plot(linspace(-0.05,0.05,3),(mean_AUC_fe_ko_hfd21+SEM_AUC_fe_ko_hfd21)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
% plot(linspace(-0.05,0.05,3),(mean_AUC_fe_ko_hfd21-SEM_AUC_fe_ko_hfd21)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
% plot(zeros(1,3), mean_AUC_ma_ko_hfd21*ones(1,3)+SEM_AUC_ma_ko_hfd21*[-1 0 1],'Color',Mcol,'LineWidth', 2)
% plot(0,mean_AUC_ma_ko_hfd21,'o','MarkerEdgeColor',Mcol,'MarkerFaceColor',Mcol,'MarkerSize', 6)
% plot(linspace(-0.05,0.05,3),(mean_AUC_ma_ko_hfd21+SEM_AUC_ma_ko_hfd21)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
% plot(linspace(-0.05,0.05,3),(mean_AUC_ma_ko_hfd21-SEM_AUC_ma_ko_hfd21)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
% xlim(limsx) 
% ylim(limsy) 
% xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
% ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
% title("HFD, 21 weeks: AUC glucose (mM*min)")
% 
% figure(3)
% clf
% [a0,b0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd25-SEM_AUC_fe_wt_hfd25], 'LineColor', Fcol,'LineStyle',':');
% [a2,b2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd25+SEM_AUC_fe_wt_hfd25], 'LineColor', Fcol,'LineStyle',':');
% [c0,d0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd25-SEM_AUC_ma_wt_hfd25], 'LineColor', Mcol,'LineStyle',':');
% [c2,d2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd25+SEM_AUC_ma_wt_hfd25], 'LineColor', Mcol,'LineStyle',':');
% contour(epsList2, SIList2, AUCglu, 'LineWidth', 2, 'ShowText', 'on')
% colormap([linspace(0,0.8,250)' linspace(0,0.8,250)' linspace(0,0.8,250)'])
% hold on
% [a1,b1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd25], 'LineColor', Fcol,'LineStyle',':', 'ShowText', 'on');
% f1=fill([a0(1,2:length(a0(1,:))-50),fliplr(a2(1,2:length(a2(1,:))))],[a0(2,2:length(a0(1,:))-50),fliplr(a2(2,2:length(a2(1,:))))],Fcol,'EdgeColor','none');
% f1.FaceAlpha=0.25;
% [c1,d1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd25], 'LineColor', Mcol,'LineStyle',':', 'ShowText', 'on');
% f2=fill([c0(1,2:length(c0(1,:))),fliplr(c2(1,2:length(c2(1,:))))],[c0(2,2:length(c0(1,:))),fliplr(c2(2,2:length(c2(1,:))))],Mcol,'EdgeColor','none');
% f2.FaceAlpha=0.25;
% plot(meanEpsilonFe*ones(1,length(SIList2)), SIList2,'Color',Fcol,'LineWidth', 2,'LineStyle',':')
% f3=fill([(meanEpsilonFe-semEpsFe)*ones(1,length(SIList2)),(meanEpsilonFe+semEpsFe)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Fcol,'EdgeColor','none');
% f3.FaceAlpha=0.1;
% plot(meanEpsilonMa*ones(1,length(SIList2)), SIList2,'Color',Mcol,'LineWidth', 2,'LineStyle',':')
% f4=fill([(meanEpsilonMa-semEpsMa)*ones(1,length(SIList2)),(meanEpsilonMa+semEpsMa)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Mcol,'EdgeColor','none');
% f4.FaceAlpha=0.1;
% plot(zeros(1,3), mean_AUC_fe_ko_hfd25*ones(1,3)+SEM_AUC_fe_ko_hfd25*[-1 0 1],'Color',Fcol,'LineWidth', 2)
% plot(0,mean_AUC_fe_ko_hfd25,'o','MarkerEdgeColor',Fcol,'MarkerFaceColor',Fcol,'MarkerSize', 6)
% plot(linspace(-0.05,0.05,3),(mean_AUC_fe_ko_hfd25+SEM_AUC_fe_ko_hfd25)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
% plot(linspace(-0.05,0.05,3),(mean_AUC_fe_ko_hfd25-SEM_AUC_fe_ko_hfd25)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
% plot(zeros(1,3), mean_AUC_ma_ko_hfd25*ones(1,3)+SEM_AUC_ma_ko_hfd25*[-1 0 1],'Color',Mcol,'LineWidth', 2)
% plot(0,mean_AUC_ma_ko_hfd25,'o','MarkerEdgeColor',Mcol,'MarkerFaceColor',Mcol,'MarkerSize', 6)
% plot(linspace(-0.05,0.05,3),(mean_AUC_ma_ko_hfd25+SEM_AUC_ma_ko_hfd25)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
% plot(linspace(-0.05,0.05,3),(mean_AUC_ma_ko_hfd25-SEM_AUC_ma_ko_hfd25)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
% xlim(limsx) 
% ylim(limsy) 
% xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
% ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
% title("HFD, 25 weeks: AUC glucose (mM*min)")
% 
% figure(4)
% clf
% [a0,b0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd39-SEM_AUC_fe_wt_hfd39], 'LineColor', Fcol,'LineStyle',':');
% [a2,b2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd39+SEM_AUC_fe_wt_hfd39], 'LineColor', Fcol,'LineStyle',':');
% [c0,d0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd39-SEM_AUC_ma_wt_hfd39], 'LineColor', Mcol,'LineStyle',':');
% [c2,d2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd39+SEM_AUC_ma_wt_hfd39], 'LineColor', Mcol,'LineStyle',':');
% contour(epsList2, SIList2, AUCglu, 'LineWidth', 2, 'ShowText', 'on')
% colormap([linspace(0,0.8,250)' linspace(0,0.8,250)' linspace(0,0.8,250)'])
% hold on
% [a1,b1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd39], 'LineColor', Fcol,'LineStyle',':', 'ShowText', 'on');
% f1=fill([a0(1,2:length(a0(1,:))-50),fliplr(a2(1,2:length(a2(1,:))))],[a0(2,2:length(a0(1,:))-50),fliplr(a2(2,2:length(a2(1,:))))],Fcol,'EdgeColor','none');
% f1.FaceAlpha=0.25;
% [c1,d1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd39], 'LineColor', Mcol,'LineStyle',':', 'ShowText', 'on');
% f2=fill([c0(1,2:length(c0(1,:))-50),fliplr(c2(1,2:length(c2(1,:))))],[c0(2,2:length(c0(1,:))-50),fliplr(c2(2,2:length(c2(1,:))))],Mcol,'EdgeColor','none');
% f2.FaceAlpha=0.25;
% plot(meanEpsilonFe*ones(1,length(SIList2)), SIList2,'Color',Fcol,'LineWidth', 2,'LineStyle',':')
% f3=fill([(meanEpsilonFe-semEpsFe)*ones(1,length(SIList2)),(meanEpsilonFe+semEpsFe)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Fcol,'EdgeColor','none');
% f3.FaceAlpha=0.1;
% plot(meanEpsilonMa*ones(1,length(SIList2)), SIList2,'Color',Mcol,'LineWidth', 2,'LineStyle',':')
% f4=fill([(meanEpsilonMa-semEpsMa)*ones(1,length(SIList2)),(meanEpsilonMa+semEpsMa)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Mcol,'EdgeColor','none');
% f4.FaceAlpha=0.1;
% plot(zeros(1,3), mean_AUC_fe_ko_hfd39*ones(1,3)+SEM_AUC_fe_ko_hfd39*[-1 0 1],'Color',Fcol,'LineWidth', 2)
% plot(0,mean_AUC_fe_ko_hfd39,'o','MarkerEdgeColor',Fcol,'MarkerFaceColor',Fcol,'MarkerSize', 6)
% plot(linspace(-0.05,0.05,3),(mean_AUC_fe_ko_hfd39+SEM_AUC_fe_ko_hfd39)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
% plot(linspace(-0.05,0.05,3),(mean_AUC_fe_ko_hfd39-SEM_AUC_fe_ko_hfd39)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
% plot(zeros(1,3), mean_AUC_ma_ko_hfd39*ones(1,3)+SEM_AUC_ma_ko_hfd39*[-1 0 1],'Color',Mcol,'LineWidth', 2)
% plot(0,mean_AUC_ma_ko_hfd39,'o','MarkerEdgeColor',Mcol,'MarkerFaceColor',Mcol,'MarkerSize', 6)
% plot(linspace(-0.05,0.05,3),(mean_AUC_ma_ko_hfd39+SEM_AUC_ma_ko_hfd39)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
% plot(linspace(-0.05,0.05,3),(mean_AUC_ma_ko_hfd39-SEM_AUC_ma_ko_hfd39)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
% xlim(limsx) 
% ylim(limsy) 
% xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
% ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
% title("HFD, 39 weeks: AUC glucose (mM*min)")
% 
% figure(5)
% clf
% [a0,b0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd54-SEM_AUC_fe_wt_hfd54], 'LineColor', Fcol,'LineStyle',':');
% [a2,b2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd54+SEM_AUC_fe_wt_hfd54], 'LineColor', Fcol,'LineStyle',':');
% [c0,d0]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd54-SEM_AUC_ma_wt_hfd54], 'LineColor', Mcol,'LineStyle',':');
% [c2,d2]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd54+SEM_AUC_ma_wt_hfd54], 'LineColor', Mcol,'LineStyle',':');
% contour(epsList2, SIList2, AUCglu, 'LineWidth', 2, 'ShowText', 'on')
% colormap([linspace(0,0.8,250)' linspace(0,0.8,250)' linspace(0,0.8,250)'])
% hold on
% [a1,b1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_fe_wt_hfd54], 'LineColor', Fcol,'LineStyle',':', 'ShowText', 'on');
% f1=fill([a0(1,2:length(a0(1,:))-50),fliplr(a2(1,2:length(a2(1,:))))],[a0(2,2:length(a0(1,:))-50),fliplr(a2(2,2:length(a2(1,:))))],Fcol,'EdgeColor','none');
% f1.FaceAlpha=0.25;
% [c1,d1]=contour(epsList2, SIList2, AUCglu,'LineWidth', 2, 'LevelList', [mean_AUC_ma_wt_hfd54], 'LineColor', Mcol,'LineStyle',':', 'ShowText', 'on');
% f2=fill([c0(1,2:length(c0(1,:))),fliplr(c2(1,2:length(c2(1,:))))],[c0(2,2:length(c0(1,:))),fliplr(c2(2,2:length(c2(1,:))))],Mcol,'EdgeColor','none');
% f2.FaceAlpha=0.25;
% plot(meanEpsilonFe*ones(1,length(SIList2)), SIList2,'Color',Fcol,'LineWidth', 2,'LineStyle',':')
% f3=fill([(meanEpsilonFe-semEpsFe)*ones(1,length(SIList2)),(meanEpsilonFe+semEpsFe)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Fcol,'EdgeColor','none');
% f3.FaceAlpha=0.1;
% plot(meanEpsilonMa*ones(1,length(SIList2)), SIList2,'Color',Mcol,'LineWidth', 2,'LineStyle',':')
% f4=fill([(meanEpsilonMa-semEpsMa)*ones(1,length(SIList2)),(meanEpsilonMa+semEpsMa)*ones(1,length(SIList2))],[SIList2, fliplr(SIList2)],Mcol,'EdgeColor','none');
% f4.FaceAlpha=0.1;
% plot(zeros(1,3), mean_AUC_fe_ko_hfd54*ones(1,3)+SEM_AUC_fe_ko_hfd54*[-1 0 1],'Color',Fcol,'LineWidth', 2)
% plot(0,mean_AUC_fe_ko_hfd54,'o','MarkerEdgeColor',Fcol,'MarkerFaceColor',Fcol,'MarkerSize', 6)
% plot(linspace(-0.05,0.05,3),(mean_AUC_fe_ko_hfd54+SEM_AUC_fe_ko_hfd54)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
% plot(linspace(-0.05,0.05,3),(mean_AUC_fe_ko_hfd54-SEM_AUC_fe_ko_hfd54)*ones(1,3),'LineWidth', 2,'Color',Fcol,'MarkerSize', 10)
% plot(zeros(1,3), mean_AUC_ma_ko_hfd54*ones(1,3)+SEM_AUC_ma_ko_hfd54*[-1 0 1],'Color',Mcol,'LineWidth', 2)
% plot(0,mean_AUC_ma_ko_hfd54,'o','MarkerEdgeColor',Mcol,'MarkerFaceColor',Mcol,'MarkerSize', 6)
% plot(linspace(-0.05,0.05,3),(mean_AUC_ma_ko_hfd54+SEM_AUC_ma_ko_hfd54)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
% plot(linspace(-0.05,0.05,3),(mean_AUC_ma_ko_hfd54-SEM_AUC_ma_ko_hfd54)*ones(1,3),'LineWidth', 2,'Color',Mcol,'MarkerSize', 10)
% xlim(limsx) 
% ylim(limsy) 
% xlabel("Beta cell insulin sensitivity, S_{\beta} (nM^{-1})")
% ylabel("Peripheral insulin sensitivity, S_P (nM^{-1}d^{-1})")
% title("HFD, 54 weeks: AUC glucose (mM*min)")
