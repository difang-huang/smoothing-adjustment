
load('data')

window=240;

[tv_beta_rolling1,se_r1,tstat_r1,r2_r1]=rollingols(dg1_month(2:end),uep1_month(1:(end-1)),window);
[tv_beta_rolling2,se_r2,tstat_r2,r2_r2]=rollingols(dg2_month(2:end),uep2_month(1:(end-1)),window);

xlswrite("dguep_month.xls",[dg1_month(2:end),uep1_month(1:end-1),dg2_month(2:end),uep2_month(1:end-1)])

res1=xlsread("results_dg1uep1.csv");
res2=xlsread("results_dg2uep2.csv");
date = datetime(1872,2,1):calmonths(1):datetime(2019,12,31); 
date = date';

plotdate=date(240:end);
figure;
subplot(2,1,1);
plot(plotdate,res1(120:end-120,2)','g-.','LineWidth',2);
ylim([0 0.06])
hold on
plot(plotdate,res2(120:end-120,2)','c-','LineWidth',2);
plot(plotdate,tv_beta_rolling1,'r:','LineWidth',2);
plot(plotdate,tv_beta_rolling2,'b--','LineWidth',2);
recessionplot;
dateFormat = 10;
datetick('x',dateFormat)
ylabel('Coefficient')
legend('Time-varying Unreinvested','Time-varying Reinvested','Rolling-window Unreinvested','Rolling-window Reinvested')
hold off
subplot(2,1,2);
tstat1_month=res1(:,5);
tstat2_month=res2(:,5);
plot(plotdate,tstat1_month(120:end-120)','g-.','LineWidth',2);
ylim([0 40])
hold on
plot(plotdate,tstat2_month(120:end-120)','c-','LineWidth',2);
plot(plotdate,tstat_r1,'r:','LineWidth',2);
plot(plotdate,tstat_r2,'b--','LineWidth',2);recessionplot;
%yline(1.96,'r--');
dateFormat = 10;
datetick('x',dateFormat)
ylabel('t-statistic')
legend('Time-varying Unreinvested','Time-varying Reinvested','Rolling-window Unreinvested','Rolling-window Reinvested')
plot(plotdate,ones(length(plotdate),1)*1.96,'r--','HandleVisibility','off');
hold off

export_fig monthly.pdf -r300


 







 





