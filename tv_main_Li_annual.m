load('data')

window=20;

[tv_beta_rolling1,se_r1,tstat_r1,r2_r1]=rollingols(dg1_annual(2:end),uep1_annual(1:(end-1)),window);
[tv_beta_rolling2,se_r2,tstat_r2,r2_r2]=rollingols(dg2_annual(2:end),uep2_annual(1:(end-1)),window);

xlswrite("dguep_annual.xls",[dg1_annual(2:end),uep1_annual(1:end-1),dg2_annual(2:end),uep2_annual(1:end-1)])

%%%%PLOT results--->
res1=xlsread("results_dg1uep1_annual.csv");
res2=xlsread("results_dg2uep2_annual.csv");
date = datetime(1873,2,1):calmonths(12):datetime(2019,12,31); 
date = date'; 

plotdate=date(20:end);
figure;
subplot(2,1,1);
plot(plotdate,res1(10:end-10,2)','g-.','LineWidth',2);
ylim([-0.1 0.8])
%xlim([1890 2020])
hold on
plot(plotdate,res2(10:end-10,2)','c-','LineWidth',2);
plot(plotdate,tv_beta_rolling1,'r:','LineWidth',2);
plot(plotdate,tv_beta_rolling2,'b--','LineWidth',2);
recessionplot;
dateFormat = 10;
datetick('x',dateFormat)
ylabel('Coefficient')
legend('Time-varying Unreinvested','Time-varying Reinvested','Rolling-window Unreinvested','Rolling-window Reinvested')
plot(plotdate,zeros(length(plotdate),1)*1.96,'r--','HandleVisibility','off');
hold off
subplot(2,1,2);
tstat1_annual=res1(:,5);
tstat2_annual=res2(:,5);
plot(plotdate,tstat1_annual(10:end-10)','g-.','LineWidth',2);
ylim([-1 13])
%xlim([1890 2020])
hold on
plot(plotdate,tstat2_annual(10:end-10)','c-','LineWidth',2);
plot(plotdate,tstat_r1,'r:','LineWidth',2);
plot(plotdate,tstat_r2,'b--','LineWidth',2);recessionplot;
dateFormat = 10;
datetick('x',dateFormat)
%yline(1.96,'r--');
ylabel('t-statistic')
legend('Time-varying Unreinvested','Time-varying Reinvested','Rolling-window Unreinvested','Rolling-window Reinvested')
plot(plotdate,ones(length(plotdate),1)*1.96,'r--','HandleVisibility','off');
hold off

export_fig annual.pdf -r300

 










