%  Shiller_data.m 
%  Properties of returns and dividend growth using Shiller's data 
%  set:  http://www.econ.yale.edu/~shiller/data.htm   
%  Data read from data sheet from Shiller and updated, headings at top
%  NB:  in his data, prices and rates are start of period, returns and 
%  growth rates are forward looking (see monthly version for comparison).

format compact
clear,clc
close all

disp(' ')
disp('We use the Shiller data') 
disp('---------------------------------------------------------------')

% mac xlwrite

% return

%% 1. Data input

% load data
D = xlsread('PredictorData2019','Monthly','c2:c1789'); % dividends
P = xlsread('PredictorData2019','Monthly','b2:b1789'); % S&P500 index 
E = xlsread('PredictorData2019','Monthly','d2:d1789'); % earnings 
r_f   = xlsread('PredictorData2019','Monthly','k2:k1789'); % risk free rate

%create uninvest and reinvest dividend 
newdata=reinvestment(D,r_f);

D12 = newdata(:,1)/12;
D12re = newdata(:,3)/12;
D12_lag = lag(D12,1); % lagged dividends 
D12re_lag = lag(D12re,1); % lagged dividends 

% create monthly variables (Chen, 2009)
dp1_month = log(D12./P);
dg1_month = log(D12./D12_lag);
dp2_month = log(D12re./P);
dg2_month = log(D12re./D12re_lag);
ep_month = log(E./P);

% computing uep  
lhs = ep_month(13:end);
rhs = [ones(size(ep_month(13:end),1),1) dp1_month(13:end)];
reg = ols(lhs,rhs);
uep1_month = reg.resid;

lhs = ep_month(13:end);
rhs = [ones(size(ep_month(13:end),1),1) dp2_month(13:end)];
reg = ols(lhs,rhs);
uep2_month = reg.resid;

% create annual variables (Chen, 2009) at the end of year (December)
dp1_annual = month2annual(dp1_month);
dg1_annual = month2annual(log(D12./lag(D12,12)));
dp2_annual = month2annual(dp2_month);
dg2_annual = month2annual(log(D12re./lag(D12re,12)));
ep_annual = month2annual(ep_month);

% computing uep  
lhs = ep_annual(2:end);
rhs = [ones(size(ep_annual(2:end),1),1) dp1_annual(2:end)];
reg = ols(lhs,rhs);
uep1_annual = reg.resid;

lhs = ep_annual(2:end);
rhs = [ones(size(ep_annual(2:end),1),1) dp2_annual(2:end)];
reg = ols(lhs,rhs);
uep2_annual = reg.resid;

% data time from 1872
dp1_month = dp1_month(13:end);
dp2_month = dp2_month(13:end);
dg1_month = dg1_month(13:end);
dg2_month = dg2_month(13:end);
ep_month = ep_month(13:end);

dp1_annual = dp1_annual(2:end);
dp2_annual = dp2_annual(2:end);
dg1_annual = dg1_annual(2:end);
dg2_annual = dg2_annual(2:end);
ep_annual = ep_annual(2:end);

% save data.mat
% 
return

%% 2. Statistics (1872-2019)

newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];

% return

%% Table 1

% date for month
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date1:date2]';
longdates3 = [date2:date4]';
longdates4 = [date1:date3]';
longdates5 = [date3:date4]';

data_long1 = newdata_month(longdates1,:);
stats_long1 = [mean(data_long1); std(data_long1)]' 

data_long2 = newdata_month(longdates2,:);
stats_long2 = [mean(data_long2); std(data_long2)]' 

data_long3 = newdata_month(longdates3,:);
stats_long3 = [mean(data_long3); std(data_long3)]' 

data_long4 = newdata_month(longdates4,:);
stats_long4 = [mean(data_long4); std(data_long4)]' 

data_long5 = newdata_month(longdates5,:);
stats_long5 = [mean(data_long5); std(data_long5)]' 

stats_month = [stats_long1 stats_long2 stats_long3 stats_long4 stats_long5];

output_file='Dividend_results.xls';
output_sheet='Summary Statistics';
xlwrite(output_file,stats_month,output_sheet,'b3');

% date for annual
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date1:date2]';
longdates3 = [date2:date4]';
longdates4 = [date1:date3]';
longdates5 = [date3:date4]';

data_long1 = newdata_annual(longdates1,:);
stats_long1 = [mean(data_long1); std(data_long1)]' 

data_long2 = newdata_annual(longdates2,:);
stats_long2 = [mean(data_long2); std(data_long2)]' 

data_long3 = newdata_annual(longdates3,:);
stats_long3 = [mean(data_long3); std(data_long3)]' 

data_long4 = newdata_annual(longdates4,:);
stats_long4 = [mean(data_long4); std(data_long4)]' 

data_long5 = newdata_annual(longdates5,:);
stats_long5 = [mean(data_long5); std(data_long5)]' 

stats_annual = [stats_long1 stats_long2 stats_long3 stats_long4 stats_long5];

output_file='Dividend_results.xls';
output_sheet='Summary Statistics';
xlwrite(output_file,stats_annual,output_sheet,'b11');

%% Table 2 Correlation Matrix

newdata_month = [dp1_month dg1_month ep_month uep1_month];
newdata_annual = [dp1_annual dg1_annual ep_annual uep1_annual];

[r_month,p_month] = corrcoef(newdata_month);
[r_annual,p_annual] = corrcoef(newdata_annual);

output_file='Dividend_results.xls';
output_sheet='Correlation Matrix';
xlwrite(output_file,r_month,output_sheet,'b3');
xlwrite(output_file,r_annual,output_sheet,'b11');
 

%% Figure 1

time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
title_vec = {'DY','DG','EY','AEY'}; 

figure; 
newdata_annual = [dp1_annual dg1_annual ep_annual uep1_annual];
for i = 1:size(newdata_annual,2)
    subplot(2,2,i)
    plot(time,newdata_annual(:,i),'r',time,mean(newdata_annual(:,i))*ones(size(newdata_annual(:,i),1),1),'b--','LineWidth',1.5);recessionplot;
    title(title_vec{i},'fontsize',10);
end
% print -depsc2 figure_annual_plot1.eps
set(gcf, 'Color', 'w') 
export_fig figure_annual_plot1.pdf -r300



title_vec = {'DY','DG','EY','AEY'}; 

figure; 
newdata_annual = [dp2_annual dg2_annual ep_annual uep2_annual];
for i = 1:size(newdata_annual,2)
    subplot(2,2,i)
    plot(time,newdata_annual(:,i),'r',time,mean(newdata_annual(:,i))*ones(size(newdata_annual(:,i),1),1),'b--','LineWidth',1.5);recessionplot;
    title(title_vec{i},'fontsize',10);
end
% print -depsc2 figure_annual_plot2.eps
set(gcf, 'Color', 'w') 
export_fig figure_annual_plot2.pdf -r300


time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
title_vec = {'DY','DG','EY','AEY'}; 
figure; 
newdata_month = [dp1_month dg1_month ep_month uep1_month];


subplot(2,2,1)
plot(time,newdata_month(:,1),'r',time,mean(newdata_month(:,1))*ones(size(newdata_month(:,1),1),1),'b--','LineWidth',1.5);recessionplot;
ylim([-5 -2])
title(title_vec{1},'fontsize',10);
subplot(2,2,2)
plot(time,newdata_month(:,2),'r',time,mean(newdata_month(:,2))*ones(size(newdata_month(:,2),1),1),'b--','LineWidth',1.5);recessionplot;
ylim([-0.04 0.04])
title(title_vec{2},'fontsize',10);
subplot(2,2,3)
plot(time,newdata_month(:,3),'r',time,mean(newdata_month(:,3))*ones(size(newdata_month(:,3),1),1),'b--','LineWidth',1.5);recessionplot;
title(title_vec{3},'fontsize',10);
subplot(2,2,4)
plot(time,newdata_month(:,4),'r',time,mean(newdata_month(:,4))*ones(size(newdata_month(:,4),1),1),'b--','LineWidth',1.5);recessionplot;
title(title_vec{4},'fontsize',10);

% for i = 1:size(newdata_month,2)
%     subplot(2,2,i)
%     plot(time,newdata_month(:,i),'r',time,mean(newdata_month(:,i))*ones(size(newdata_month(:,i),1),1),'b--','LineWidth',1.5);recessionplot;
%     title(title_vec{i},'fontsize',10);
% end
% print -depsc2 figure_month_plot1.eps
set(gcf, 'Color', 'w') 
export_fig figure_month_plot1.pdf -r300


title_vec = {'DY','DG','EY','AEY'}; 
figure; 
newdata_month = [dp2_month dg2_month ep_month uep2_month];


subplot(2,2,1)
plot(time,newdata_month(:,1),'r',time,mean(newdata_month(:,1))*ones(size(newdata_month(:,1),1),1),'b--','LineWidth',1.5);recessionplot;
title(title_vec{1},'fontsize',10);
subplot(2,2,2)
plot(time,newdata_month(:,2),'r',time,mean(newdata_month(:,2))*ones(size(newdata_month(:,2),1),1),'b--','LineWidth',1.5);recessionplot;
title(title_vec{2},'fontsize',10);
subplot(2,2,3)
plot(time,newdata_month(:,3),'r',time,mean(newdata_month(:,3))*ones(size(newdata_month(:,3),1),1),'b--','LineWidth',1.5);recessionplot;
title(title_vec{3},'fontsize',10);
subplot(2,2,4)
plot(time,newdata_month(:,4),'r',time,mean(newdata_month(:,4))*ones(size(newdata_month(:,4),1),1),'b--','LineWidth',1.5);recessionplot;
title(title_vec{4},'fontsize',10);
% print -depsc2 figure_month_plot2.eps
set(gcf, 'Color', 'w') 
export_fig figure_month_plot2.pdf -r300

% for i = 1:size(newdata_month,2)
%     subplot(2,2,i)
%     plot(time,newdata_month(:,i),'r',time,mean(newdata_month(:,i))*ones(size(newdata_month(:,i),1),1),'b--','LineWidth',1.5);recessionplot;
%     title(title_vec{i},'fontsize',10);
% end
% print -depsc2 figure_month_plot2.eps



%% 3 Short-Run Predictability + IVX/EM

%% Table 3

% (1) month results + uninvestment
 
newdata_month = [dg1_month uep1_month dg2_month  uep2_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date1:date2]';
longdates3 = [date2:date4]';
longdates4 = [date1:date3]';
longdates5 = [date3:date4]';
 
data_long1 = newdata_month(longdates1,:);
data_long2 = newdata_month(longdates2,:);
data_long3 = newdata_month(longdates3,:);
data_long4 = newdata_month(longdates4,:);
data_long5 = newdata_month(longdates5,:);
 
% (1.1) 
r = data_long1(:,1);
ECON_AR = data_long1(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test + IVX 2015 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j]; 
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b3');
xlwrite(output_file,qLL_hat,output_sheet,'f3');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b3');
 
 
 
% (1.2) 
r = data_long2(:,1);
ECON_AR = data_long2(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b4');
xlwrite(output_file,qLL_hat,output_sheet,'f4');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b4');
 
% (1.3) 
r = data_long3(:,1);
ECON_AR = data_long3(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));

IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b5');
xlwrite(output_file,qLL_hat,output_sheet,'f5');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b5');
 
% (1.4) 
r = data_long4(:,1);
ECON_AR = data_long4(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));

IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b6');
xlwrite(output_file,qLL_hat,output_sheet,'f6');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b6');
 
% (1.5) 
r = data_long5(:,1);
ECON_AR = data_long5(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));

IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b7');
xlwrite(output_file,qLL_hat,output_sheet,'f7');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b7');
 
 
% (2) month results + reinvestment
 
% (2.1)
r = data_long1(:,3);
ECON_AR = data_long1(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));

IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h3');
xlwrite(output_file,qLL_hat,output_sheet,'l3');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h3');
 
% (2.2) 
r = data_long2(:,3);
ECON_AR = data_long2(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));

IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h4');
xlwrite(output_file,qLL_hat,output_sheet,'l4');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h4');
 
% (2.3) 
r = data_long3(:,3);
ECON_AR = data_long3(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));

IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h5');
xlwrite(output_file,qLL_hat,output_sheet,'l5');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h5');
 
% (2.4) 
r = data_long4(:,3);
ECON_AR = data_long4(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));

IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h6');
xlwrite(output_file,qLL_hat,output_sheet,'l6');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h6');
 
% (2.5) 
r = data_long5(:,3);
ECON_AR = data_long5(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));

IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h7');
xlwrite(output_file,qLL_hat,output_sheet,'l7');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h7');
 
% (3) annual results + uninvestment
 
newdata_annual = [dg1_annual uep1_annual dg2_annual  uep2_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date1:date2]';
longdates3 = [date2:date4]';
longdates4 = [date1:date3]';
longdates5 = [date3:date4]';
 
data_long1 = newdata_annual(longdates1,:);
data_long2 = newdata_annual(longdates2,:);
data_long3 = newdata_annual(longdates3,:);
data_long4 = newdata_annual(longdates4,:);
data_long5 = newdata_annual(longdates5,:);
 
% (3.1) 
r = data_long1(:,1);
ECON_AR = data_long1(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b9');
xlwrite(output_file,qLL_hat,output_sheet,'f9');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b9');
 
% (3.2) 
r = data_long2(:,1);
ECON_AR = data_long2(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b10');
xlwrite(output_file,qLL_hat,output_sheet,'f10');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b10');
 
% (3.3) 
r = data_long3(:,1);
ECON_AR = data_long3(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b11');
xlwrite(output_file,qLL_hat,output_sheet,'f11');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b11');
 
% (3.4) 
r = data_long4(:,1);
ECON_AR = data_long4(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));

IVX_Wald=nan(size(ECON_AR,2),2,length(h));


 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b12');
xlwrite(output_file,qLL_hat,output_sheet,'f12');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b12');
 
% (3.5) 
r = data_long5(:,1);
ECON_AR = data_long5(:,2);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));
 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b13');
xlwrite(output_file,qLL_hat,output_sheet,'f13');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'b13');
 
 
% (4) annual results + reinvestment
 
% (4.1) 
 
r = data_long1(:,3);
ECON_AR = data_long1(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h9');
xlwrite(output_file,qLL_hat,output_sheet,'l9');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h9');
 
% (4.2) 
r = data_long2(:,3);
ECON_AR = data_long2(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h10');
xlwrite(output_file,qLL_hat,output_sheet,'l10');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h10');
 
% (4.3) 
r = data_long3(:,3);
ECON_AR = data_long3(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h11');
xlwrite(output_file,qLL_hat,output_sheet,'l11');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h11');
 
% (4.4) 
r = data_long4(:,3);
ECON_AR = data_long4(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h12');
xlwrite(output_file,qLL_hat,output_sheet,'l12');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h12');
 
% (4.5) 
r = data_long5(:,3);
ECON_AR = data_long5(:,4);
 
h=1; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat=nan(size(ECON_AR,2),4,length(h));
qLL_hat=nan(size(ECON_AR,2),length(h));
IVX_Wald=nan(size(ECON_AR,2),2,length(h));
 
for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
 
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            
            % EM 2006 test
            y_i_j=r_h(2:end-(h(j)-1),j);
            X_i_j=ECON_AR(1:end-h(j),i);
            Z_i_j=ones(length(r_h)-h(j),1);            
            [~,IVX_Wald_i_j,pval_i_j]=Compute_IVX_Wald(r,ECON_AR,h(j),0,0.99);
            IVX_Wald(i,:,j)=[IVX_Wald_i_j pval_i_j];  
            [qLL_hat_i_j]=Compute_qLL_hat(y_i_j,X_i_j,Z_i_j,h(j));
            qLL_hat(i,j)=qLL_hat_i_j; 
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=ECON_AR; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star(b,i,j)=results_i_j_star_b.tstat(2);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat(i,3,j)=sum(beta_hat_tstat_star(:,i,j)>...
            beta_hat(i,2,j))/B;
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h13');
xlwrite(output_file,qLL_hat,output_sheet,'l13');
output_sheet='IVX';
xlwrite(output_file,IVX_Wald,output_sheet,'h13');
  

%%  Figure 2   ROLLING WINDOW FIGURES

date = datetime(1872,1,1):calyears(1):datetime(2019,12,31); date = date';

% (1) annual + uninvestment 

r = dg1_annual;
win=20;
count=0;
for t=win+1:1:size(dg1_annual,1)
    count=count+1;
    One=ones(win,1);
    [b_coef(:,count),sb_coef(:,count),R2_coef(:,count)] = olsgmm(r(t-win+1:t),[One uep1_annual(t-win:t-1)],1,1); 
    [b_bs(:,count),se_bs(:,count)] = wild_bootstrap(r(t-win+1:t), uep1_annual(t-win:t-1),1000); 
end

t_coef = b_bs./se_bs;

% (2) annual + reinvestment 

r = dg2_annual;
win=20;
count=0;
% rolling estimates with one break dp
for t=win+1:1:size(dg2_annual,1)
    count=count+1;
    One=ones(win,1);
    [b_roladj(:,count),sb_roladj(:,count),R2_roladj(:,count)] = olsgmm(r(t-win+1:t),[One uep2_annual(t-win:t-1)],1,1);    
    [b_bs(:,count),se_roladj(:,count)] = wild_bootstrap(r(t-win+1:t), uep2_annual(t-win:t-1),1000); 
end

t_roladj = b_bs./se_roladj;
 
 
% make plots
figure
subplot(3,1,1)
plot(date(win:end-1),b_coef(2,:)','r:','LineWidth',2);
hold on
plot(date(win:end-1),b_roladj(2,:)','b--','LineWidth',2);recessionplot;
ylim([0 0.6])
% hold on
% plot(date(win:end-1),b_coef(2,:)'-se_bs(2,:)','r--','LineWidth',1);
% plot(date(win:end-1),b_coef(2,:)'+se_bs(2,:)','r--','LineWidth',1);
% hold on
% plot(date(win:end-1),b_roladj(2,:)'-se_roladj(2,:)','b--','LineWidth',1);
% plot(date(win:end-1),b_roladj(2,:)'+se_roladj(2,:)','b--','LineWidth',1);recessionplot;
ylabel('Coefficient')
legend('Unreinvested','Reinvested')
subplot(3,1,2)
plot(date(win:end-1),t_coef(2,:),'r:','LineWidth',2);
hold on
plot(date(win:end-1),t_roladj(2,:),'b--','LineWidth',2);recessionplot;
hold on
plot(date(win:end-1),1.96*ones(size(date(win:end-1))),'m:','LineWidth',2);
ylabel('T−statistic')
legend('Unreinvested','Reinvested')
subplot(3,1,3)
plot(date(win:end-1),R2_coef,'r:','LineWidth',2);
hold on
plot(date(win:end-1),R2_roladj,'b--','LineWidth',2);recessionplot;
ylabel('R−squared')
legend('Unreinvested','Reinvested') 
print -depsc2 figure_rolling_annual.eps
 
 
date = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); date = date';

% (1) month + uninvestment 

r = dg1_month;
win=20*12;
count=0;
for t=win+1:1:size(dg1_month,1)
    count=count+1;
    One=ones(win,1);
    [b_coef(:,count),sb_coef(:,count),R2_coef(:,count)] = olsgmm(r(t-win+1:t),[One uep1_month(t-win:t-1)],12,1);    
    [b_bs(:,count),se_bs(:,count)] = wild_bootstrap(r(t-win+1:t), uep1_month(t-win:t-1),1000); 
end

t_coef = b_bs./se_bs;

% (2) month + reinvestment 

r = dg2_month;
win=20*12;
count=0;
% rolling estimates with one break dp
for t=win+1:1:size(dg2_month,1)
    count=count+1;
    One=ones(win,1);
    [b_roladj(:,count),sb_roladj(:,count),R2_roladj(:,count)] = olsgmm(r(t-win+1:t),[One uep2_month(t-win:t-1)],12,1);    
    [b_bs(:,count),se_roladj(:,count)] = wild_bootstrap(r(t-win+1:t), uep2_month(t-win:t-1),1000); 
end

t_roladj = b_bs./se_roladj;
 
% make plots
figure
subplot(3,1,1)
plot(date(win:end-1),b_coef(2,:)','r:','LineWidth',2);
hold on
plot(date(win:end-1),b_roladj(2,:)','b--','LineWidth',2);recessionplot;
% hold on
% plot(date(win:end-1),b_coef(2,:)'-se_bs(2,:)','r--','LineWidth',1);
% plot(date(win:end-1),b_coef(2,:)'+se_bs(2,:)','r--','LineWidth',1);
% hold on
% plot(date(win:end-1),b_roladj(2,:)'-se_roladj(2,:)','b--','LineWidth',1);
% plot(date(win:end-1),b_roladj(2,:)'+se_roladj(2,:)','b--','LineWidth',1);recessionplot;
ylabel('Coefficient')
legend('Unreinvested','Reinvested')
subplot(3,1,2)
plot(date(win:end-1),t_coef(2,:),'r:','LineWidth',2);
hold on
plot(date(win:end-1),t_roladj(2,:),'b--','LineWidth',2);recessionplot;
hold on
plot(date(win:end-1),1.96*ones(size(date(win:end-1))),'m:','LineWidth',2);
ylabel('T−statistic')
legend('Unreinvested','Reinvested')
subplot(3,1,3)
plot(date(win:end-1),R2_coef,'r:','LineWidth',2);
hold on
plot(date(win:end-1),R2_roladj,'b--','LineWidth',2);recessionplot;
ylabel('R−squared')
legend('Unreinvested','Reinvested') 
print -depsc2 figure_rolling_month.eps
 
%% Table 4 Good Time/Bad Time (Plot)
 
% (1) Annual + uninvestment

load Data_Recessions
NBERRec = datetime(Recessions,'ConvertFrom',"datenum",...
    "Format","yyyy-MM-dd");

time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
isrecession = @(x)any(isbetween(x,NBERRec(:,1),NBERRec(:,2)));
idxrecession = arrayfun(isrecession,time);

data_annual = [dg1_annual uep1_annual dg2_annual  uep2_annual];
r = data_annual(:,1);
ECON_AR = data_annual(:,2);
ECON_AR1 = ECON_AR.*idxrecession;
ECON_AR2 = ECON_AR.*(1-idxrecession);
 
h=[1 2 3 4 5]; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=sum(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat_good=nan(size(ECON_AR,2),4,length(h));
beta_hat_bad=nan(size(ECON_AR,2),4,length(h));

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Regression models
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR1(1:end-h(j),i) ECON_AR2(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat_good(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            beta_hat_bad(i,:,j)=[results_i_j.beta(3) results_i_j.tstat(3) ...
                nan(1,1) 100*results_i_j.rsqr];
    end;
end;


% anova

% d1 = ones(length(r_h)-h(j),1);
% d2 = ECON_AR1(1:end-h(j),i);
% d3 = ECON_AR2(1:end-h(j),i);
% d4 = r_h(2:end-(h(j)-1),j);
% t = table(d1,d2,d3,d4);
% mdl12 = fitlm(t,'d4~d2+d3');
% a = anova(mdl12)
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=[ECON_AR1 ECON_AR2]; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR1(1:end-1-(h(j)-1),i) ECON_AR2(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star1(b,i,j)=results_i_j_star_b.tstat(2);
                beta_hat_tstat_star2(b,i,j)=results_i_j_star_b.tstat(3);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat_good(i,3,j)=sum(beta_hat_tstat_star1(:,i,j)>...
            beta_hat_good(i,2,j))/B;
        beta_hat_bad(i,3,j)=sum(beta_hat_tstat_star2(:,i,j)>...
            beta_hat_bad(i,2,j))/B;        
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat_good);
disp(beta_hat_bad);

beta_good = reshape(beta_hat_good(1,1,:),1,5);
beta_bad = reshape(beta_hat_bad(1,1,:),1,5);
t_good = reshape(beta_hat_good(1,2,:),1,5);
t_bad = reshape(beta_hat_bad(1,2,:),1,5);

% figures
figure
subplot(1,2,1)
plot(h,beta_good,'-s','MarkerSize',15)
xticks([1 2 3 4 5])
hold on
plot(h,beta_bad,'--ro','MarkerSize',15)
title('(a) Coefficient of AEY for Good and Bad Times ','Interpreter','latex');
ylabel('$b_{K}$','Interpreter','latex')
xlabel('Horizon K')
legend('Bad Times','Good Times','Location','NorthEast')
grid on;

subplot(1,2,2)
plot(h,t_good,'-s','MarkerSize',15)
xticks([1 2 3 4 5])
hold on
plot(h,t_bad,'--ro','MarkerSize',15)
hold on
plot(h,1.96*ones(size(h)),'m:','LineWidth',2);
title('(b) AEY Test Statistics for Good and Bad Times','Interpreter','latex');
ylabel('Test Statistics')
xlabel('Horizon K')
legend('Bad Times','Good Times','Location','NorthEast')
grid on;

% print -depsc2 figure_goodbad_annual1.eps
set(gcf, 'Color', 'w')
export_fig figure_goodbad_annual1.pdf -r300


% (2) annual + reinvestment
load Data_Recessions
NBERRec = datetime(Recessions,'ConvertFrom',"datenum",...
    "Format","yyyy-MM-dd");

time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
isrecession = @(x)any(isbetween(x,NBERRec(:,1),NBERRec(:,2)));
idxrecession = arrayfun(isrecession,time);

data_annual = [dg1_annual uep1_annual dg2_annual  uep2_annual];
r = data_annual(:,3);
ECON_AR = data_annual(:,4);
ECON_AR1 = ECON_AR.*idxrecession;
ECON_AR2 = ECON_AR.*(1-idxrecession);
 
h=[1 2 3 4 5]; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=sum(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat_good=nan(size(ECON_AR,2),4,length(h));
beta_hat_bad=nan(size(ECON_AR,2),4,length(h));

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Regression models
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR1(1:end-h(j),i) ECON_AR2(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat_good(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            beta_hat_bad(i,:,j)=[results_i_j.beta(3) results_i_j.tstat(3) ...
                nan(1,1) 100*results_i_j.rsqr];
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=[ECON_AR1 ECON_AR2]; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR1(1:end-1-(h(j)-1),i) ECON_AR2(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star1(b,i,j)=results_i_j_star_b.tstat(2);
                beta_hat_tstat_star2(b,i,j)=results_i_j_star_b.tstat(3);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat_good(i,3,j)=sum(beta_hat_tstat_star1(:,i,j)>...
            beta_hat_good(i,2,j))/B;
        beta_hat_bad(i,3,j)=sum(beta_hat_tstat_star2(:,i,j)>...
            beta_hat_bad(i,2,j))/B;        
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat_good);
disp(beta_hat_bad);

beta_good = reshape(beta_hat_good(1,1,:),1,5);
beta_bad = reshape(beta_hat_bad(1,1,:),1,5);
t_good = reshape(beta_hat_good(1,2,:),1,5);
t_bad = reshape(beta_hat_bad(1,2,:),1,5);

% figures
figure
subplot(1,2,1)
plot(h,beta_good,'-s','MarkerSize',15)
xticks([1 2 3 4 5])
hold on
plot(h,beta_bad,'--ro','MarkerSize',15)
title('(a) Coefficient of AEY for Good and Bad Times ','Interpreter','latex');
ylabel('$b_{K}$','Interpreter','latex')
xlabel('Horizon K')
legend('Bad Times','Good Times','Location','NorthEast')
grid on;

subplot(1,2,2)
plot(h,t_good,'-s','MarkerSize',15)
xticks([1 2 3 4 5])
hold on
plot(h,t_bad,'--ro','MarkerSize',15)
hold on
plot(h,1.96*ones(size(h)),'m:','LineWidth',2);
title('(b) AEY Test Statistics for Good and Bad Times','Interpreter','latex');
ylabel('Test Statistics')
xlabel('Horizon K')
legend('Bad Times','Good Times','Location','NorthEast')
grid on;

% print -depsc2 figure_goodbad_annual2.eps
set(gcf, 'Color', 'w')
export_fig figure_goodbad_annual2.pdf -r300


% (3) Month + uninvestment

load Data_Recessions
NBERRec = datetime(Recessions,'ConvertFrom',"datenum",...
    "Format","yyyy-MM-dd");
date = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); date = date';
isrecession = @(x)any(isbetween(x,NBERRec(:,1),NBERRec(:,2)));
idxrecession = arrayfun(isrecession,date);

data_month = [dg1_month uep1_month dg2_month  uep2_month];
 
r = data_month(:,1);
ECON_AR = data_month(:,2);
ECON_AR1 = ECON_AR.*idxrecession;
ECON_AR2 = ECON_AR.*(1-idxrecession);
 
h=[1 4 8 12 16 20 24 36 48 60]; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=sum(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat_good=nan(size(ECON_AR,2),4,length(h));
beta_hat_bad=nan(size(ECON_AR,2),4,length(h));

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Regression models
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR1(1:end-h(j),i) ECON_AR2(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat_good(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            beta_hat_bad(i,:,j)=[results_i_j.beta(3) results_i_j.tstat(3) ...
                nan(1,1) 100*results_i_j.rsqr];
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=[ECON_AR1 ECON_AR2]; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR1(1:end-1-(h(j)-1),i) ECON_AR2(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star1(b,i,j)=results_i_j_star_b.tstat(2);
                beta_hat_tstat_star2(b,i,j)=results_i_j_star_b.tstat(3);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat_good(i,3,j)=sum(beta_hat_tstat_star1(:,i,j)>...
            beta_hat_good(i,2,j))/B;
        beta_hat_bad(i,3,j)=sum(beta_hat_tstat_star2(:,i,j)>...
            beta_hat_bad(i,2,j))/B;        
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat_good);
disp(beta_hat_bad);


beta_good = reshape(beta_hat_good(1,1,:),1,10);
beta_bad = reshape(beta_hat_bad(1,1,:),1,10);
t_good = reshape(beta_hat_good(1,2,:),1,10);
t_bad = reshape(beta_hat_bad(1,2,:),1,10);

% figures
figure
subplot(1,2,1)
plot(h,beta_good,'-s','MarkerSize',15)
xticks([1 4 8 12 16 20 24 36 48 60])
hold on
plot(h,beta_bad,'--ro','MarkerSize',15)
title('(a) Coefficient of AEY for Good and Bad Times ','Interpreter','latex');
ylabel('$b_{K}$','Interpreter','latex')
xlabel('Horizon K')
legend('Bad Times','Good Times','Location','NorthEast')
grid on;

subplot(1,2,2)
plot(h,t_good,'-s','MarkerSize',15)
xticks([1 4 8 12 16 20 24 36 48 60])
hold on
plot(h,t_bad,'--ro','MarkerSize',15)
hold on
plot(h,1.96*ones(size(h)),'m:','LineWidth',2);
title('(b) AEY Test Statistics for Good and Bad Times','Interpreter','latex');
ylabel('Test Statistics')
xlabel('Horizon K')
legend('Bad Times','Good Times','Location','NorthEast')
grid on; 

% print -depsc2 figure_goodbad_month1.eps
set(gcf, 'Color', 'w')
export_fig figure_goodbad_month1.pdf -r300


% (4) Month + reinvestment

load Data_Recessions
NBERRec = datetime(Recessions,'ConvertFrom',"datenum",...
    "Format","yyyy-MM-dd");
date = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); date = date';
isrecession = @(x)any(isbetween(x,NBERRec(:,1),NBERRec(:,2)));
idxrecession = arrayfun(isrecession,date);

data_month = [dg1_month uep1_month dg2_month  uep2_month];
 
r = data_month(:,3);
ECON_AR = data_month(:,4);
ECON_AR1 = ECON_AR.*idxrecession;
ECON_AR2 = ECON_AR.*(1-idxrecession);
 
h=[1 4 8 12 16 20 24 36 48 60]; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=sum(r(t:t+(h(j)-1)));
    end;
end;
 
% Compute in-sample results
 
beta_hat_good=nan(size(ECON_AR,2),4,length(h));
beta_hat_bad=nan(size(ECON_AR,2),4,length(h));

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Regression models
            X_i_j=[ones(length(r_h)-h(j),1) ECON_AR1(1:end-h(j),i) ECON_AR2(1:end-h(j),i)];
            results_i_j=nwest(r_h(2:end-(h(j)-1),j),X_i_j,h(j));
            beta_hat_good(i,:,j)=[results_i_j.beta(2) results_i_j.tstat(2) ...
                nan(1,1) 100*results_i_j.rsqr];
            beta_hat_bad(i,:,j)=[results_i_j.beta(3) results_i_j.tstat(3) ...
                nan(1,1) 100*results_i_j.rsqr];
    end;
end;
 
% Compute fixed-regressor wild bootstrap p-values
 
X_sink=[ECON_AR1 ECON_AR2]; % REASON IN HERE: ALL INFORMATION ALL CONTAINTED 
X_sink=[ones(length(r_h),1) X_sink];
results_sink=ols(r_h(2:end,1),X_sink(1:end-1,:));
epsilon_hat=results_sink.resid;
B=1000;
beta_hat_tstat_star=nan(B,size(ECON_AR,2),length(h));
rng('default'); % for reproducability
for b=1:B;
    disp(b);
    u_star_b=randn(length(r_h)-1,1);
    r_star_b=[r(1) ; mean(r)+epsilon_hat.*u_star_b];
    r_h_star_b=nan(length(r),length(h));
    for j=1:length(h);
        for t=1:length(r)-(h(j)-1);
            r_h_star_b(t,j)=mean(r_star_b(t:t+(h(j)-1)));
        end;
    end;
    for j=1:length(h);
        for i=1:size(ECON_AR,2);
                X_i_j=[ones(length(r_h)-1-(h(j)-1),1) ...
                    ECON_AR1(1:end-1-(h(j)-1),i) ECON_AR2(1:end-1-(h(j)-1),i)];
                results_i_j_star_b=nwest(...
                    r_h_star_b(2:end-(h(j)-1),j),X_i_j,h(j));
                beta_hat_tstat_star1(b,i,j)=results_i_j_star_b.tstat(2);
                beta_hat_tstat_star2(b,i,j)=results_i_j_star_b.tstat(3);
        end;
    end;
end;
for j=1:length(h);
    for i=1:size(ECON_AR,2);
        beta_hat_good(i,3,j)=sum(beta_hat_tstat_star1(:,i,j)>...
            beta_hat_good(i,2,j))/B;
        beta_hat_bad(i,3,j)=sum(beta_hat_tstat_star2(:,i,j)>...
            beta_hat_bad(i,2,j))/B;        
    end;
end;
disp('In-sample results: individual predictors');
disp(beta_hat_good);
disp(beta_hat_bad);


beta_good = reshape(beta_hat_good(1,1,:),1,10);
beta_bad = reshape(beta_hat_bad(1,1,:),1,10);
t_good = reshape(beta_hat_good(1,2,:),1,10);
t_bad = reshape(beta_hat_bad(1,2,:),1,10);

% figures
figure
subplot(1,2,1)
plot(h,beta_good,'-s','MarkerSize',15)
xticks([1 4 8 12 16 20 24 36 48 60])
hold on
plot(h,beta_bad,'--ro','MarkerSize',15)
title('(a) Coefficient of AEY for Good and Bad Times ','Interpreter','latex');
ylabel('$b_{K}$','Interpreter','latex')
xlabel('Horizon K')
legend('Bad Times','Good Times','Location','NorthEast')
grid on;

subplot(1,2,2)
plot(h,t_good,'-s','MarkerSize',15)
xticks([1 4 8 12 16 20 24 36 48 60])
hold on
plot(h,t_bad,'--ro','MarkerSize',15)
hold on
plot(h,1.96*ones(size(h)),'m:','LineWidth',2);
title('(b) AEY Test Statistics for Good and Bad Times','Interpreter','latex');
ylabel('Test Statistics')
xlabel('Horizon K')
legend('Bad Times','Good Times','Location','NorthEast')
grid on; 

% print -depsc2 figure_goodbad_month2.eps
set(gcf, 'Color', 'w')
export_fig figure_goodbad_month2.pdf -r300
