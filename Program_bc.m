
format compact
clear,clc
close all

disp(' ')
disp('We use the Shiller data') 
disp('---------------------------------------------------------------')

% mac xlwrite

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

%% 2. Statistics (1872-2019)

newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];

%% 3 Short-Run Predictability + BC

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

for j=1:length(h);
    for i=1:size(ECON_AR,2);
   
        % Bivariate regression models
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;

disp('In-sample results: individual bc predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b3');
 
 
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b4');
 
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b5');
 
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b6');
 
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b7');
 
 
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h3');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h4');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h5');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h6');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h7');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
disp('In-sample results: individual predictors');
disp(beta_hat);
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b9');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
 
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b10');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b11');


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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
 
 
% Write in-sample results to Excel file
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b12');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'b13');
 
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
 
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h9');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
% Write in-sample results to Excel file
 
output_file='Dividend_results.xls';
output_sheet='Short-Run Predictability';

output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h10');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
% Write in-sample results to Excel file
 

output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h11');


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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
  
% Write in-sample results to Excel file
 

output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h12');
 
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
            X_i_j=ECON_AR(1:end,i);
            [b,t,Tbeta2] =ols_bc(r_h(1:end-(h(j)-1),j),X_i_j);
            beta_hat(i,:,j)=[b(2) Tbeta2 b(3) t(3)]';
    end;
end;
 
 
% Write in-sample results to Excel file
 

output_file='Dividend_results.xls';
output_sheet='Short-Run BC Predictability';
xlwrite(output_file,beta_hat,output_sheet,'h13');
  