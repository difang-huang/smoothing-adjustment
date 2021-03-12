format compact
clear,clc
close all
 
disp(' ')
disp('We use the Shiller data') 
disp('---------------------------------------------------------------')
 
% mac xlwrite
 return
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
ep_month = (E./P);
 
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
 
 
%% 2 Out of Sample Month
 
%% (1) month data: start from 1872/R=10*12
newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_month(longdates1,:);
 
r = data_long(:,2);
h = [1 3 6 12 15];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 20*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Month';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b3');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e3');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h3');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k3');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n3');
 
 
 
 
%% (2) month data: start from 1872/R=20*12
newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_month(longdates1,:);
 
r = data_long(:,2);
h = [1 3 6 12 15];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 30*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Month';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b4');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e4');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h4');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k4');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n4');
 
 
 
 
%% (3) month data: start from 1872/R=30*12
newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_month(longdates1,:);
 
r = data_long(:,2);
h = [1 3 6 12 15];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 40*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Month';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b5');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e5');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h5');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k5');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n5');
 
 
 
 
 
 
 
 
%% (1) month data: start from 1927/R=10*12
newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_month(longdates2,:);
 
r = data_long(:,2);
h = [1 3 6 12 15];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 20*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Month';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b7');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e7');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h7');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k7');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n7');
 
 
 
 
%% (2) month data: start from 1927/R=20*12
newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_month(longdates2,:);
 
r = data_long(:,2);
h = [1 3 6 12 15];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 30*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Month';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b8');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e8');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h8');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k8');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n8');
 
 
 
 
%% (3) month data: start from 1927/R=30*12
newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_month(longdates2,:);
 
r = data_long(:,2);
h = [1 3 6 12 15];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 40*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Month';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b9');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e9');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h9');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k9');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n9');
 
 
 
 
 
 
 
 
 
 
 
%% (1) month data: start from 1945/R=10*12
newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_month(longdates3,:);
 
r = data_long(:,2);
h = [1 3 6 12 15];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 20*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Month';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b11');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e11');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h11');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k11');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n11');
 
 
 
 
%% (2) month data: start from 1945/R=20*12
newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_month(longdates3,:);
 
r = data_long(:,2);
h = [1 4 8 12 16];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 30*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Month';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b12');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e12');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h12');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k12');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n12');
 
 
 
 
%% (3) month data: start from 1927/R=30*12
newdata_month = [dp1_month dg1_month uep1_month dp2_month  dg2_month  uep2_month ep_month];
time = datetime(1872,1,1):calmonths(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,12,1));
date3 = find(time== datetime(1945,12,1));
date4 = find(time== datetime(2019,12,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_month(longdates3,:);
 
r = data_long(:,2);
h = [1 4 8 12 16];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 40*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Month';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b13');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e13');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h13');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k13');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n13');
 
 
 
 
 
 
 
 
 
%% 3 Out of Sample Annual
 
%% (1) month data: start from 1872/R=10
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_annual(longdates1,:);
 
r = data_long(:,2);
h = [1 2 3 4 5];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 40; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b3');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e3');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h3');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k3');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n3');
 
 
 
 
%% (2) month data: start from 1872/R=20
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_annual(longdates1,:);
 
r = data_long(:,2);
h = [1 2 3 4 5];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 20; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b4');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e4');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h4');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k4');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n4');
 
 
 
 
%% (3) month data: start from 1872/R=30
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_annual(longdates1,:);
 
r = data_long(:,2);
h = [1 2 3 4 5];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 30; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b5');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e5');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h5');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k5');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n5');
 
 
 
 
 
 
 
 
%% (1) month data: start from 1927/R=10
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_annual(longdates2,:);
 
r = data_long(:,2);
h = [1 2 3 4 5];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 40; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b7');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e7');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h7');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k7');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n7');
 
 
 
 
%% (2) month data: start from 1927/R=20
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_annual(longdates2,:);
 
r = data_long(:,2);
h = [1 2 3 4 5];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 20; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b8');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e8');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h8');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k8');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n8');
 
 
 
 
%% (3) month data: start from 1927/R=30
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_annual(longdates2,:);
 
r = data_long(:,2);
h = [1 2 3 4 5];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 30; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b9');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e9');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h9');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k9');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n9');
 
 
 
 
 
 
 
 
 
 
 
%% (1) month data: start from 1945/R=10
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_annual(longdates3,:);
 
r = data_long(:,2);
h = [1 2 3 4 5];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 40; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b11');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e11');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h11');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k11');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n11');
 
 
 
 
%% (2) month data: start from 1945/R=20
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_annual(longdates3,:);
 
r = data_long(:,2);
h = [1 2 3 4 5];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 20; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b12');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e12');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h12');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k12');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n12');
 
 
 
 
%% (3) month data: start from 1927/R=30
newdata_annual = [dp1_annual dg1_annual uep1_annual dp2_annual dg2_annual  uep2_annual ep_annual];
time = datetime(1872,1,1):calyears(1):datetime(2019,12,31); time = time';
date1 = find(time== datetime(1872,1,1));
date2 = find(time== datetime(1926,1,1));
date3 = find(time== datetime(1945,1,1));
date4 = find(time== datetime(2019,1,1));
longdates1 = [date1:date4]';
longdates2 = [date2:date4]';
longdates3 = [date3:date4]';
 
data_long = newdata_annual(longdates3,:);
 
r = data_long(:,2);
h = [1 2 3 4 5];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(r,1);
R = 30; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,1,length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon 
           for i = 1:size(FC_PR,2) ;
           % Compute residuals
            b=ols(data_long(1:R+(p-1),7), [ones(R+(p-1),1) data_long(1:R+(p-1),1) ]); 
            UEP=data_long(1:R+(p-1),7) - [ones(R+(p-1),1) data_long(1:R+(p-1),1)]*b.beta;
            
            % Predictive regressions
            X_i_j_p = [ones(R+(p-1)-h(j),1) UEP(1:end-h(j)) ];           
            results4=ols(r_h(2:R+p-h(j),j),X_i_j_p);
 
            % Out-of-sample forecasts
            FC_PR(p,1,j)=[1 UEP(end)]*results4.beta;
           end
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(FC_PR,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(FC_PR,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b13');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e13');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h13');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k13');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n13');
 
 
 
 
 
%% OTHER PREDICTORS GW 
format compact
clear,clc
close all
 
disp(' ')
disp('We use the Goyal data') 
disp('---------------------------------------------------------------')


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
ep_month = (E./P);
 
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

R         = xlsread('PredictorData2019','Monthly','q674:q1789'); % S&P 500 VW returns with dividends (1927:1-2019:12)
r_f_lag   = xlsread('PredictorData2019','Monthly','k673:k1788'); % risk-free rate, one period lagged (1926:12-2019:11)
D4        = xlsread('PredictorData2019','Monthly','c674:c1789'); % dividends (12-month moving sum of S&P 500 dividends)
D4_lag    = xlsread('PredictorData2019','Monthly','c673:c1788');
SP500     = xlsread('PredictorData2019','Monthly','b674:b1789'); % S&P500 index (1927:1-2019:12)
SP500_1   = xlsread('PredictorData2019','Monthly','b673:b1789'); % lagged S&P500 index (1926:12-2019:11)
SP500_lag = xlsread('PredictorData2019','Monthly','b673:b1788'); % S&P500 index (1927:1-2019:12)
EQP      = log( 1+ R(1:end)) - log(1+r_f_lag(1:end)); % log excess return 1
SR       = log( 1+ R(1:end)) ;
DP       = log(D4)-log(SP500);     % log dividend-price ratio 2
DY       = log(D4)-log(SP500_lag); % log dividend-price ratio 3
DG       = log( D4./D4_lag ) ;
E4       = xlsread('PredictorData2019','Monthly','d674:d1789'); % earnings (12-month moving sum of S&P 500 earnings)
EP       = log(E4)-log(SP500);  % log earnings-price ratio 4 
DE       = log(D4)-log(E4);     % log dividend-payout ratio 5
SVAR     = xlsread('PredictorData2019','Monthly','o674:o1789');  % volatility (Monthly sum of squared daily returns on S&P 500 index)6
BM       = xlsread('PredictorData2019','Monthly','e674:e1789');  % book-to-market ratio7
NTIS     = xlsread('PredictorData2019','Monthly','j674:j1789');  % net equity issuing activity 8
TBL      = xlsread('PredictorData2019','Monthly','f674:f1789');  % T-bill rate (3-month Treasury bill yield (secondary market)) 9
LTY      = xlsread('PredictorData2019','Monthly', 'i674:i1789'); % long-term government bond yield 10
LTR      = xlsread('PredictorData2019','Monthly','m674:m1789');  % long-term government bond return 11
 
AAA      = xlsread('PredictorData2019','Monthly','g674:g1789'); % Mondys AAA-rated corporate bond yield 13
BAA      = xlsread('PredictorData2019','Monthly','h674:h1789'); % Moodys BAA-rated corporate bond yield 14
DFY      = BAA-AAA; % default yield spread 15
CORPR    = xlsread('PredictorData2019','Monthly','n674:n1789'); % long-term corporate bond return16
DFR      = CORPR-LTR; % default return spread 17
INF      = xlsread('PredictorData2019','Monthly','l674:l1789'); % inflation, lagged (1926:12-2019:11) 18
% cay      = xlsread('PredictorData2019','Monthly','j326:j1788'); 
E4_lag   = xlsread('PredictorData2019','Monthly','d673:d1788'); % earnings, lagged (1926:3-2010:3)
EG      = log((1/4)*E4)- log((1/4)*E4_lag); % log earnings growth 
DP_SOP   = log(1+(1/4)*D4./SP500); % log (1+D/P)
r_f      = xlsread('PredictorData2019','Monthly','l674:l1789'); % risk-free rate
r_f      = log(1+r_f); % log risk-free rate 
T   = size(EQP,1) ; 
d4  = log(D4) ;
e4  = log(E4) ;
p4  = log(SP500) ; 

% 
ECON = [DP DY EP DE SVAR BM NTIS TBL LTR DFY DFR INF];
 
% 2018.12
 
ECON2 = ECON(1:1104,:);
TMS = LTY - ECON2(:,8); % term spread of long-term goverment bond 12
 
ECON3 = [ECON2 TMS LTY];
ECON = ECON3;
 
% DG = DG(1:1104);

DG = dg1_month(661:1764);
 
r = DG;
h = [1 3 6 12 15];
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
T = size(DG,1);
R = 20*12; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,size(ECON,2),length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon        
        for i=1:size(ECON,2); % loop of Goyal-Welch predictors
            X_i_j_p=[ones(R+(p-1)-h(j),1) ECON(1:R+(p-1)-h(j),i)];
            results_i_j_p=ols(r_h(2:R+p-h(j),j),X_i_j_p);
            FC_PR(p,i,j)=[1 ECON(R+(p-1),i)]*results_i_j_p.beta;
        end;
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(ECON,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(ECON,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS GW';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b3');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e3');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h3');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k3');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n3');
 
 
%% Annual GW
 
format compact
clear,clc
close all
 
disp(' ')
disp('We use the Goyal data') 
disp('---------------------------------------------------------------')


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
ep_month = (E./P);
 
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

R         = xlsread('PredictorData2019','Annual','t58:t150'); % S&P 500 VW returns with dividends (1927:1-2019:4)
r_f_lag   = xlsread('PredictorData2019','Annual','l57:l149'); % risk-free rate, one period lagged (1926:4-2019:3)
D4        = xlsread('PredictorData2019','Annual','c58:c150'); % dividends (12-month moving sum of S&P 500 dividends)
D4_lag    = xlsread('PredictorData2019','Annual','c57:c149'); % dividends (12-month moving sum of S&P 500 dividends)
SP500     = xlsread('PredictorData2019','Annual','b58:b150'); % S&P500 index (1927:1-2019:4)
SP500_1   = xlsread('PredictorData2019','Annual','b57:b150'); % lagged S&P500 index (1926:4-2019:3)
SP500_lag = xlsread('PredictorData2019','Annual','b57:b149'); % S&P500 index (1926:1-2019:4)
EQP    = log( 1+ R(1:end)) - log(1+r_f_lag(1:end)); % log excess return 1
SR     = log( 1+ R(1:end)) ;
rf     = log(1+r_f_lag(1:end)) ;
DG     = log( D4./D4_lag ) ;
DGP    = log( D4./D4_lag )- log(1+r_f_lag(1:end)) ;
DP     = log(D4)-log(SP500); % log dividend-price ratio 2
DPP    = log(D4)-log(SP500) - log(1+r_f_lag(1:end)) ;
DY     = log(D4)-log(SP500_lag); % log dividend-price ratio 3
E4     = xlsread('PredictorData2019','Annual','d58:d150'); % earnings (12-month moving sum of S&P 500 earnings)
EP     = log(E4)-log(SP500); % log earnings-price ratio 4 
DE     = log(D4)-log(E4); % log dividend-payout ratio 5
DEG    = diff(DE);
SVAR   = xlsread('PredictorData2019','Annual','p58:p150'); % volatility (Monthly sum of squared daily returns on S&P 500 index)6
BM     = xlsread('PredictorData2019','Annual','e58:e150'); % book-to-market ratio7
NTIS   = xlsread('PredictorData2019','Annual','k58:k150'); % net equity issuing activity 8
TBL    = xlsread('PredictorData2019','Annual','f58:f150'); % T-bill rate (3-month Treasury bill yield (secondary market)) 9
LTY   = xlsread('PredictorData2019','Annual', 'i58:i150'); % long-term government bond yield 10
LTR    = xlsread('PredictorData2019','Annual','n58:n150'); % long-term government bond return 11
% TMS   = LTY-TBL; % term spread of long-term goverment bond 12
AAA    = xlsread('PredictorData2019','Annual','g58:g150'); % Mondys AAA-rated corporate bond yield 13
BAA    = xlsread('PredictorData2019','Annual','h58:h150'); % Moodys BAA-rated corporate bond yield 14
DFY    = BAA-AAA; % default yield spread 15
CORPR  = xlsread('PredictorData2019','Annual','o58:o150'); % long-term corporate bond return16
DFR     = CORPR-LTR; % default return spread 17
INF     = xlsread('PredictorData2019','Annual','m58:m150'); % inflation, lagged (1926:4-2019:3) 18
cay     = xlsread('PredictorData2019','Annual','j76:j149'); 
REQP    = R - INF ;
E4_lag  = xlsread('PredictorData2019','Annual','d57:d149'); % earnings, lagged (1926:3-2010:3)
GM      = log( (E4./SP500)./(E4_lag./SP500) ); % log earnings growth 
EG      = log(E4./E4_lag) ;
DP_SOP  = log(1+D4./SP500); % log (1+D/P)
r_f     = xlsread('PredictorData2019','Annual','l58:l150'); % risk-free rate
r_f     = log(1+r_f); % log risk-free rate 
time  =   datetime(1927,1,31):calmonths(12):datetime(2019,12,31); time = time';  
 N    =   time(2:end) ; 
 T    =   size(EQP(1:end),1) ;  
 d4   = log(D4) ;
 e4   = log(E4) ;
 p4  = log(SP500) ; 
 k0 = 0.2; k1 = 0.96;   
  
  
 
%% Goyal and Welch Preditors
 
 
 
%
ECON = [DP DY EP DE SVAR BM NTIS TBL LTR DFY DFR INF];
 
% 2018.12
 
ECON2 = ECON(1:92,:);
TMS = LTY - ECON2(:,8); % term spread of long-term goverment bond 12
 
ECON3 = [ECON2 TMS LTY];
ECON = ECON3;
 
DG = dg1_annual(56:147);
 
 
r = DG;
h = 1:5; 
r_h=nan(length(r),length(h));
for j=1:length(h);
    for t=1:length(r)-(h(j)-1);
        r_h(t,j)=mean(r(t:t+(h(j)-1)));
    end;
end;
 
 
T = size(DG,1);
R = 20; % in-sample period
P=T-R;    % out-of-sample period
FC_PM=nan(P,1);
FC_PR=nan(P,size(ECON,2),length(h));
 
% Compute out-of-sample forecasts
 
for p=1:P; % loop of out-of-sample process  
    disp(p);  
    FC_PM(p)=mean(r(1:R+(p-1)));  % Prevailing mean benchmark forecast
    for j=1:length(h); % loop of forecast horizon        
        for i=1:size(ECON,2); % loop of Goyal-Welch predictors
            X_i_j_p=[ones(R+(p-1)-h(j),1) ECON(1:R+(p-1)-h(j),i)];
            results_i_j_p=ols(r_h(2:R+p-h(j),j),X_i_j_p);
            FC_PR(p,i,j)=[1 ECON(R+(p-1),i)]*results_i_j_p.beta;
        end;
    end;
end;
 
% Evaluate forecasts
R2OS_PR=nan(size(ECON,2),2,length(h));
for j=1:length(h);
    actual_j=r_h(R+1:end-(h(j)-1),j);
    u_PM_j=actual_j-FC_PM(1:end-(h(j)-1));
    u_PR_j=kron(ones(1,size(FC_PR,2)),actual_j)-FC_PR(1:end-(h(j)-1),:,j);
    MSFE_PM_j=mean(u_PM_j.^2);
    MSFE_PR_j=mean(u_PR_j.^2);
    R2OS_PR_j=100*(1-MSFE_PR_j/MSFE_PM_j);
    R2OS_PR(:,1,j)=R2OS_PR_j';
    L_j=h(j);
    for i=1:size(ECON,2);
        f_CW_i_j=u_PM_j.^2-u_PR_j(:,i).^2+(FC_PM(1:end-(h(j)-1))-...
            FC_PR(1:end-(h(j)-1),i,j)).^2;
        results_CW_i_j=nwest(f_CW_i_j,ones(length(f_CW_i_j),1),h(j));
        R2OS_PR(i,2,j)=1-normcdf(results_CW_i_j.tstat,0,1);
    end;     
end;
 
disp(R2OS_PR)
 
output_file='Dividend_results.xls';
output_sheet='R2OS GW Annual';
xlwrite(output_file,R2OS_PR(:,:,1),output_sheet,'b3');
xlwrite(output_file,R2OS_PR(:,:,2),output_sheet,'e3');
xlwrite(output_file,R2OS_PR(:,:,3),output_sheet,'h3');
xlwrite(output_file,R2OS_PR(:,:,4),output_sheet,'k3');
xlwrite(output_file,R2OS_PR(:,:,5),output_sheet,'n3');
%  
 
 


