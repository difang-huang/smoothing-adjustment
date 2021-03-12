function [b,t,Tbeta2] = ols_bc(lhv,rhv);
 
% function ols_bc does ols regressions with AH2004 corrected standard errors
% Inputs:
%  lhv T x N vector, left hand variable data 
%  rhv T x K matrix, right hand variable data
%  If N > 1, this runs N regressions of the left hand columns on all the (same) right hand variables. 

% Output:
%  b: regression coefficients K x 1 vector of coefficients
%  t: standard t statistics
%  Tbeta2: corrected t statistics
%   Note: program checks whether first is a constant and ignores that one for test
 
if size(rhv,1) ~= size(lhv,1);
   disp('olsgmm_bc: left and right sides must have same number of rows. Current rows are');
   size(lhv)
   size(rhv)
end;


%% AR(1) for rhv
lags = 0;
weight = 0;
T = size(rhv,1);
[b,se] = olsgmm(rhv(2:end), [ones(T-1,1) rhv(1:end-1)], lags, weight) ;
rho2c  = b(2) + (1+3*b(2))/T + 3*(1+3*b(2))/(T^2) ;    % bias-corrected estimate of rho 
vc_hat = rhv(2:end) - b(1).*ones(T-1,1) - rho2c.*rhv(1:end-1) ;  % bias-corrected residual v_t
VARrho2  = (1 + 3/T + 9/T^2 )* (se(2)).^2 ;         % bias-corrected variance of rho

%% BC
[b,se] = olsgmm(lhv(2:end), [ones(T-1,1) rhv(1:end-1) vc_hat], lags, weight) ;
t = b./se ;
SEbeta2 = sqrt( (b(3))^2 * VARrho2 + (se(2))^2 ) ; % bias-corrected s.e. of beta 
Tbeta2  = b(2)/SEbeta2 ; 
end