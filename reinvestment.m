%  reinvestment.m
%  Convert monthly data to annul data
%  set:  http://www.econ.yale.edu/~shiller/data.htm   
%  Material originally used in Chen, "On the reversal of return and dividend growth predictability: A tale of two periods," 
%  JFE, 2009.  
function newdata=reinvestment(data,rate)

% from monthly dividend data to annual dividend data with rf 

[nobs,nvar] = size(data);

newdata = nan(nobs,3);

for i = 12:nobs
    newdata(i,1) = sum(data(i-11:i,1));
end

for ii = 12:nobs
    for i = 1:11
        newdata(ii,2) = data(ii,1);
        newdata(ii-i,2) = data(ii-i,1).*exp(sum(rate(ii-i+1:ii,1)));
        newdata(ii,3) = sum(newdata(ii-11:ii,2));
    end
end

end




