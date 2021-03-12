function newdata=month2annual(data)

% from monthly dividend data to annual dividend data with rf 

[nobs,nvar] = size(data);

T = nobs/12;

newdata = nan(T,1);

for i = 1:T
    newdata(i,1) = data(12*i,1);
end


end