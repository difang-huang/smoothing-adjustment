function [beta,se]=resid_bootstrap(y,x,bn)
n=length(y);
xx=[ones(n,1),x];
b=regress(y,xx);
xhat=xx*b;
ehat=y-xhat;
bb2=ones(bn,2);
for i=1:bn
v=randi([1,n],n,1);
ee=ehat(v);
yy=xhat+ee;
bb2(i,:)=regress(yy,xx);
end
beta=b;
se=std(bb2);
end



