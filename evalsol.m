function [f,med,dev,ys]=evalsol(w,m,n,datx);

par=w;
md=par(1,1:m*n);
dv=par(1,m*n+1:2*m*n);
ys=par(1,2*m*n+1:2*m*n+m)';

for i= 1:m,
med(i,:)=md(1,i*n-n+1:i*n);
dev(i,:)=dv(1,i*n-n+1:i*n);
end;


v1=datx;
for i=1:length(v1),
v=v1(i,1);
f(1,i)=efbd(v,ys,med,dev); 
end;



