function [z,z1,z2,z3]=rmsemg(w,m,n,sc,datx);

par=w;
md=par(1,1:m*n);
dv=par(1,m*n+1:2*m*n);
ys=par(1,2*m*n+1:2*m*n+m)';
for i= 1:m,
med(i,:)=md(1,i*n-n+1:i*n);
dev(i,:)=dv(1,i*n-n+1:i*n);
end;


v1=datx(:,1)';
fy=datx(:,2)';
for i=1:length(fy),
v=v1(1,i);
f(1,i)=efbd(v,ys,med,dev); 
end;

z=sqrt((1/length(fy))*sum((fy-f).^2,'omitnan')); 
z1 = sqrt((1/length(fy))*sum((fy-f).^2,'omitnan')); 
z2 = max(abs(fy-f)); % maxima diferencia 
z3 = (1-(var(fy-f)/var(fy))); %varianza relativa del modelo y los datos y solo los datos
z  = (z1*z2)*(1/z3);
%z2 muy pequeña respecto a z3
