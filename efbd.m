 %función para calcular una EBDF de M términos
% x vector de n entradas
% y vector de m centros
% med matriz con medias de las funciones de base difusas gaussianas
% dev matriz con desviaciones estándar de las funciones de base difusa gaussianas


function [f,b,rul] = efbd1(x,y,med,dev)
m=length(y);
n=length(x);
xp = meshgrid(x(1,(1:n)),1:m);
mp = med((1:m), (1:n));
dp = dev((1:m), (1:n));
rul = prod(exp((-0.5.*(xp-mp).^2)./(dp.^2)),2); 
b = sum(rul);
a = sum(rul.*y);
f = a / max([b, eps]);
