%Algoritmo gen�tico simple 
%Codificaci�n real
%Selecci�n estoc�stica por ruleta
%Problema: Aproximaci�n de la funci�n SINC con una EFBD

clc;
clear;
close all;

load('spotify_pro_5.mat');

datz=data(:,2:5);

indz=randperm(length(datz)); % permutación aleatoria datos
datx=datz(indz,:);

T=ceil(0.7*length(datz));
V=T+1;

datt=datx(1:T,:); %  70% datos de evolución
datv=datx(T+1:end,:); %30%  %datos de valiación

%salidad de datav
plot(datv(:,1),datv(:,2),'.b');
hold
%iterador en x y 
plot(data(:,1),data(:,3),'r');

 
% Par�metros de la EFBD
rul=24; % No de reglas
inpt=4; %No de entradas 
totpar=2*inpt*rul+rul; % No. total de par�metros
scale=max(data(:,2)); % escala  de los parametros (+/-) scale/2


% Par�metros del algoritmo gen�tico
Nrun=10; %10  % No. de corridas independientes
Ngen=50;% 50 en data %1e4  % No. de generaciones e4
npop=30; % Tama�o de la poblaci�n
press=0.01; % Presi�n selectiva
ipop=60; % Tama�o de la poblaci�n intermedia (parejas)
pcross=0.7; % Probabilidad de cruce
pmut=0.03;% Probabilidad de mutaci�n

for r =1: Nrun
    tic
% inicializaci�n
pop=2*scale*(rand(npop,totpar)-0.5);

% Ciclo evolutivo
for gen =1:Ngen
 
%evaluaci�n de individuos
for i =1:npop
    w=pop(i,:);
    err=rmsemg(w,rul,inpt,scale,datt); % Calculo de la aproximaci�n funci�n SINC
    fobj(i,1)=err;
end
  
fobjMax = max(fobj);
fobjMin = min(fobj);
fobjAv  = mean (fobj);

outpop(gen,1:3,r)=[ fobjMax fobjMin fobjAv];
gen

%calificaci�n de individuos
Efficiency = (1 - press) * (fobjMax - fobj)/max([fobjMax - fobjMin, eps]) + press;

%Selecci�n por ruleta
Wheel = cumsum(Efficiency);
  for j=1:ipop
    %Selecci�n del primer individuo de la pareja
    Shoot = rand(1,1)*max(Wheel);
    Index = 1;
    while((Wheel(Index)<Shoot)&(Index<length(Wheel))) 
      Index = Index + 1;
    end
    indiv1(j,:) = pop(Index,:);
    findiv1(j) = fobj(Index);
    Iindiv1(j) = Index;
    % Selecci�n del segundo individuo de la pareja
    Shoot = rand(1,1)*max(Wheel);
    Index = 1;
    while((Wheel(Index)<Shoot)&(Index<length(Wheel))) 
      Index = Index + 1;
    end
    indiv2(j,:)= pop(Index,:);
    findiv2(j)= fobj(Index);
    Iindiv2(j)= Index;
  end
 
 % Cruce
 
  for j=1:ipop
    if (pcross>rand(1,1))
      pind = ceil((totpar-1)*rand(1,1) + 1);
      x1 = [indiv1(j,1:pind)  indiv2(j,pind+1:totpar)];
      x2 = [indiv2(j,1:pind)  indiv1(j,pind+1:totpar)];
      indiv1(j,:) = x1;
      indiv2(j,:) = x2;
    end
  end
  
  %Mutaci�n
  
  for j=1:ipop
   if (pmut>rand(1,1))
      pind= ceil((totpar-1)*rand(1,1) + 1);
      vind=2*(rand(1,1)-0.5)*scale;
      indiv1(j,pind)=vind;
   end
    if (pmut>rand(1,1))
       pind= ceil((totpar-1)*rand(1,1) + 1);
      vind=2*(rand(1,1)-0.5)*scale;
      indiv2(j,pind)=vind;
    end
  end 
 
 poplast=pop;
 
 %  Nueva poblaci�n 
 for j=1:npop
   indexsel = ceil((ipop-1)*rand(1,1) + 1);
   if rand(1,1) > 0.5
     pop(j,:)=indiv1(indexsel,:);
   else
     pop(j,:)=indiv2(indexsel,:);
   end
 end
 
 Wheel = 0;
 
end

%Obtenci�n del mejor individuo de la poblaci�n final
 for i =1:npop,
     w=poplast(i,:);
     err=rmsemg(w,rul,inpt,scale,datv);
     fobj(i,1)=err;
 end
% 
 [sobj k]=sort(fobj);
% 
 indx=k(1,1);
 m=rul;
 n=inpt;
 par(r,:)=poplast(indx,:);

 clear fobj 
 save(['data1.1/corrida',num2str(r),'.mat'])
 toc
end

sol=1; % corrida de interés

evop=outpop(:,:,sol);
f=evalsol(par(sol,:),m,n,data(:,1));
save('prueba1Data');

%%
clc;clear all;close all;
load('prueba3Data2.mat')
close all
sol=5  %47 %5; % corrida de interés

evop=outpop(:,:,sol);

[f,med,dev,ys] = evalsol(par(sol,:),m,n,data(:,6));
min(min(dev))
min(min(med))
min(ys)

f=f-35;


figure;
plot(data(:,1),data(:,6),'b');
hold;
plot(data(:,1),f,'r');

        grid on;
        set(gca,'FontSize',14);
        set(gcf,'Color','white');
        %title('Generaciones','FontSize',16);
        xlabel('Entrada unica','FontSize',16,'Interpreter','latex');
        ylabel('Valor HIT','FontSize',16,'Interpreter','latex'); 
        legend('Datos validación','Datos algoritmo genético');
        %saveas(gcf,'inp1.png')
        
figure
hist(f);
hold;
h = histogram(data(:,6));
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
        grid on;
        set(gca,'FontSize',14);
        set(gcf,'Color','white');
        title('Histograma','FontSize',16);
        xlabel('Valor','FontSize',16,'Interpreter','latex');
        ylabel('Numero de repeticiones','FontSize',16,'Interpreter','latex'); 
        legend('Datos algoritmo genético','Datos validación');
        %saveas(gcf,'gen1.png')

figure;
plot(evop,'LineWidth',2);

        grid on;
        set(gca,'FontSize',14);
        set(gcf,'Color','white');
        title('Generaciones','FontSize',16);
        xlabel('N Generaciones','FontSize',16,'Interpreter','latex');
        ylabel('Error','FontSize',16,'Interpreter','latex');      
        legend('fobjMax', 'fobjMin', 'fobjAv');
 
        saveas(gcf,'Evop1.png') 
        
[err,z1,z2,z3]=rmsemg(w,rul,inpt,scale,datv);
err
z1
z2
z3
%writematrix([data(:,6),f'], ['proyecto31','.csv']);
%%
load('prueba1Data.mat')
close all
sol=9;%6 % corrida de interés

evop=outpop(:,:,sol);

[f,med,dev,ys] = evalsol(par(sol,:),m,n,data(:,6));
f=f+10;

figure;
plot(data(:,1),data(:,6),'b');
hold;
plot(data(:,1),f,'r');

        grid on;
        set(gca,'FontSize',14);
        set(gcf,'Color','white');
        %title('Generaciones','FontSize',16);
        xlabel('Entrada unica','FontSize',16,'Interpreter','latex');
        ylabel('Valor HIT','FontSize',16,'Interpreter','latex'); 
        legend('Datos validación','Datos algoritmo genético');
        %saveas(gcf,'inp.png')
        
figure
hist(f);
hold;
h = histogram(data(:,6));
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
        grid on;
        set(gca,'FontSize',14);
        set(gcf,'Color','white');
        title('Histograma','FontSize',16);
        xlabel('Valor','FontSize',16,'Interpreter','latex');
        ylabel('Numero de repeticiones','FontSize',16,'Interpreter','latex'); 
        legend('Datos algoritmo genético','Datos validación');
        %saveas(gcf,'gen.png')

figure;
plot(evop,'LineWidth',2);

        grid on;
        set(gca,'FontSize',14);
        set(gcf,'Color','white');
        title('Generaciones','FontSize',16);
        xlabel('N Generaciones','FontSize',16,'Interpreter','latex');
        ylabel('Error','FontSize',16,'Interpreter','latex');      
        legend('fobjMax', 'fobjMin', 'fobjAv');
 
        %saveas(gcf,'Evop.png')       
        
 [err,z1,z2,z3]=rmsemg(w,rul,inpt,scale,datv);
 err
 z1
 z2
 z3
 writematrix([data(:,6),f'], ['proyecto3','.csv']);
%%
% load('prueba2Data1.mat')
% close all
% sol=1; % corrida de interés
% 
% evop=outpop(:,:,sol);
% 
% [f,med,dev,ys] = evalsol(par(sol,:),m,n,data(:,6));
% %f=f-50;
% figure;
% plot(evop);
% 
% figure;
% plot(data(:,1),data(:,6),'b');
% hold;
% plot(data(:,1),f,'r');
% 
% figure
% hist(f);
% hold;
% h = histogram(data(:,6));
% h.FaceColor = [0 0.5 0.5];
% h.EdgeColor = 'w';
% 
% [err,z1,z2,z3]=rmsemg(w,rul,inpt,scale,datv);
% err
%%
clc;clear all;close all;

load('spotify_pro_5.mat');

figure;
plot(data(:,1),data(:,6),'b');

figure;
h = histogram(data(:,6));
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
