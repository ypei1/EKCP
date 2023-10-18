,%% 
clear all; clear session; close all; clc;warning off all; 
global maxit tol vfactor N S theta eta A gama_S gama_t_S 
global x0 tech_SN kappa_SIN labor_N deficit_N    
vfactor = 0.5; tol = 1E-10; maxit = 1E+10; N = 2; S = 2; theta = 4; eta=2;   
A=(gamma((theta+1-eta)./theta))^(1./(1-eta));    %1X1
gama_S = ones(S,1);  gama_S = [ 1 ; 0.5 ];%SX1
gama_t_S = ones(S,1);gama_t_S = [ 0.5 ; 1 ];%SXN
% gama_s or gama_tilda_s = 
%   ...     
%  stage1     [ x ]
%  stage2     [ y ]
deficit_N = zeros(1,N); %1XN
labor_N = [10 10]; %1XN
tech_SN = ones(S,N);  tech_SN = [ 1 1 ; 1 1];%SXN
% technology
% tech_SN = 
%   ...      cty1  cty2
%  stage1     [ x1  z2]
%  stage2     [ y1  w2]
kappa_SIN =ones(S*N ,N); 
% trade cost
% kappa_SIN = 
%  stage1     [ 1  d1]
%        ...  [ d1  1]
%  stage2     [ 1  d2]
%        ...  [ d2  1]
x0 = [1 ,1]; x0 = x0/sum(labor_N.*x0) ; % 
sum(x0 .* labor_N)
%%
step = 1000;
High_kappa = 999;
Middle_kappa = 3;
kappa_vector1 = linspace(1,Middle_kappa ,step);
kappa_vector2 = exp(linspace(log(Middle_kappa) ,log(High_kappa),30*step + 1));
kappa_vector = [kappa_vector1(1:step-1) kappa_vector2 ];
welfare_vector = zeros(size(kappa_vector));
tradesharegdp_vector = zeros(size(kappa_vector));
for i = 1 : 31*step
%%
kappa_SIN =kappa_vector(i)*ones(S*N ,N) 
 %---------------------------------------------------------------------------------------------------- 
[ImportShareGDP,labor_level, wage_level,SecPrice_level,pie_level,...
   Realwage_level,X_level,worldGDP_level,worldGDP_hat,itwage,wfmax, Dev] ...
    = SYZ_snake(x0,theta,gama_S, gama_t_S,A,tech_SN,kappa_SIN,...
    labor_N, deficit_N,maxit,tol,vfactor,N);
welfare_vector(i) = Realwage_level(1);
tradesharegdp_vector(i) = ImportShareGDP(1);
 % ---------------------------------------------------------------------------------------------------- 
 if i/10 == fix(i/10 )   disp('i = :');disp([i itwage]);    else end
end

%% Autarky 
kappa_SIN = 999999999999*ones(S*N ,N) 
 %---------------------------------------------------------------------------------------------------- 
[ImportShareGDP,labor_level, wage_level,SecPrice_level,pie_level,...
   Realwage_level,X_level,worldGDP_level,worldGDP_hat,itwage,wfmax, Dev] ...
    = SYZ_snake(x0,theta,gama_S, gama_t_S,A,tech_SN,kappa_SIN,...
    labor_N, deficit_N,maxit,tol,vfactor,N);
welfare_Autarky = Realwage_level(1);
tradesharegdp_Autarky = ImportShareGDP(1);
%%
if wealfare < 1
welfare_vector_absolute = welfare_vector/welfare_Autarky;
sz = 30;size=2;figure(98)          % define figure
%subplot(2,3,1);     % subplot(x,y,n)x表示显示的行数，y表示列数，n表示第几幅图片
plot(kappa_vector(1:5*step),welfare_vector_absolute(1:5*step),'-','Color','b','LineWidth',2);
yline(1,'--');
axis([1 5.5  0.95 1.5]);
xlabel('Trade cost of all stages','FontSize',16)
ylabel('Welfare reletive to Autarky','FontSize',16 )
legend('Welfare reletive to Autarky','NumColumns',3,'FontSize',16) 


%subplot(2,3,1);     % subplot(x,y,n)x表示显示的行数，y表示列数，n表示第几幅图片
plot(kappa_vector(1:step),welfare_vector_absolute(1:step),'-','Color','b','LineWidth',2);
yline(1,'--');
axis([1 3  0.95 1.5]);
xlabel('Trade cost of all stages','FontSize',16)
ylabel('Welfare reletive to Autarky','FontSize',16 )
legend('Welfare reletive to Autarky','NumColumns',3,'FontSize',16) 
%%
end
%%
if opentradeshare < 1
tradesharegdp_vector_absolute = tradesharegdp_vector/tradesharegdp_Autarky;
tradesharegdp_vector_absolute = tradesharegdp_vector ;
sz = 30;size=2;figure(98)          % define figure
%subplot(2,3,1);     % subplot(x,y,n)x表示显示的行数，y表示列数，n表示第几幅图片
plot(kappa_vector(1:6*step),tradesharegdp_vector_absolute(1:6*step),'-','Color','b','LineWidth',2);
yline(1.5,'--');
yline(0,'--');
axis([1 8  -0.2 1.7]);
xlabel('Trade cost of all stages','FontSize',16)
ylabel('Import share of GDP','FontSize',16 )
legend('Import share of GDP','NumColumns',3,'FontSize',16) 


%subplot(2,3,1);     % subplot(x,y,n)x表示显示的行数，y表示列数，n表示第几幅图片
plot(kappa_vector(1:step),tradesharegdp_vector_absolute(1:step),'-','Color','b','LineWidth',2);
yline(1,'--');
axis([1 3   -0.2 1.7]);
xlabel('Trade cost of all stages','FontSize',16)
ylabel('Import share of GDP','FontSize',16 )
legend('Import share of GDP','NumColumns',3,'FontSize',16) 
end