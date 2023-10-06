%% %%  Yang Pei. University of Houston. ypei1.work@gmail.com
clear
global vfactor tol tolp maxit x0 y0 alpha_c alpha_i delta inves
global N J  gama_labor gama_capital gama_va gama_njk tech kappa sigma  theta_all  deficit labor capital w r
N=4; J=3;  vfactor  = 0.5;                    %Number of countries.
maxit    = 1E+9;vfactor  = 0.1; % maximum iteration
tol      = 1E-07;maxit    = 1E+10;  tolp = tol*1e-5;
gama_labor=ones(J,N);      % the share of value added in production function
gama_capital = ones(J,N);      % the share of value added in production function
gama_va =ones(J,N);  
gama_njk=ones(N*J,J);      % each column k of country n contains production function coefficient of intermediate good k of country n.
tech=ones(J,N);            % mean papameter of productivity in sector-country dimension
theta_all=ones(J,1);                  % Sectoral trade elasticities parameter 
kappa=ones(J*N,N);   % iceberg cost,>1 %tariff, <1. %Trade cost=1+tariff * iceberg cost
deficit=zeros(1,N);        % Trade deficit
labor=ones(1,N);           % labor in country level
capital = ones(1,N);
w = ones(1,N);               % wage in country level
r = ones(1,N);
price=ones(J,N);           % prices in country-sector level
P=ones(1,N);               % prices in country level
sigma=2*ones(J,1);         % ELS within sector (sectoral intermediate goods)     
alpha=ones(J,N);        % the share of each sector in final demand by country (parameter of untility function) 
 alpha_c =alpha; 
 alpha_i =alpha;
delta =0.05;
%%  
theta_all=4*ones(J,1);  % JX1
sigma=2*ones(J,1);      % JX1
alpha_c =repmat(1/3,J,N); %JXN
alpha_i =repmat(1/4,J,N); %JXN
gama_njk = repmat( repmat([2/9; 2/9; 2/9],1,J) ,N,1); % NJJ
gama_va=(1/3)*ones(J,N); % JXN
gama_labor = [1/2 1/3 1/4 1/2; 1/2.5 1/3 1/3 1/2;1/3 2/3 3/5 1/2];
gama_capital =ones(J,N) - gama_labor;
labor=  [5  5  10  10]; % NX1
capital =  [15  25  5  10]; % NX1
deficit = zeros(1, N);  % NX1
delta = 0.05; 
inves = capital.*delta;
tech_new= [1 5  5  10; 2 3 4 6; 1 2 4 5]; % JXN
kappa_new = repmat(ones(N,N),J,1);
%%
 world_gdp = 100;
 x0 =  world_gdp* ones(1,N)./(ones(1,N)*labor' + ones(1,N)*capital'); 
 y0 =  world_gdp* ones(1,N)./(ones(1,N)*labor' + ones(1,N)*capital');  % x0*labor' +  y0*capital'
 
[wage_level,rent_level, SecPrice_level,pie_level,...
    CouPrice_level_c,CouPrice_level_i,Realwage_level,Xjn_level,cost_level,worldGDP_level,itppie,itwage,wfmax,pfmax] ...
    = fcp2(x0,y0, tech_new, kappa_new,theta_all, N,J, ...
    gama_labor, gama_capital, gama_va,gama_njk,sigma,...
    labor,capital,inves, deficit,alpha_c, alpha_i, maxit,tol,tolp,vfactor)