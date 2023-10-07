%% %%  Yang Pei. University of Houston. ypei1.work@gmail.com
clear 
global x0 tech kappa  deficit labor w
global  theta_all N J   gama_labor  gama_njk sigma alpha maxit tol tolp vfactor 
N=4; J=3;  vfactor  = 0.5;                    %Number of countries.
maxit    = 1E+9;vfactor  = 0.1; % maximum iteration
tol      = 1E-07;maxit    = 1E+10;  tolp = tol*1e-5;
gama_labor=ones(J,N);      % the share of value added in production function
gama_njk=ones(N*J,J);      % each column k of country n contains production function coefficient of intermediate good k of country n.
tech=ones(J,N);            % mean papameter of productivity in sector-country dimension
theta_all=ones(J,1);                  % Sectoral trade elasticities parameter 
kappa=ones(J*N,N);   % iceberg cost,>1 %tariff, <1. %Trade cost=1+tariff * iceberg cost
deficit=zeros(N,1);        % Trade deficit
labor=ones(N,1);           % labor in country level
w=ones(N,1);               % wage in country level
price=ones(J,N);           % prices in country-sector level
P=ones(N,1);               % prices in country level
sigma=2*ones(J,1);         % ELS within sector (sectoral intermediate goods)     
alpha=ones(J,N);           % the share of each sector in final demand by country (parameter of untility function) 
%%  
theta_all=4*ones(J,1);  % JX1
sigma=2*ones(J,1);      % JX1
alpha =repmat(1/3,J,N); %JXN
gama_labor=(1/3)*ones(J,N); % JXN
gama_njk = repmat( repmat([2/9; 2/9; 2/9],1,J) ,N,1); % NJJ
labor=  [5 ;5 ;10 ;10]; % NX1
deficit = zeros(N, 1);  % NX1
tech= [1 5  5  10; 2 3 4 6; 1 2 4 5]; % JXN
kappa = repmat(ones(N,N),J,1);
% global  theta_all N J   gama_labor  gama_njk sigma alpha maxit tol tolp vfactor 
%%
 world_gdp = 100;
 x0= world_gdp*ones(1,N)/ (ones(1,N)*labor); % sum(.*labor)
 
[wage_level,SecPrice_level,pie_level,...
    CouPrice_level,Realwage_level,Xjn_level,RealIncome_level,...
    RealIncome_perlevel,cost_level,worldGDP_level,itppie,itwage,wfmax,pfmax] ...
    = fcp1(x0,tech, kappa, labor,deficit) ;