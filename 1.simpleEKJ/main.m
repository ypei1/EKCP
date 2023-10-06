%% %%  Yang Pei. University of Houston. ypei1.work@gmail.com
clear
global vfactor tol tolp maxit x0  alpha  
global N J  gama_labor gama_njk tech kappa sigma  theta_all  deficit labor w
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
% other variables setup
Q=ones(J,N);               % # of sectoral intermediate good 
Xjn=ones(J,N);             % sectorcal spending,$ in sector-country level
pie=ones(J*N,N);           % expenditure share matrix for each tradable sector, (destination country, source country)
E=ones(J*N,N);             % Ejni=btra_jin (source country, destination country),row export to column;Value Country n exports to i
M=ones(J*N,N);             % Mjni=btra_jni(destination country, source country),row import from column;Value Country n imports from i
D=ones(J,N);               % trade deficit, country-sector level Dnj=sum i : Mjni-Ejni
GO=ones(J,N);              % gross output by sector and country, $
IO=ones(N*J,J);            % input-output value for each country,(source country,destination country)
Biltrade=ones(J*N,N);     % bilateral trade matrices,for each sector tj , (destination country,source country)
Biltrade_new=ones(J*N,N);
VA=ones(J,N);              % valued added in sector-country level,$, J*N
VA=GO.*gama_labor  ;        % valued added in sector-country level,$, J*N
VAn=sum(VA)'  ;             % valued added in country level,$,N*1
%%  
theta_all=4*ones(J,1);  % JX1
sigma=2*ones(J,1);      % JX1
alpha =repmat(1/3,J,N); %JXN
gama_labor=(1/3)*ones(J,N); % JXN
gama_njk = repmat( repmat([2/9; 2/9; 2/9],1,J) ,N,1); % NJJ
labor=  [5 ;5 ;10 ;10]; % NX1
deficit = zeros(N, 1);  % NX1

tech_new= [1 5  5  10; 2 3 4 6; 1 2 4 5]; % JXN
kappa_new = repmat(ones(N,N),J,1);
%%
 x0= 10*ones(1,N)/ (ones(1,N)*labor); % sum(.*labor)
 
   
 
  
[wage_level,SecPrice_level,pie_level,...
    CouPrice_level,Realwage_level,Xjn_level,RealIncome_level,...
    RealIncome_perlevel,cost_level,worldGDP_level,itppie,itwage,wfmax,pfmax] ...
    = fcp1(x0,tech_new, kappa_new,theta_all, N,J, ...
    gama_labor,gama_njk,sigma,labor,deficit,alpha,maxit,tol,tolp,vfactor)