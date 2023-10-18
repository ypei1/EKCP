function [ImportShareGDP, labor_level, wage_level,SecPrice_level,pie_level,...
   Realwage_level,X_level,worldGDP_level,worldGDP_hat,itwage,wfmax, Dev] ...
    = SYZ_snake(x0,theta,gama_S, gama_t_S,A,tech_SN,kappa_SIN,...
    labor_N, deficit_N,maxit,tol,vfactor,N)
global w0 labor_N gama_S gama_t_S pie0
%%
w0 = x0; 
%%
p0=ones(1,N); wfmax = 1;  itwage   = 1; 
while (itwage <= maxit) && (wfmax > tol)
%%       
[p0,pie0] = PHL(w0,theta,gama_S, gama_t_S,A,tech_SN,kappa_SIN);
X00 = labor_N.*w0 + deficit_N;
%%
[LaborS_Nnew, Dev] = Labor_Stage(X00, pie0 ,w0, labor_N, gama_S, gama_t_S);
%% 
X = X00;
E_n1L1 = X(2)*pie0(1,2)*gama_t_S(2) ; 
E_n2L1 = 0 ;
E_n1L2 = X(2)*pie0(2,2)*gama_t_S(1) + X(1)*pie0(2,1)*gama_t_S(1); 
E_n2L2 = X(1)*pie0(2,1)*gama_t_S(2) ;
E_n1L3 = X(2)*pie0(3,2)*gama_t_S(2); 
E_n2L3 = X(1)*pie0(3,1)*gama_t_S(1) + X(2)*pie0(3,2)*gama_t_S(1) ;
E_n1L4 = 0 ; 
E_n2L4 = X(1)*pie0(4,1)*gama_t_S(2);
%
M_n1L1 = 0 ; 
M_n2L1 = X(2)*pie0(1,2)*gama_t_S(2) ;
M_n1L2 = X(1)*pie0(2,1)*gama_t_S(2); 
M_n2L2 = X(2)*pie0(2,2)*gama_t_S(1) + X(1)*pie0(2,1)*gama_t_S(1);
M_n1L3 = X(1)*pie0(3,1)*gama_t_S(1) + X(2)*pie0(3,2)*gama_t_S(1); 
M_n2L3 = X(2)*pie0(3,2)*gama_t_S(2);
M_n1L4 = X(1)*pie0(4,1)*gama_t_S(2); 
M_n2L4 = 0 ;
%
 IMn1 = M_n1L1 + M_n1L2 + M_n1L3 + M_n1L4; 
 IMn2 = M_n2L1 + M_n2L2 + M_n2L3 + M_n2L4; 
 EXn1 = E_n1L1 + E_n1L2 + E_n1L3 + E_n1L4;  
 EXn2 = E_n2L1 + E_n2L2 + E_n2L3 + E_n2L4;   
%
 IM = [IMn1 IMn2 ]';
 EX = [EXn1 EXn2 ]';
 deficit = deficit_N';
% 
 SnpNJ=EX-IM+deficit; % sum(SnpNJ)
 w1 = w0 + vfactor*SnpNJ'./labor_N; 
%%
 wfmax=sum(abs(SnpNJ),[1 2]);
 itwage       = itwage + 1;
 w00 = w0;
 w0 = w1; 
 %%
 wage_level = w0;
 SecPrice_level = p0; Realwage_level = w0./p0;
 pie_level = pie0;
 X_level = X;
 worldGDP_old = sum(w00.*labor_N ,[1 2]); 
 worldGDP_level = sum(w1.*labor_N ,[1 2]); 
 worldGDP_hat=worldGDP_level./worldGDP_old ;
 labor_level = LaborS_Nnew;
if itwage/2 == fix(itwage/2 )   disp('wfmax itwage worldGDP_hat:');disp([wfmax itwage worldGDP_hat]);    else end
end

 ImportShareGDP = IM' ./ (labor_level.* wage_level);
end
