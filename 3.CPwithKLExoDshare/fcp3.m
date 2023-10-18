function  [wage_level,rent_level, SecPrice_level,pie_level,...
    CouPrice_level_c,CouPrice_level_i,Realwage_level,Xjn_level,cost_level,deficit,worldGDP_level,itppie,itwage,wfmax,pfmax] ...
    = fcp3(x0,y0, tech, kappa, Nlabor_n, Lbar_n, K_n_guess,inves, deficit_phi_n , theta_all, N, J,   gama_labor, ...
    gama_capital, gama_va, gama_njk, sigma, alpha_c, alpha_i,  maxit, tol ,tolp, vfactor) 
%% Yang Pei. University of Houston. ypei1.work@gmail.com  
AA=gamma((theta_all+ones(J,1)-sigma)./theta_all); 
A=AA.^(ones(J,1)./(ones(J,1)-sigma));                 
Gammar=ones(J,N) ;  
for j   = 1:1:J
     for n   = 1:1:N
         test=[gama_labor(j,n)*gama_va(j,n) gama_capital(j,n)*gama_va(j,n)  gama_njk(J*(n-1)+1:1:J*(n-1)+J,j)'];
         test=(test).^(-test);
         Gammar(j,n)=prod(test,[1 2]);
     end
end        
%%
worldG = x0*Nlabor_n' +  y0*K_n_guess';
 w0  = x0;       r0 = y0;
%%
wfmax = 1;  itwage   = 1;  
while (itwage <= maxit) && (wfmax > tol)
pf0=ones(J,N); 
[pf0,cost_level,pie1,pfmax,it] = PH3(w0,r0,theta_all,gama_labor, gama_capital, gama_va, gama_njk,kappa,J,N,maxit,tolp,Gammar,A,tech) ;%pf0,cost1,pie1,pfmax,it
%%
   CouPrice_level_i = prod((pf0./alpha_i).^alpha_i,1);
 %  I_va =  CouPrice_level_i.* inves;
 %  C_va = w0.*labor +r0.*capital +deficit - I_va;
 % delt_i = (reshape(alpha_i,N*J,1)).*(kron((I_va)',ones(J,1)));
 % delt_c = (reshape(alpha_c,N*J,1)).*(kron((C_va)',ones(J,1)));
 % delt = delt_i + delt_c;
 % pie_trans1=zeros(N,J*N);    
 deficit = -deficit_phi_n.*(r0.*K_n_guess + w0.*Nlabor_n )+ (Lbar_n ./ sum(Lbar_n)).*sum(deficit_phi_n.*(r0.*K_n_guess + w0.*Nlabor_n )); % N X 1 
 
   C_va = w0.*Nlabor_n +r0.*K_n_guess +deficit - ( CouPrice_level_i.* inves);
  delt_i = (reshape(alpha_i,N*J,1)).*(kron(( CouPrice_level_i.* inves)',ones(J,1)));
  delt_c = (reshape(alpha_c,N*J,1)).*(kron((w0.*Nlabor_n +r0.*K_n_guess +deficit - ( CouPrice_level_i.* inves))',ones(J,1)));
 % delt = delt_i + delt_c;
  pie_trans1=zeros(N,J*N);    
    for j=1:1:J
    for n=1:1:N
        pie_trans1(:,J*(n-1)+j)=pie1(N*(j-1)+n,:);
    end
    end
 Xjn_trans = ((eye(N*J)-kron(ones(1,N),gama_njk).*kron(pie_trans1,ones(J,1)))^-1)*(delt_i + delt_c);
 Xjn = reshape(Xjn_trans,J,N);
 %% 
Xjni=zeros(J*N,N);  Xjin=zeros(J*N,N);  
Xjni = pie1.* repmat(reshape(Xjn',N*J,1),1,N) ; 
 for j =1:1:J 
        Xjin(1+N*(j-1):N+N*(j-1), : )= Xjni(1+N*(j-1):N+N*(j-1), : )';
 end

 IM=sum(reshape(sum(Xjni,2),N,J),2)';%1 N  n import
 EX=sum(reshape(sum(Xjin,2),N,J),2)';%1 N  n output
  %ZW=(w1-w0);
 SnpN=EX-IM+deficit;
  Znj1 = (EX-IM+deficit) ; %1XN sum(Znj)
  Znj2_labor = sum( gama_va.*gama_labor.* reshape(sum( reshape( Xjni ,N,N*J ), [1]),J,N ),[1] )- Nlabor_n.*w0;
  Znj2_capital = sum( gama_va.*gama_capital.* reshape(sum( reshape( Xjni ,N,N*J ), [1]),J,N ),[1] )- K_n_guess.*r0;
  Znj2 = Znj2_labor + Znj2_capital;
 %%
   w11= (w0).*(ones(1,N) + vfactor * Znj2_labor  ./ Nlabor_n );%JXN sum(WLnj) sum(w1.*labor')
   r11= (r0).*(ones(1,N) + vfactor * Znj2_capital./ K_n_guess);%JXN sum(WLnj) sum(w1.*labor')
   w1 =  worldG * w11 ./( w11*Nlabor_n' +  r11*K_n_guess'); 
   r1 =  worldG * r11 ./( w11*Nlabor_n' +  r11*K_n_guess'); %w1*labor' +  r1*capital'
   r0 = r1;  
   w0 = w1; 
%%
 wfmax=sum(abs(SnpN./EX),[1 2]);
 itwage       = itwage + 1;
 w0 = w1; 
 r0 = r1; 
%%
worldGDP = w1*Nlabor_n' +  r1*K_n_guess';
if itwage/50 == fix(itwage/50 )   disp(' wfmax   :');disp([ wfmax   ]); ...
  disp('  itwage  :');disp([  itwage  ]); ...
  disp('  worldGDP:');disp([   worldGDP]); 
else end
end

%% 
wage_level=w0;%1XN
rent_level=r0;
SecPrice_level=pf0;%JXN CouPrice_level;
 CouPrice_level_c = prod((pf0./alpha_c).^alpha_c,1);
% CouPrice_level_i = prod((pf0./alpha_i).^alpha_i,1);
 
pie_level=pie1;%JXNXN 
Xjn_level=Xjn ;%expenditure change t1 t0
worldGDP_level=worldGDP; 

 Realwage_level = w0./CouPrice_level_c;
 %RealIncome_level = (w0 .*labor' ) ./CouPrice_level_c;
 %RealIncome_perlevel =   RealIncome_level./labor';
itppie=it;
%%
  
end
