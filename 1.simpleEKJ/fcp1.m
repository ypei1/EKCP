function [wage_level,SecPrice_level,pie_level,...
    CouPrice_level,Realwage_level,Xjn_level,RealIncome_level,...
    RealIncome_perlevel,cost_level,worldGDP_level,itppie,itwage,wfmax,pfmax] ...
    = fcp1(x0,tech, kappa, labor,deficit) 
%% Yang Pei. University of Houston. ypei1.work@gmail.com
global  theta_all N J   gama_labor  gama_njk sigma alpha maxit tol tolp vfactor 
global A Gammar
AA=gamma((theta_all+ones(J,1)-sigma)./theta_all); 
A=AA.^(ones(J,1)./(ones(J,1)-sigma));                 
Gammar=ones(J,N) ;  
for j   = 1:1:J
     for n   = 1:1:N
         test=[gama_labor(j,n)  gama_njk(J*(n-1)+1:1:J*(n-1)+J,j)'];
         test=(test).^(-test);
         Gammar(j,n)=prod(test,[1 2]);
     end
end        
%%
worldG = sum(x0.*labor');
w0  = x0;     w0 = [  0.2514    0.3681    0.3110    0.3793]
%%
pf0=ones(J,N); %nx0=ones(J,N);
wfmax = 1;  itwage   = 1;  
while (itwage <= maxit) && (wfmax > tol)
%%        w0 = {  w1 ...  wN  } J X N
[pf0,cost_level,pie1,pfmax,it] = PH1(w0,theta_all,gama_labor,gama_njk,kappa,J,N,maxit,tolp,Gammar,A,tech) ;%pf0,cost1,pie1,pfmax,it
%%
  delt=(reshape(alpha,N*J,1)).*(kron((w0'.*labor+deficit),ones(J,1)));
  pie_trans1=zeros(N,J*N);        
    for j=1:1:J
    for n=1:1:N
        pie_trans1(:,J*(n-1)+j)=pie1(N*(j-1)+n,:);
    end
    end
 Xjn_trans = ((eye(N*J)-kron(ones(1,N),gama_njk).*kron(pie_trans1,ones(J,1)))^-1)*delt;
 Xjn = reshape(Xjn_trans,J,N);
 %% 
Xjni=zeros(J*N,N);  Xjin=zeros(J*N,N);  
Xjni = pie1.* repmat(reshape(Xjn',N*J,1),1,N) ; 
 for j =1:1:J 
        Xjin(1+N*(j-1):N+N*(j-1), : )= Xjni(1+N*(j-1):N+N*(j-1), : )';
 end

 IM=sum(reshape(sum(Xjni,2),N,J),2);%NX1  n import
 EX=sum(reshape(sum(Xjin,2),N,J),2);%NX1 n output
  %ZW=(w1-w0);
 SnpN=EX-IM+deficit;
  Znj1 = (EX-IM+deficit)'; %JXN sum(Znj)
  Znj2 = sum( gama_labor.* reshape(sum( reshape( Xjni ,N,N*J ), [1]),J,N ),[1] )- labor'.*w0;
 %%
%  rightmatrix = gama_labor.* reshape(sum( reshape( Xjni ,N,N*J ), [1]),J,N );  % rightmatrix = zeros(J,N);
%  Znj2 = sum( rightmatrix,[1] )- labor'.*w0;
 %%
 w1= (w0).*(ones(1,N)+vfactor*Znj2./( labor'));%JXN sum(WLnj) sum(w1.*labor')
 w1 = worldG  * w1/ sum(w1.*labor'); %
%%
 wfmax=sum(abs(SnpN./EX),[1 2]);
 itwage       = itwage + 1;
 w0 = w1; 
%%
  
 
worldGDP = sum(w1.*labor' ,[1 2]);  

if itwage/50 == fix(itwage/50 )   disp(' wfmax   :');disp([ wfmax   ]); ...
  disp('  itwage  :');disp([  itwage  ]); ...
  disp('  worldGDP:');disp([   worldGDP]); 
else end


end

%% 
wage_level=w0;%1XN
SecPrice_level=pf0;%JXN CouPrice_level;
pie_level=pie1;%JXNXN 
Xjn_level=Xjn ;%expenditure change t1 t0
worldGDP_level=sum(wage_level.*labor',[1 2]); 
 CouPrice_level = prod((pf0./alpha).^alpha,1);
 Realwage_level = w0./CouPrice_level;
 RealIncome_level = (w0 .*labor' ) ./CouPrice_level;
 RealIncome_perlevel =   RealIncome_level./labor';
itppie=it;
%%
  
end
