function [pf0,cost1,pie1,pfmax,it]= PH1(w0,theta_all,gama_labor,gama_njk,kappa,J,N,maxit,tolp,Gammar,A,tech)                                 
%%
wf0      =  w0; %(J,N)  -----------> (1,N)wf0      =  nw0
pf0      =  ones(J,N);%(J,N)
pf1      =  ones(J,N);%(J,N)
pfmax = 1; % deviations
it   = 1; %time of iteration
%%   %% Yang Pei. University of Houston. ypei1.work@gmail.com          %                w0 = {  w1 ...  wN  } 1 X N     w0 = {  w1 ...  wN  } J X N
while (it <= maxit) && (pfmax > tolp)
    cost1  = Gammar.*(repmat(wf0,J,1).^gama_labor).*reshape(prod(reshape( repmat(reshape(pf0,N*J,1),1,J).^gama_njk  ,J,J*N),[1]),N,J)';  
    pf1= repmat(A,1,N) .* (reshape((     sum(  kron(tech,ones(N,1)).*(kappa.*kron(cost1,ones(N,1))).^(-kron(theta_all,ones(N,N)))   ,2)        ),N,J)').^(repmat(-theta_all.^(-1),1,N));
%%
    pfdev    = abs((pf1 - pf0)./pf0); % Checking tolerance
    pf0      = pf1;
    pfmax    = max(max(pfdev));
    it       = it + 1;
end    
    pie1=(   ones(J*N,N) ./ kron(A.^theta_all, ones(N,N))  ).* kron(tech,ones(N,1)).*(kappa.*kron(cost1,ones(N,1))./     repmat(reshape(pf1',J*N,1),1,N) ).^(-kron(theta_all,ones(N,N)));    
end
