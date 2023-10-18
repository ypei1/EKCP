function [p0,pie0] = PHL(w0,theta,gama_S, gama_t_S,A,tech_SN,kappa_SIN) 
%%
w = w0(1); TFP = tech_SN(1,1); gama = gama_S(1); 
D = 1; gama_t_prod = -theta*gama_t_S(1);
P0_1_1 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(1); TFP = tech_SN(2,1); gama = gama_S(2); 
D = 1; gama_t_prod = -theta*gama_t_S(2);
P0_1_2 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(1); TFP = tech_SN(1,1); gama = gama_S(1); 
D = kappa_SIN(2,1); gama_t_prod = -theta*gama_t_S(1);
P0_1_3 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(2); TFP = tech_SN(2,2); gama = gama_S(2); 
D = kappa_SIN(3,2); gama_t_prod = -theta*gama_t_S(2);
P0_1_4 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(2); TFP = tech_SN(1,2); gama = gama_S(1); 
D = kappa_SIN(1,2); gama_t_prod = -theta*gama_t_S(1);
P0_1_5 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(1); TFP = tech_SN(2,1); gama = gama_S(2); 
D = 1; gama_t_prod = -theta*gama_t_S(2);
P0_1_6 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(2); TFP = tech_SN(1,2); gama = gama_S(1); 
D = 1; gama_t_prod = -theta*gama_t_S(1);
P0_1_7 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(2); TFP = tech_SN(2,2); gama = gama_S(2); 
D = kappa_SIN(3,2); gama_t_prod = -theta*gama_t_S(2);
P0_1_8 = ((w/TFP)^gama*D)^gama_t_prod ;
%%
 p0_1 = A*(P0_1_1*P0_1_2 + P0_1_3*P0_1_4 + P0_1_5*P0_1_6 + P0_1_7*P0_1_8)^(-1/theta);
%%
w = w0(1); TFP = tech_SN(1,1); gama = gama_S(1); 
D = 1; gama_t_prod = -theta*gama_t_S(1);
P0_2_1 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(1); TFP = tech_SN(2,1); gama = gama_S(2); 
D = kappa_SIN(4,1); gama_t_prod = -theta*gama_t_S(2);
P0_2_2 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(1); TFP = tech_SN(1,1); gama = gama_S(1); 
D = kappa_SIN(2,1); gama_t_prod = -theta*gama_t_S(1);
P0_2_3 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(2); TFP = tech_SN(2,2); gama = gama_S(2); 
D = 1; gama_t_prod = -theta*gama_t_S(2);
P0_2_4 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(2); TFP = tech_SN(1,2); gama = gama_S(1); 
D = kappa_SIN(1,2); gama_t_prod = -theta*gama_t_S(1);
P0_2_5 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(1); TFP = tech_SN(2,1); gama = gama_S(2); 
D = kappa_SIN(4,1); gama_t_prod = -theta*gama_t_S(2);
P0_2_6 = ((w/TFP)^gama*D)^gama_t_prod ;

w = w0(2); TFP = tech_SN(1,2); gama = gama_S(1); 
D = 1; gama_t_prod = -theta*gama_t_S(1);
P0_2_7 = ((w/TFP)^gama*D)^gama_t_prod ;

w = w0(2); TFP = tech_SN(2,2); gama = gama_S(2); 
D = 1; gama_t_prod = -theta*gama_t_S(2);
P0_2_8 = ((w/TFP)^gama*D)^gama_t_prod ;
%%
 p0_2 = A*(P0_2_1*P0_2_2 + P0_2_3*P0_2_4 + P0_2_5*P0_2_6 + P0_2_7*P0_2_8)^(-1/theta);
%%
 p0 = [p0_1 p0_2];
%% 
pie0 = ones(4,2);   
%%
pp = p0(1); n = 1;
%%
w = w0(1); TFP = tech_SN(1,1); gama = gama_S(1); 
D = 1; gama_t_prod = -theta*gama_t_S(1);
P0_1_1 = ((w/TFP)^gama*D)^gama_t_prod ; 
%
w = w0(1); TFP = tech_SN(2,1); gama = gama_S(2); 
D = 1; gama_t_prod = -theta*gama_t_S(2);
P0_1_2 = ((w/TFP)^gama*D)^gama_t_prod ; 
%
Chain1_n1 = P0_1_1 * P0_1_2 * (pp/A)^theta ;
%%
w = w0(1); TFP = tech_SN(1,1); gama = gama_S(1); 
D = kappa_SIN(2,1); gama_t_prod = -theta*gama_t_S(1);
P0_1_3 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(2); TFP = tech_SN(2,2); gama = gama_S(2); 
D = kappa_SIN(3,2); gama_t_prod = -theta*gama_t_S(2);
P0_1_4 = ((w/TFP)^gama*D)^gama_t_prod ;
%
Chain2_n1 = P0_1_3 * P0_1_4 * (pp/A)^theta ;
%%
%
w = w0(2); TFP = tech_SN(1,2); gama = gama_S(1); 
D = kappa_SIN(1,2); gama_t_prod = -theta*gama_t_S(1);
P0_1_5 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(1); TFP = tech_SN(2,1); gama = gama_S(2); 
D = 1; gama_t_prod = -theta*gama_t_S(2);
P0_1_6 = ((w/TFP)^gama*D)^gama_t_prod ;
%
Chain3_n1 = P0_1_5 * P0_1_6 * (pp/A)^theta ;

%%
%
w = w0(2); TFP = tech_SN(1,2); gama = gama_S(1); 
D = 1; gama_t_prod = -theta*gama_t_S(1);
P0_1_7 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(2); TFP = tech_SN(2,2); gama = gama_S(2); 
D = kappa_SIN(3,2); gama_t_prod = -theta*gama_t_S(2);
P0_1_8 = ((w/TFP)^gama*D)^gama_t_prod ;
%
Chain4_n1 = P0_1_7 * P0_1_8 * (pp/A)^theta ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pp = p0(2); n = 2;
%%
w = w0(1); TFP = tech_SN(1,1); gama = gama_S(1); 
D = 1; gama_t_prod = -theta*gama_t_S(1);
P0_2_1 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(1); TFP = tech_SN(2,1); gama = gama_S(2); 
D = kappa_SIN(4,1); gama_t_prod = -theta*gama_t_S(2);
P0_2_2 = ((w/TFP)^gama*D)^gama_t_prod ;
%
Chain1_n2 = P0_2_1 * P0_2_2 * (pp/A)^theta ;
%%
w = w0(1); TFP = tech_SN(1,1); gama = gama_S(1); 
D = kappa_SIN(2,1); gama_t_prod = -theta*gama_t_S(1);
P0_2_3 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(2); TFP = tech_SN(2,2); gama = gama_S(2); 
D = 1; gama_t_prod = -theta*gama_t_S(2);
P0_2_4 = ((w/TFP)^gama*D)^gama_t_prod ;
%
Chain2_n2 = P0_2_3 * P0_2_4 * (pp/A)^theta ;
%%
%
w = w0(2); TFP = tech_SN(1,2); gama = gama_S(1); 
D = kappa_SIN(1,2); gama_t_prod = -theta*gama_t_S(1);
P0_2_5 = ((w/TFP)^gama*D)^gama_t_prod ;
%
w = w0(1); TFP = tech_SN(2,1); gama = gama_S(2); 
D = kappa_SIN(4,1); gama_t_prod = -theta*gama_t_S(2);
P0_2_6 = ((w/TFP)^gama*D)^gama_t_prod ;
%
Chain3_n2 = P0_2_5 * P0_2_6 * (pp/A)^theta ;

%%
w = w0(2); TFP = tech_SN(1,2); gama = gama_S(1); 
D = 1; gama_t_prod = -theta*gama_t_S(1);
P0_2_7 = ((w/TFP)^gama*D)^gama_t_prod ;

w = w0(2); TFP = tech_SN(2,2); gama = gama_S(2); 
D = 1; gama_t_prod = -theta*gama_t_S(2);
P0_2_8 = ((w/TFP)^gama*D)^gama_t_prod ;
%
Chain4_n2 = P0_2_7 * P0_2_8 * (pp/A)^theta ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pie0 = [Chain1_n1 Chain1_n2;Chain2_n1 Chain2_n2;Chain3_n1 Chain3_n2;Chain4_n1 Chain4_n2]; 
end
