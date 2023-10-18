function [LaborS_N, dev] = Labor_Stage(X00, pie0 ,w0, labor_N, gama_S, gama_t_S)
LaborS_N = ones(2,2);
MatGamaa1 =  gama_S(1).*gama_t_S(1);
MatGamaa2 =  gama_S(2).*gama_t_S(2);
LaborS_N(1,1) = (MatGamaa1*[pie0(1,1)+pie0(2,1), pie0(1,2)+pie0(2,2)]*X00')/w0(1);
LaborS_N(2,1) = (MatGamaa2*[pie0(1,1)+pie0(3,1), pie0(1,2)+pie0(3,2)]*X00')/w0(1);

LaborS_N(1,2) = (MatGamaa1*[pie0(3,1)+pie0(4,1), pie0(3,2)+pie0(4,2)]*X00')/w0(1);
LaborS_N(2,2) = (MatGamaa2*[pie0(2,1)+pie0(4,1), pie0(2,2)+pie0(4,2)]*X00')/w0(1);
dev = sum(LaborS_N, 1) - labor_N; 
end
