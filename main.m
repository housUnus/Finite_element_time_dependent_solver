
hsps = [0.2 0.1 0.08 0.07 0.05 0.01 0.008 0.006];
htime =[0.2 0.1 0.08 0.05 0.02 0.01 0.007 0.003 0.001];
u = @(x,t) x.*(x-4).*(100.*t).*exp(x).*cos(x);
cl1 = [0 0];
cl2 = [0 0];
% ********************* Spatail Validaiton  **************************
 sll = Solver1D_P2_Ref_T_EDP_a(0,4,0.1,0,1,.1,u,1,5,1);
header = [00000 sll.tInt'];
dlmwrite('spatial10^-3.xls',header,'delimiter','\t','precision',3)
for h = hsps
     sll = Solver1D_P2_Ref_T_EDP_a(0,4,h,0,1,.1,u,1,5,1);
    sll= sll.CLS(cl1,cl2);
    sll= sll.initialize(0);
    sll= sll.LoopTime();
    sll = sll.erreur();
   mrs = [sll.h mean(sll.er)*10^3];
  dlmwrite('spatial10^-3.xls',mrs,'delimiter','\t','-append','newline', 'pc');
 end


% ********************* Spataile-temporel Validaiton  **************************

header = [00000 hsps];
dlmwrite('tmix10^-3.xls',header,'delimiter','\t','precision',3)
disp(header)
for ht = htime
    mrs =[ht];
for h = hsps
    sll = Solver1D_P2_Ref_T_EDP_a(0,4,h,0,1,ht,u,1,5,1);
    sll= sll.CLS(cl1,cl2);
    sll= sll.initialize(0);
    sll= sll.LoopTime();
    sll = sll.erreur();
    mrs = [mrs mean(sll.mrs)*10^3];
end
    disp(mrs);
    dlmwrite('tmix10^-3.xls',mrs,'delimiter','\t','-append','newline', 'pc');

end


