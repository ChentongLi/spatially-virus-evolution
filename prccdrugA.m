function prccdrugA()

p=[];
% by this loop, we can get a more accurate result. But this code need serveral days to run.
for i=1:100
    p=[prcc();p];
end
% save prcc.mat p

bar(mean(p))
% if just want to get the quickly similar result as we paper showed, using bar(prcc()) is OK.

title('Sensitivity of R_0 of with Asunaprevir model with respect to model parameters')
set(gca,'xTick',1:7,'xTicklabel',{'CL_{int}','Q','T','f_u','CL','D','IC_{50}'})

function PRCC=prcc()

runs=1000;

clint_lhs=LHS_Call(1.83*0.8,1.83,1.83*1.2,0.04,runs,'unif')*1e-2;
Q_lhs=LHS_Call(1.2*0.8,1.2,1.2*1.2,0.04,runs,'unif');
T_lhs=LHS_Call(12*0.8,12,12*1.2,0.04,runs,'unif');
fu_lhs=LHS_Call(1*0.8,1,1*1.2,0.04,runs,'unif')*1e-2;
cl_lhs=LHS_Call(49.5*0.8,49.5,49.5*1.2,0.4,runs,'unif');
D_lhs=LHS_Call(2*0.8,2,2*1.2,0.04,runs,'unif')*1e5;
IC50_lhs=LHS_Call(2.45*0.8,2.45,2.45*1.2,0.04,runs,'unif');

y=[];
for i=1:runs
    y(i)=eigdrug(clint_lhs(i),Q_lhs(i),T_lhs(i),fu_lhs(i),cl_lhs(i),D_lhs(i),IC50_lhs(i));
end

R0=y';
PM=[clint_lhs,Q_lhs,T_lhs,fu_lhs,cl_lhs,D_lhs,IC50_lhs];
rank_PM=ranking(PM);
rank_R0=ranking(R0);
N=size(PM,2);

PRCC=[];
for j=1:N
    PMt=rank_PM;
    rank_x=PMt(:,j);
    PMt(:,j)=[];
    s1=regstats(rank_R0,PMt,'linear','r');
    y_res=s1.r;
    s2=regstats(rank_x,PMt,'linear','r');
    x_res=s2.r;
    PRCC=[PRCC,corr(y_res,x_res)];
end

% bar(PRCC)
% title('Sensitivity of R_0 of with Asunaprevir model with respect to model parameters')
% set(gca,'xTick',1:7,'xTicklabel',{'CL_{int}','Q','T','f_u','CL','D','IC_{50}'})

function r=ranking(x)

[a,b]=size(x);
for j=1:b
   [s,i]=sort(x(:,j));
   r(i,j)=[1:a]';
end

function R0=eigdrug(clint,Q,T,fu,cl,D,IC50)

dx=0.1;L=7;
Lm=floor(L/dx);
x=0:dx:dx*(Lm-1);

lambda=1.02*1e4/24;
ds=1.2*1e-3/24;
di=1.5/24;
dv=11.4/24;
beta=2.02*1e-7/24;
alpha=22.55/24;
Dv=15.8*1e-8*60*60;

% clint=1.83*1e-2;
% Q=1.2;
% T=12;
% fu=0.01;
% cl=49.5;
vd=194;
% D=200*1e3;
% IC50=2.45;

g=clint*fu/L;
d=cl/vd;
A=D/vd*exp(-g*x/Q)/(exp(d*T)-1);
p=alpha/d*log((IC50*exp(d*T)+A)./(IC50+A)); % result of \int_0^T p(x,t) dt

diffusionv=diag(ones(Lm,1)*(-2*Dv/dx^2))+diag(ones(Lm-1,1)*(Dv/dx^2),1)+diag(ones(Lm-1,1)*(Dv/dx^2),-1);
diffusionv(1,1)=-2*Dv/dx^2;diffusionv(1,2)=2*Dv/dx^2;
diffusionv(Lm,Lm-1)=2*Dv/dx^2;diffusionv(Lm,Lm)=-2*Dv/dx^2;
S0=lambda/ds;

F=[zeros(Lm,Lm),diag(S0*beta*ones(Lm,1));zeros(Lm,Lm),zeros(Lm,Lm)];
V=[di*eye(Lm),zeros(Lm,Lm);zeros(Lm,Lm),dv*eye(Lm)-diffusionv];
B=-V*T+[zeros(Lm,Lm),zeros(Lm,Lm);diag(p),zeros(Lm,Lm)];

for k=0.01:0.01:20
   A=expm(F*T/k+B);
   le=max(eig(A));
   if abs(le-1)<1e-1
        R0=k;
        break;
   end
end

% if k==20
%     [clint,Q,T,fu,cl,IC50]
%     D
%     le
%     A=expm(F*T/0.01+B);
%     max(eig(A))
% end


