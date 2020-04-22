function prcc

runs=1000;

%beta_LHS=lhsnorm(2.02,1,runs)*1e-7/24;
%alpha_LHS=lhsnorm(22.55,10,runs)/24;
beta_LHS=LHS_Call(0.8*2.02,2.02,1.2*2.02,0.04,runs,'unif')*1e-7/24;
alpha_LHS=LHS_Call(0.8*22.55,22.55,1.2*22.55,0.04,runs,'unif')/24;
Dv_LHS=LHS_Call(0.8*15.8,15.8,1.2*15.8,0.04,runs,'unif')*1e-8*3600;
di_LHS=LHS_Call(0.8*1.5,1.5,1.2*1.5,0.04,runs,'unif')/24;
dv_LHS=LHS_Call(0.8*11.4,11.4,11.4*1.2,0.04,runs,'unif')/24;

y=[];
for i=1:runs
    y(i)=eigvirus(beta_LHS(i),alpha_LHS(i),Dv_LHS(i),di_LHS(i),dv_LHS(i));
end

R0=y';
PM=[beta_LHS,alpha_LHS,Dv_LHS,di_LHS,dv_LHS];
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

bar(PRCC)
ylabel('Sensitivity indexes of R_0')
title('Sensitivity of R_0 of without drug model with respect to model parameters')
ylim([-1,1])
set(gca,'xTick',1:5,'xTicklabel',{'\beta','\alpha','D_v','d_I','d_v'})
print prcc1.eps -depsc2 -r600

function r=ranking(x)

[a,b]=size(x);
for j=1:b
   [s,i]=sort(x(:,j));
   r(i,j)=[1:a]';
end

function R0=eigvirus(beta,alpha,Dv,di,dv)
dx=0.1;L=7;
Lm=floor(L/dx);
x=0:dx:dx*(Lm-1);

lambda=1.02*1e4/24;
ds=1.2*1e-3/24;
% di=1.5/24;
% dv=11.4/24;
% beta=2.02*1e-7/24;
% alpha=22.55/24;
% Dv=15.8*1e-8*60*60;

diffusionv=diag(ones(Lm,1)*(-2*Dv/dx^2))+diag(ones(Lm-1,1)*(Dv/dx^2),1)+diag(ones(Lm-1,1)*(Dv/dx^2),-1);
diffusionv(1,1)=-2*Dv/dx^2;diffusionv(1,2)=2*Dv/dx^2;
diffusionv(Lm,Lm-1)=2*Dv/dx^2;diffusionv(Lm,Lm)=-2*Dv/dx^2;
S0=lambda/ds;

F=[zeros(Lm,Lm),diag(S0*beta*ones(Lm,1));zeros(Lm,Lm),zeros(Lm,Lm)];
V=[di*eye(Lm),zeros(Lm,Lm);-diag(ones(Lm,1)*alpha),dv*eye(Lm)-diffusionv];

R0=max(eig(F/V));

