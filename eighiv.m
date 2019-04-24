function eighiv
record=[];
t=5:0.1:35;
for i=t
    record=[record,eigvirus(i)];
end
plot(t,record,'LineWidth',5)
xlabel('blood vessel length')
ylabel('R_0')

function R0=eigvirus(yy)
% dx=0.1;L=10;C0=1;
% Lm=floor(L/dx);
% x=0:dx:dx*(Lm-1);
% D=0.02;Dv=0.1;
% AD=[0.9,0.85];
% IC50=[0.1,0.1];
% lambda=(sin(x)+1)*2;lambda=lambda';% birth distribution
% ds=0.01;di=0.5;dv=0.8;beta=0.001;q=0.1;
% a0=15;decay=0.1;

dx=0.1;L=yy;C0=1;
Lm=floor(L/dx);
x=0:dx:dx*(Lm-1);
D=0.02;Dv=0.1;
AD=[0.9,0.85];
IC50=[0.1,0.1];
lambda=(sin(x)+1)*2;lambda=lambda';% birth distribution
ds=0.01;di=0.5;dv=0.8;beta=0.001;q=0.05;
a0=15;decay=0.1;


alpha=zeros(Lm,4);
alpha(:,1)=a0*(1-AD(1)*C0*exp(-decay*x)./(IC50(1)+C0*exp(-decay*x))).*(1-AD(2)*C0*exp(-decay*x)./(IC50(2)+C0*exp(-decay*x)));
alpha(:,2)=a0*(1-AD(1)*C0*exp(-decay*x)./(IC50(1)+C0*exp(-decay*x)));
alpha(:,3)=a0*(1-AD(2)*C0*exp(-decay*x)./(IC50(2)+C0*exp(-decay*x)));
alpha(:,4)=a0;

diffusion=diag(ones(Lm,1)*(-2*D/dx^2))+diag(ones(Lm-1,1)*(D/dx^2-q/(2*dx)),1)+diag(ones(Lm-1,1)*(D/dx^2+q/(2*dx)),-1);
diffusion(1,Lm)=D/dx^2+q/(2*dx);diffusion(Lm,1)=D/dx^2-q/(2*dx);
diffusionv=diag(ones(Lm,1)*(-2*Dv/dx^2))+diag(ones(Lm-1,1)*(Dv/dx^2-q/(2*dx)),1)+diag(ones(Lm-1,1)*(Dv/dx^2+q/(2*dx)),-1);
diffusionv(1,Lm)=Dv/dx^2+q/(2*dx);diffusionv(Lm,1)=Dv/dx^2-q/(2*dx);

A=diag(ones(Lm,1)*ds)-diffusion;
S0=linsolve(A,lambda);

F=[zeros(Lm,Lm),diag(S0)*beta;zeros(Lm,Lm),zeros(Lm,Lm)];
V=[di*eye(Lm),zeros(Lm,Lm);-di*diag(alpha(:,2)),dv*eye(Lm)];
K=[-diffusion,zeros(Lm,Lm);zeros(Lm,Lm),-diffusionv];
%plot(S0)
R0=max(eig(F/(K+V)));


