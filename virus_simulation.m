%function [xxA,xxD]= virus2()
function virus_simulation
dx=0.1;L=7;
Lm=floor(L/dx);

lambda=1.02*1e4/24;
ds=1.2*1e-3/24;
Sm=lambda/ds;

S=ones(Lm,1)*Sm*0.5;
I=[ones(Lm,1)*Sm*0.7,ones(Lm,1)*Sm*0.2,ones(Lm,1)*Sm*0.15,zeros(Lm,1)];
v=[ones(Lm,1)*Sm*0.3,ones(Lm,1)*Sm*0.1,ones(Lm,1)*Sm*0.1,zeros(Lm,1)];

xxA=[];
xxD=[];
N=200;

for i=1:N
    [xrA,xrD]=simuhcv(S,I,v);
    xxA=[xxA,xrA];
    xxD=[xxD,xrD];
end
subplot(2,1,1)
histfit(xxA,25,'kernel')
title('(a)')
xlabel('Resist Asunaprevir mutant virus appear positions on liver')
ylabel('Frequency')
xlim([0,7])
subplot(2,1,2)
histfit(xxD,25,'kernel')
title('(b)')
xlabel('Resist Daclatasvir mutant virus appear positions on liver ')
ylabel('Frequency')
xlim([0,7])
print pdis1.eps -depsc2 -r600

function [xrA,xrD]=simuhcv(S0,I0,v0)
xrA=[];xrD=[];
dx=0.1;L=7;
Lm=floor(L/dx);
x=0:dx:dx*(Lm-1);
dt=0.1;

u1=2.51e-5; u2=u1; % mutations per nucleotide per genome replication
lambda=1.02*1e4/24;
ds=1.2*1e-3/24;
di=1.5/24;
dv=11.4/24;
beta=2.02*1e-7/24;
alpha=22.55/24;
Dv=15.8*1e-8*60*60;
r1=1.55/24;
r2=5.5/24;
Sm=lambda/ds;
diffusionv=diag(ones(Lm,1)*(-2*Dv/dx^2))+diag(ones(Lm-1,1)*(Dv/dx^2),1)+diag(ones(Lm-1,1)*(Dv/dx^2),-1);
diffusionv(1,1)=-Dv/dx^2;diffusionv(Lm,Lm)=-Dv/dx^2;

rho=2e3;
clintA=1.83*1e-2;
Q=1.2;TA=12;fuA=0.01;
clA=49.5;vdA=194;
DA=200*1e3;
IC50A=2.45;

clintD=4.34*1e-3;
TD=24;fuD=0.006;
clD=4.2;vdD=47;
DD=50*1e3;
IC50D=0.041;

gA=clintA*fuA/L;
dA=clA/vdA;
DA=DA/vdA;

gD=clintD*fuD/L;
dD=clD/vdD;
DD=DD/vdD;
% v1=[];v2=[];v3=[];v4=[];
S=S0;
I=I0;
v=v0;
N=10000;
for k=1:N
    % refresh S,I
    temp=(S+sum(I,2))/Sm-1;
    S=(lambda+S/dt)./(1/dt+r1*temp+ds+beta*sum(v,2));
    temp=(S+sum(I,2))/Sm-1;
    for i=1:4
        I(:,i)=(beta*S.*v(:,i)+I(:,i)/dt)./(1/dt+r2*temp+di);
    end

    DnA=DA/(exp(dA*TA)-1)*exp(-gA*x/Q)*exp(-dA*rem(k*dt,TA));
    DnA=DnA';
    PA=IC50A./(IC50A+DnA);
    PAR=IC50A*rho./(IC50A*rho+DnA);

    DnD=DD/(exp(dD*TD)-1)*exp(-gD*x/Q)*exp(-dD*rem(k*dt,TD));
    DnD=DnD';
    PD=IC50D./(IC50D+DnD);
    PDR=IC50D*rho./(IC50D*rho+DnD);

    P=[PA.*PD,PA.*PDR,PD.*PAR,PAR.*PDR]*alpha;

    for i=1:2
        % reaction 
        newv=poissrnd(P.*I*dt/2);
        [newv,xA,xD]=mutation(newv,Lm,u1,u2,dx);
        xrA=[xrA,xA];
        xrD=[xrD,xD];
        v=v+newv;
        v=dey(v,Lm,dv,dt/2);
        % diffusion
        v=(eye(Lm)+diffusionv*dt/2)*v;
    end
%    v1(:,k)=v(:,1);
%    v2(:,k)=v(:,2);
%    v3(:,k)=v(:,3);
%    v4(:,k)=v(:,4);
end

function v=dey(v,Lm,dv,dt)
for i=1:Lm
    for j=1:4
        tmp=floor(v(i,j));
        dp=v(i,j)-tmp;
        if rand()<dp*(1-exp(-dv*v(i,j)*dt))
            v(i,j)=v(i,j)-dp;
        end
        for m=1:tmp
            if rand()<1-exp(-dv*v(i,j)*dt) 
                v(i,j)=v(i,j)-1;
            end
        end
    end
end

function [newv,xA,xD]=mutation(newv,Lm,u1,u2,dx)
xA=[];xD=[];
for i=1:Lm
    tmp=floor(newv(i,1));
    for m=1:tmp
        tmpr=rand();
        if tmpr<=u1
            newv(i,1)=newv(i,1)-1;
            newv(i,2)=newv(i,2)+1;
	    xD=[xD,dx*i];
        end  
        if tmpr>u1 && tmpr<=u1+u2
            newv(i,1)=newv(i,1)-1;
            newv(i,3)=newv(i,3)+1;
	    xA=[xA,dx*i];
        end
    end
    tmp=floor(newv(i,2));
    for m=1:tmp
        if rand()<=u2
            newv(i,2)=newv(i,2)-1;
            newv(i,4)=newv(i,4)+1;
	    xA=[xA,dx*i];
        end
    end
    tmp=floor(newv(i,3));
    for m=1:tmp
        if rand()<=u1
            newv(i,3)=newv(i,3)-1;
            newv(i,4)=newv(i,4)+1;
	    xD=[xD,dx*i];
        end
    end
end
