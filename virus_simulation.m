function virus_simulation
t=[];
y=[];
N=1000;
for i=1:N
    [tr,yr]=simulations();
    t=[t,tr];
    y=[y,yr];
end
subplot(3,2,1)
ytmp=y(1,:);
ytmp(ytmp==0)=[];
histfit(ytmp,20,'kernel')
title('(a)')
xlabel('The first resist drug one mutant virus appear positions on blood vessel ')
ylabel('Frequency')
subplot(3,2,2)
ytmp=y(2,:);
ytmp(ytmp==0)=[];
histfit(ytmp,20,'kernel')
title('(b)')
xlabel('The first resist drug two mutant virus appear positions on blood vessel ')
ylabel('Frequency')
subplot(3,2,3)
ytmp=y(3,:);
ytmp(ytmp==0)=[];
histfit(ytmp,20,'kernel')
title('(c)')
xlabel('The first resist two drugs mutant virus appear positions on blood vessel ')
ylabel('Frequency')
subplot(3,2,4)
ttmp=t(1,:);
ttmp(ttmp==0)=[];
histfit(ttmp,20,'kernel')
title('(d)')
xlabel('The first resist drug one mutant virus appear time')
ylabel('Frequency')
subplot(3,2,5)
ttmp=t(2,:);
ttmp(ttmp==0)=[];
histfit(ttmp,20,'kernel')
title('(e)')
xlabel('The first resist drug two mutant virus appear time')
ylabel('Frequency')
subplot(3,2,6)
ttmp=t(3,:);
ttmp(ttmp==0)=[];
histfit(ttmp,20,'kernel')
title('(f)')
xlabel('The first resist two drugs mutant virus appear time')
ylabel('Frequency')

function [tr,yr]=simulations()
tr=zeros(3,1);
yr=zeros(3,1);
u1=5e-4;
u2=6e-4;
dx=0.1;L=10;dt=0.1;
Lm=floor(L/dx);
x=0:dx:dx*(Lm-1);
D=0.02;Dv=0.1;
lambda=(sin(x)+1)*2;lambda=lambda';% birth distribution
ds=0.01;di=0.5;dv=0.8;
a0=40;decay=0.15;
beta=0.001*ones(4,1);
q=0.05;
diffusion=diag(ones(Lm,1)*(-2*D/dx^2))+diag(ones(Lm-1,1)*(D/dx^2-q/(2*dx)),1)+diag(ones(Lm-1,1)*(D/dx^2+q/(2*dx)),-1);
diffusion(1,Lm)=D/dx^2+q/(2*dx);diffusion(Lm,1)=D/dx^2-q/(2*dx);
diffusionv=diag(ones(Lm,1)*(-2*Dv/dx^2))+diag(ones(Lm-1,1)*(Dv/dx^2-q/(2*dx)),1)+diag(ones(Lm-1,1)*(Dv/dx^2+q/(2*dx)),-1);
diffusionv(1,Lm)=Dv/dx^2+q/(2*dx);diffusionv(Lm,1)=Dv/dx^2-q/(2*dx);
B=diag(ones(Lm,1)*di)-diffusion;
S=zeros(Lm,1);
I=zeros(Lm,4);
v=zeros(Lm,4);
v(1,1)=20;
alpha=zeros(Lm,4);
AD=[0.95,0.9];
C0=1;
IC50=[0.1,0.1];
alpha(:,1)=a0*(1-AD(1)*C0*exp(-decay*x)./(IC50(1)+C0*exp(-decay*x))).*(1-AD(2)*C0*exp(-decay*x)./(IC50(2)+C0*exp(-decay*x)));
alpha(:,2)=a0*(1-AD(1)*C0*exp(-decay*x)./(IC50(1)+C0*exp(-decay*x)));
alpha(:,3)=a0*(1-AD(2)*C0*exp(-decay*x)./(IC50(2)+C0*exp(-decay*x)));
alpha(:,4)=a0;
for k=1:2000
    % refresh S,I
    A=diag(v*beta+ones(Lm,1)*ds)-diffusion;
    S=linsolve(A,lambda);
    for i=1:4
        I(:,i)=linsolve(B,beta(i)*v(:,i).*S);
    end
    for i=1:2
        % reaction 
        newv=poissrnd(alpha.*I*di*dt/2);
        [tr,yr,newv]=mutation(newv,Lm,dx,u1,u2,(k-1)*dt+i*dt/2,tr,yr);
        v=v+newv;
        v=dey(v,Lm,dv,dt/2);
        % diffusion
        v=(eye(Lm)+diffusionv*dt/2)*v;
    end
end

function v=dey(v,Lm,dv,dt)
for i=1:Lm
    for j=1:4;
        tmp=floor(v(i,j));
        for m=1:tmp
            if rand()<1-exp(-dv*v(i,j)*dt) 
                v(i,j)=v(i,j)-1;
            end
        end
        dp=v(i,j)-tmp;
        if rand()<dp*(1-exp(-dv*v(i,j)*dt))
            v(i,j)=v(i,j)-dp;
        end
    end
end

function [tr,yr,newv]=mutation(newv,Lm,dx,u1,u2,t,tr,yr)
for i=1:Lm
    tmp=floor(newv(i,1));
    for m=1:tmp
        tmpr=rand();
        if tmpr<=u1
            newv(i,1)=newv(i,1)-1;
            newv(i,2)=newv(i,2)+1;
            if yr(1)==0
                yr(1)=i*dx;
                tr(1)=t;
            end
        end  
        if tmpr>u1 && tmpr<u1+u2
            newv(i,1)=newv(i,1)-1;
            newv(i,3)=newv(i,3)+1;
            if yr(2)==0
                yr(2)=i*dx;
                tr(2)=t;
            end
        end
    end
    tmp=floor(newv(i,2));
    for m=1:tmp
        if rand()<=u2
            newv(i,2)=newv(i,2)-1;
            newv(i,4)=newv(i,4)+1;
            if yr(3)==0
                yr(3)=i*dx;
                tr(3)=t;
            end
        end
    end
    tmp=floor(newv(i,3));
    for m=1:tmp
        if rand()<=u1
            newv(i,3)=newv(i,3)-1;
            newv(i,4)=newv(i,4)+1;
            if yr(3)==0
                yr(3)=i*dx;
                tr(3)=t;
            end
        end
    end
end
