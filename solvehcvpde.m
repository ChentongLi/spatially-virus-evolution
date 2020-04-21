function solvehcvpde
dx=0.1;L=7;
Lm=floor(L/dx);
x=0:dx:dx*(Lm-1);
dt=0.1;
N=20000;
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
diffusionv(1,1)=Dv/dx;diffusionv(1,2)=-Dv/dx;
diffusionv(Lm,Lm-1)=Dv/dx;diffusionv(Lm,Lm)=-Dv/dx;S=ones(Lm,1)*lambda/ds;
I=zeros(Lm,1);
v=zeros(Lm,1);
I(1)=1;
v(1)=10;
S0=S;I0=I;v0=v;
for k=1:N
    temp=(S+I)/Sm-1;
    S=(lambda+S/dt)./(1/dt+r1*temp+ds+beta*v);
    temp=(S+I)/Sm-1;
    I=(beta*S.*v+I/dt)./(1/dt+r2*temp+di);
    v=(eye(Lm)*(dv+1/dt)-diffusionv)\(alpha*I+v/dt);
    S0=[S0,S];
    I0=[I0,I];
    v0=[v0,v];
end

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

S=ones(Lm,1)*lambda/ds-10;
I=ones(Lm,1)*10;
v=ones(Lm,1)*10;
SA=S;IA=I;vA=v;
N2=15000;
for k=1:N2
    % refresh S,I
    temp=(S+I)/Sm-1;
    S=(lambda+S/dt)./(1/dt+r1*temp+ds+beta*v);
    temp=(S+I)/Sm-1;
    I=(beta*S.*v+I/dt)./(1/dt+r2*temp+di);

    DnA=DA/(exp(dA*TA)-1)*exp(-gA*x/Q)*exp(-dA*rem(k*dt,TA));
    DnA=DnA';
    PA=IC50A./(IC50A+DnA)*alpha;

    v=(eye(Lm)*(dv+1/dt)-diffusionv)\(PA.*I+v/dt);
    SA=[SA,S];
    IA=[IA,I];
    vA=[vA,v];
end

S=ones(Lm,1)*lambda/ds-10;
I=ones(Lm,1)*10;
v=ones(Lm,1)*10;
SD=S;ID=I;vD=v;
N3=10000;
for k=1:N3
    % refresh S,I
    temp=(S+I)/Sm-1;
    S=(lambda+S/dt)./(1/dt+r1*temp+ds+beta*v);
    temp=(S+I)/Sm-1;
    I=(beta*S.*v+I/dt)./(1/dt+r2*temp+di);

    DnD=DD/(exp(dD*TD)-1)*exp(-gD*x/Q)*exp(-dD*rem(k*dt,TD));
    DnD=DnD';
    PD=IC50D./(IC50D+DnD)*alpha;

    v=(eye(Lm)*(dv+1/dt)-diffusionv)\(PD.*I+v/dt);
    SD=[SD,S];
    ID=[ID,I];
    vD=[vD,v];
end
t=0:dt:N*dt;
t2=0:dt:N2*dt;
t3=0:dt:N3*dt;

subplot(3,3,1)
surf(t,x,S0)
title('(a)')
xlabel('Time(hour)')
ylabel('Positions(cm)')
zlabel('S(x,t)')
ylim([0,7])
shading interp
subplot(3,3,2)
surf(t,x,I0)
title('(b)')
xlabel('Time(hour)')
ylabel('Positions(cm)')
zlabel('I(x,t)')
ylim([0,7])
shading interp
subplot(3,3,3)
surf(t,x,v0)
title('(c)')
xlabel('Time(hour)')
ylabel('Positions(cm)')
ylim([0,7])
zlabel('v(x,t)')
shading interp
subplot(3,3,4)
surf(t2,x,SA)
title('(d)')
xlabel('Time(hour)')
ylabel('Positions(cm)')
zlabel('S(x,t)')
ylim([0,7])
shading interp
subplot(3,3,5)
surf(t2,x,IA)
title('(e)')
xlabel('Time(hour)')
ylabel('Positions(cm)')
zlabel('I(x,t)')
ylim([0,7])
shading interp
subplot(3,3,6)
surf(t2,x,vA)
title('(f)')
xlabel('Time(hour)')
ylabel('Positions(cm)')
xlim([0,150])
zlabel('v(x,t)')
ylim([0,7])
shading interp
subplot(3,3,7)
surf(t3,x,SD)
title('(g)')
xlabel('Time(hour)')
ylabel('Positions(cm)')
zlabel('S(x,t)')
ylim([0,7])
shading interp
subplot(3,3,8)
surf(t3,x,ID)
title('(h)')
xlabel('Time(hour)')
ylabel('Positions(cm)')
zlabel('I(x,t)')
ylim([0,7])
shading interp
subplot(3,3,9)
surf(t3,x,vD)
title('(i)')
xlabel('Time(hour)')
ylabel('Positions(cm)')
zlabel('v(x,t)')
xlim([0,150])
ylim([0,7])
shading interp
