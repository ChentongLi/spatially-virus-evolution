% dt=0.05;
% this is the newton method to find the root.
% this method is worked for the model that the spectral radius go to 1 as the parameter go to infinity (which works bad for line search method). 
dx=0.1;
x=-10:dx:10;
Lm=length(x);
lambda = 100*x.^2+10000;
k1 = 0.000025*(0.01*x.^2 + 1);
k2 = 0.000065*(0.01*x.^2 + 1);
T0s = 1000 * x.^2 - 85091.8 * exp(-0.1 * x) - 85091.8 * exp(0.1 * x) + 300000;
% y(x) = 1000 * x^2 - 85091.8 * e^(-0.1 * x) - 85091.8 * e^(0.1 * x) + 300000;

D1 = 0.09648; D2 = 0.05; D3 = 0.08; D4 = 0.17; p = 0.3; dT = 0.01; dL = 0.001; dI = 0.2;
dV = 23; N = 2000; alpha= 0.01; alpha1 = 0.5; alpha2 = 0.4;

T0 = 10000; L0 = 0; I0 = 0; V0 = 100000;

Delta=diag(-2*ones(Lm,1)/dx^2)+diag(ones(Lm-1,1)/dx^2,1)+diag(ones(Lm-1,1)/dx^2,-1);
Delta(1,1)=-2/dx^2;Delta(1,2)=2/dx^2;
Delta(Lm,Lm-1)=2/dx^2;Delta(Lm,Lm)=-2/dx^2;

% for k=16:0.01:20 
%   k=1/k;
%   FV = [zeros(Lm,Lm)   ,  0.28*k*p* diag(k2.*T0s)+ D2*Delta-diag(ones(Lm,1)*(alpha+dL)), 0.4*k*p*diag(k1.*T0s);
%       diag(ones(Lm,1)*alpha*k), 0.28*k*(1-p)* diag(k2.*T0s)+ D3*Delta-diag(ones(Lm,1)*dI), 0.4*k*(1-p)*diag(k1.*T0s);
%       zeros(Lm,Lm)  ,   diag(ones(Lm,1)*N*k*0.3*dI), D4*Delta-diag(ones(Lm,1)*dV)] ;
%   le= max(eig(expm(FV)))
%   if abs(le - 1) <= 1e-2
%       R0=1/k;
%       break;
%   end
% end

k=20;
le=10;
while abs(le) > 1e-10
    FV = [zeros(Lm,Lm)   ,  0.28/k*p* diag(k2.*T0s)+ D2*Delta-diag(ones(Lm,1)*(alpha+dL)), 0.4/k*p*diag(k1.*T0s);
          diag(ones(Lm,1)*alpha/k), 0.28/k*(1-p)* diag(k2.*T0s)+ D3*Delta-diag(ones(Lm,1)*dI), 0.4/k*(1-p)*diag(k1.*T0s);
          zeros(Lm,Lm)  ,   diag(ones(Lm,1)*N/k*0.3*dI), D4*Delta-diag(ones(Lm,1)*dV)] ;
    le= max(eig(expm(FV)))-1;
    kk=k+0.1;
    FV = [zeros(Lm,Lm)   ,  0.28/kk*p* diag(k2.*T0s)+ D2*Delta-diag(ones(Lm,1)*(alpha+dL)), 0.4/kk*p*diag(k1.*T0s);
          diag(ones(Lm,1)*alpha/kk), 0.28/kk*(1-p)* diag(k2.*T0s)+ D3*Delta-diag(ones(Lm,1)*dI), 0.4/kk*(1-p)*diag(k1.*T0s);
          zeros(Lm,Lm)  ,   diag(ones(Lm,1)*N/kk*0.3*dI), D4*Delta-diag(ones(Lm,1)*dV)] ;    
    le1= max(eig(expm(FV)))-1;
    k=k-0.1*le/(le1-le);
end
R0=k


