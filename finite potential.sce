clf;clc;clear;  //Finite Potential well
h=1973;m=0.511e6;e=3.795;
a=2;n1=4;V1=0;
n=800;V0=-20;E=linspace(V0+0.1,V1+0.2,1000)
r=linspace(-n1*a,n1*a,n);
d=r(2)-r(1);
t=-(h^2)/(2*m*(d^2));
V=zeros(n,n);
K=-2*eye(n,n)+diag(ones(n-1,1),1) +diag(ones(n-1,1),-1);
for i=(0.5*n*(n1-1)/n1)+1:0.5*n*(n1+1)/n1
    V(i,i)=V0;
end
H=t*K+V;
[U1,EV]=spec(H);
En=diag(EV);
U1=-U1;
plotModes = sum(diag(En)>V0 & diag(En)<V1);
disp("Ground state energy : "+string(En(1))+" eV")
for i=2:plotModes
    disp(string(i-1)+" excited state energy : "+string(En(i))+" eV")
end
figure(0);scf(0);clf(0);
plot(r',diag(V)*ones(1,n))
title('Finite potential well')
xlabel('x (A)');ylabel('V (eV)')
for plotMode = 1:plotModes
    rescale=(EV(1,1)-V0)/(max(abs(U1(:,plotMode)))+min(abs(U1(:,plotMode))));
    U1_rescale(:,i)=U1(:,plotMode)*rescale+En(plotMode);
  plot(r', U1_rescale(:,i), 'r')
  plot(r', En(plotMode)*ones(n,1),'k--')
end
a1=gca()
a1.tight_limits=["on","on"];
a1.data_bounds=[-n1*a,V0-1;n1*a,1]
e=sqrt(-E./(E-V0*ones(E)));
theta=a*sqrt(2*m*(E-V0))/h;
figure(1);scf(1);clf(1);
plot(E,tan(theta),'r');plot(E,e,'k--');plot(E,-(ones(E)./tan(theta)),'b')
title("Plot of Trancendental equation")
xlabel('z');ylabel('z*tan(z)  ,z*cot(z)')
a2=gca();
a2.tight_limits=["on","on"];
a2.x_location="origin";
a2.data_bounds=[V0,-10;V1+0.5, 10]
