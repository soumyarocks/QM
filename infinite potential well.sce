//Name - Soumya Swaroop       //roll no.- 21HPH2657


clf;clc;clear;  //Infinite Potential well
h=1973;m=0.511e6;e=3.795;
a=2;n1=1.5;
n=1000;V0=0;V1=9.999999999e8;
r=linspace(-n1*a,n1*a,n);E=linspace(V0+0.1,V1+0.2,1000)
d=r(2)-r(1);
t=-(h^2)/(2*m*(d^2));
V=V1*eye(n,n);
K=-2*eye(n,n)+diag(ones(n-1,1),1) +diag(ones(n-1,1),-1);
for i=(0.5*n*(n1-1)/n1)+1:0.5*n*(n1+1)/n1
    V(i,i)=V0;
end
H=t*K+V;
[U1,EV]=spec(H);
En=diag(EV);
disp("Ground state energy : "+string(En(2))+" eV")
disp("First excited state energy : "+string(En(3))+" eV")
disp("Second excited state energy : "+string(En(4))+" eV")
for i=1:4
subplot(2,2,i)
plot(r',-U1(:,i),'r')
title("Plot of wavefunction "+string(i))
xlabel("x(A)")
ylabel("Wavefunction")
a=gca();
a.x_location="origin";
a.y_location="origin";
end
