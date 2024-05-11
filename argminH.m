function [ v ] = argminH(y,u,p,A,B,epsilon,OCP )
beta=OCP.beta;
gamma=OCP.gamma;
s=OCP.s;

u(1)=min(max((-p'*B*y-beta-gamma)/OCP.nu,s),OCP.umax);
u(2)=min(max((-p'*B*y+beta+gamma)/OCP.nu,OCP.umin),-s);
u(3)=min(max((-p'*B*y)/OCP.nu+gamma,-s),0);
u(4)=min(max((-p'*B*y)/OCP.nu-gamma,0),s);


H(1)=OCP.nu*0.5*u(1)^2 + OCP.gamma*abs(u(1)) ... 
    + OCP.beta*(abs(u(1))>s).*abs(u(1))+p'*B*y*u(1);
H(2)=OCP.nu*0.5*u(2)^2 + OCP.gamma*abs(u(2))... 
    + OCP.beta*(abs(u(2))>s).*abs(u(2))+p'*B*y*u(2);
H(3)=OCP.nu*0.5*u(3)^2 +p'*B*y*u(3) + OCP.gamma*abs(u(3));
H(4)=OCP.nu*0.5*u(4)^2 +p'*B*y*u(4) + OCP.gamma*abs(u(4));
[~,pos]=min(H);
v=u(pos);


end


