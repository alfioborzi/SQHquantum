function [J] = get_J(y,u,yd,OCP)
nu=OCP.nu;
dt=OCP.dt;
beta=OCP.beta;
gamma=OCP.gamma;
s=OCP.s;

J=0.5*(norm(y-yd)^2) + (nu/2)*u*(u')*dt + gamma*sum(abs(u))*dt + beta*sum((abs(u)>s).*abs(u))*dt;
end

