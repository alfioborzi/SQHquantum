%
% SQH Method for optimal control of 2 uncoupled spin 1/2 
% Model y' = (A + u B ) y, y(0)=y_0
% Cost functional with L2, L1, and discontinuous L1 costs
% the weights of these costs are nu, gamma, and beta (threshold s)
% Target y_d
%
% Main references:
%
% T. Breitenbach, A. Borzì,
% A sequential quadratic Hamiltonian scheme for solving non-smooth 
% quantum control problems with sparsity,
% Journal of Computational and Applied Mathematics, 369 (2019) 112583.
%
% A. Borzì, The Sequential Quadratic Hamiltonian Method, Chapman & Hall, 2023
%
disp(' ');
close all; 
clear all;

ltype = {'b-','r--','m-','b-*','r:','m-.'};

OCP=struct('dt',10^-5,'T',0.008,'umin',-60,'umax',60,'nu',10^-3,'beta',1,'gamma',0.5,'s',10);

% Initial condition and system matrices
y0=[0;0;1;0;0;1];
y0=y0/norm(y0);
A=[0,-1,0,0,0,0;1,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,1,0;0,0,0,-1,0,0;0,0,0,0,0,0]*2*pi*483;
B=[0,0,0,0,0,0;0,0,-1,0,0,0;0,1,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,-1;0,0,0,0,1,0]*2*pi;

%Algorithm SQH paramters
kappa=10^-8;                
eta=10^-5;
zeta=0.9;
sigma=1.1;
kmax=10000;        %Maximum iteration number
epsilon=1.0;      %Initial guess for epsilon

% target
yd=[0;1;0;0;1;0];
yd=yd/norm(yd);


beta=OCP.beta;
gamma=OCP.gamma;
s=OCP.s;

Nt=round(OCP.T/OCP.dt);

u=ones(1,Nt+1);
u_old=u;

y=forward(A,B,y0,u,OCP);
y_old=y;

p=backward(A,B,y(:,end),yd,u,OCP);

J_plot(1,1)=get_J(y(:,end),u,yd,OCP);
J_plot(1,1)
Jk(1)=J_plot(1,1);



count_updates=1;
ctime=cputime;                  %Start time measurement

% SQH loop 

for k=1:kmax
    
    epsik(k) = epsilon ; 
    
    for i=1:Nt+1
                
     uz(1)=min(max(s,(2*epsilon*u_old(1,i)-(p(:,i)')*B*y(:,i)-beta-gamma) ... 
         /(2*epsilon+OCP.nu)),OCP.umax);
     uz(2)=min(max(OCP.umin,(2*epsilon*u_old(1,i)-(p(:,i)')*B*y(:,i)+beta+gamma) ... 
         /(2*epsilon+OCP.nu)),-s);

     uz(3)=min(max(-s,(2*epsilon*u_old(1,i)-(p(:,i)')*B*y(:,i) +gamma) ... 
         /(2*epsilon+OCP.nu)),0);

     uz(4)=min(max(0,(2*epsilon*u_old(1,i)-(p(:,i)')*B*y(:,i) -gamma) ... 
         /(2*epsilon+OCP.nu)),s);


        H(1)=(OCP.nu/2)*uz(1)^2 + OCP.beta*(abs(uz(1))>s).*abs(uz(1)) ... 
            + OCP.gamma*abs(uz(1)) ...
            + (p(:,i)')*(A+uz(1)*B)*y(:,i) + epsilon*(uz(1)-u_old(1,i))^2;
        H(2)=(OCP.nu/2)*uz(2)^2 + OCP.beta*(abs(uz(2))>s).*abs(uz(2)) ... 
            + OCP.gamma*abs(uz(2)) ...
            + (p(:,i)')*(A+uz(2)*B)*y(:,i) + epsilon*(uz(2)-u_old(1,i))^2;

        H(3)=(OCP.nu/2)*uz(3)^2                                       ... 
            + OCP.gamma*abs(uz(3)) ...
            + (p(:,i)')*(A+uz(3)*B)*y(:,i) + epsilon*(uz(3)-u_old(1,i))^2;

        H(4)=(OCP.nu/2)*uz(4)^2                                       ... 
            + OCP.gamma*abs(uz(4)) ...
            + (p(:,i)')*(A+uz(4)*B)*y(:,i) + epsilon*(uz(4)-u_old(1,i))^2;


        [~,pos]=min(H);
        u(1,i)=uz(pos);

    end 
  
   du=sum((u-u_old).^2)*OCP.dt;
   
   y=forward(A,B,y0,u,OCP);
   J_int=get_J(y(:,end),u,yd,OCP);
   if(J_int-Jk(count_updates)>-eta*du)       
        % If the improvemnt to the cost functional is not sufficiently large, 
        % take the old values of y and u, increase epsilon
       u=u_old;
       y=y_old;
       epsilon=epsilon*sigma;                
   else
        % If the improvment to the cost functional is sufficently large, 
        % take the new values of y, u, J, calculate p and decrease epsilon
       count_updates=count_updates+1;
       p=backward(A,B,y(:,end),yd,u,OCP);
       epsilon=epsilon*zeta;                         
       Jk(count_updates)=J_int;
       
       u_old=u;
       y_old=y;
   % Print current values
        fprintf('k %i | k up %i | J %e | eps %e | tau %e\n',k,... 
            count_updates-1,J_int,epsilon,du) 

   end
   if(du<kappa && k > 1 )
        % If the norm square of the update is too small, then return the latest 
        % variable u that caused sufficient decrease of the cost functional 
        fprintf('Convergence achieved\n');
        u=u_old;
        break;
   end 
end 



ct=cputime-ctime; %Stop time measurment for the calculation 
fprintf('CPU-time in sec: %f\n',ct)


figure(1)
plot(0:OCP.dt/(OCP.T):1,u/100,ltype{1},'Linewidth',2)
axis([-inf inf -1 1])
xlabel('t/T')
ylabel('u/100')
  print('-depsc2', 'qmdu01.eps','-b0'); 
% print('-dpdf', 'qmdu01.pdf','-b0');



figure(2)
plot(0:OCP.dt/(OCP.T):1,y(1,:),ltype{1},'Linewidth',2), hold on; 
plot(0:OCP.dt/(OCP.T):1,y(2,:),ltype{2},'Linewidth',2),
plot(0:OCP.dt/(OCP.T):1,y(3,:),ltype{3},'Linewidth',2),
%plot(0:OCP.dt/(OCP.T):1,y(4,:),ltype{4},'Linewidth',1),
%plot(0:OCP.dt/(OCP.T):1,y(5,:),ltype{5},'Linewidth',1),
%plot(0:OCP.dt/(OCP.T):1,y(6,:),ltype{6},'Linewidth',1),
axis([-inf inf -1 1])
xlabel('x/T')
ylabel('y')
legend({'$y_1$','$y_2$','$y_3$'},'Interpreter','Latex','FontSize',12)

print('-depsc2', 'trajectory123qm.eps','-b0'); 
% print('-dpdf', 'trajectory123qm.pdf','-b0');

figure(3)
plot(0:OCP.dt/(OCP.T):1,y(4,:),ltype{1},'Linewidth',2), hold on; 
plot(0:OCP.dt/(OCP.T):1,y(5,:),ltype{2},'Linewidth',2),
plot(0:OCP.dt/(OCP.T):1,y(6,:),ltype{3},'Linewidth',2),
%plot(0:OCP.dt/(OCP.T):1,y(4,:),ltype{4},'Linewidth',1),
%plot(0:OCP.dt/(OCP.T):1,y(5,:),ltype{5},'Linewidth',1),
%plot(0:OCP.dt/(OCP.T):1,y(6,:),ltype{6},'Linewidth',1),
axis([-inf inf -1 1]);
xlabel('x/T')
ylabel('y')
legend({'$y_4$','$y_5$','$y_6$'},'Interpreter','Latex','FontSize',12)

print('-depsc2', 'trajectory456qm.eps','-b0'); 
% print('-dpdf', 'trajectory456qm.pdf','-b0');


figure(4)
plot(Jk,ltype{4},'Linewidth',2)
xlabel('SQH iterations')
ylabel('J')

print('-depsc2', 'histJqmd01.eps','-b0'); 
% print('-dpdf', 'histJqmd01.pdf','-b0');

figure(5)
semilogy(epsik,ltype{1},'Linewidth',2);
xlabel('SQH iterations')
ylabel('epsilon')

print('-depsc2', 'epsqmd01.eps','-b0'); 
% print('-dpdf', 'epsqmd01.pdf','-b0');


Ce=((y(:,end)-yd)'*(y(:,end)-yd));
fprintf('Eucl distance to target: %e\n',Ce) 

% Check solution for optimality;
fprintf('\n')
fprintf('Check numerical optimality of the SQH solution\n')
numeric_optimality(u,y0,yd,A,B,OCP)


