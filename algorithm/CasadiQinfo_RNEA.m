function [Hx, Hu,Qxx, Quu, Qux] = CasadiQinfo_RNEA(xi,ui, Vxi, Vxxi,params)

% Some definitions    
x_size = length(xi);
q =xi(1:length(xi)/2); 
qd = xi(length(xi)/2+1:end);

ui_un = ui;
ui = [ui;0];

lambda = Vxi;
lambda1 = Vxi(1:length(xi)/2); %for simplification sake
lambda2 = Vxi(length(xi)/2+1:end);
sz = params.sz; Nb= sz(1);
sztau = sz; sztau(end)= sz(end)-1;
dt=params.dt;

%%
% qdd = full( params.F.qddABA(q,qd,ui_un) );
% qdd_org = qdd;
% 
% M = full(params.F.all_first(q,qd,qdd));
% M_org = M;
% 
% q_val   = M(:,1:sz); 
% qd_val  = M(:,sz+1:sz*2);
% tau_val = M(:,sz*2+1:end);
% 
% qd_x=[zeros([length(q) length(q)]) eye(length(q))];
% AJ_2 = [q_val';qd_val']; 
% 
% fx = eye(x_size,x_size) + [qd_x;AJ_2']*dt;
% fu = [zeros(length(ui_un)+1,length(ui_un)); tau_val(:,1:length(ui_un))*dt];
% fx_org = fx;
% fu_org = fu;
% 
% mu_val = full( params.F.InvM_Mult(q,lambda2') ); 
% mu_org = mu_val;
% 
% [qq1, qd_qd1, q_qd1, q_tau, H_org, J_org, o_org] = params.F.mod_second_all(q,qd,qdd,mu_val,q_val,qd_val,tau_val);
% t3 = full([qq1 q_qd1;q_qd1' qd_qd1]*dt); 
% 
% t3_org = t3;
% 
% q_tau = full(q_tau(:,1:end-1));
% q_tau_org = q_tau;

%%
if params.modRNEA
    %use modRNEA
    [M, t3, q_tau] =  params.F.mod_all(q,qd,ui , lambda2);
else
    [M, t3, q_tau] =  params.F.All_second(q,qd,ui , lambda2);
end

M = full(M);
t3 = full(t3)*dt;

q_val   = M(:,1:sz); 
qd_val  = M(:,sz+1:sz*2);
tau_val = M(:,sz*2+1:end);

qd_x = [zeros([length(q) length(q)]) eye(length(q))];
AJ_2 = [q_val';qd_val']; 
fx   = eye(x_size,x_size) + [qd_x;AJ_2']*dt;
fu   = [zeros(length(ui_un)+1,length(ui_un)); tau_val(:,1:length(ui_un))*dt];

q_tau = full(q_tau(:,1:end-1));

lx = full(params.L.Lx(xi)); 
lu = full(params.L.Lu(ui(1:end-1)));

Hx = (lx + lambda'*fx)'; 
Hu = (lu + lambda'*fu)';

lxx = full(params.L.Lxx(xi));
Hxx = lxx + t3; 

Hux =  [q_tau*dt; zeros(size(q_tau))]; 
Huu = full(params.L.Luu(xi));

Qxx = Hxx  + fx'*Vxxi*fx; Qxx = .5 *(Qxx + Qxx'); 
Qux = Hux' + fu'*Vxxi*fx;
Quu = Huu  + fu'*Vxxi*fu; Quu = 0.5*(Quu+Quu');  

end