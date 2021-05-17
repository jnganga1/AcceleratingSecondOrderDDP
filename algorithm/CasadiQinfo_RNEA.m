function [Hx, Hu,HamFxx,Qxx, Quu, Qux,Fxx] = CasadiQinfoRNEA(xi,ui, Vxi, Vxxi,params)

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

%you're here - you're tracking if all_first is first than individual ones

qdd =params.F.qddABA(q,qd,ui_un);
% q_val= params.F.q(q,qd,qdd);
% qd_val =params.F.qd(q,qd); 
% tau_val =params.F.tau(q);



M =full(params.F.all_first(q,qd,qdd));
q_val = M(:,1:sz); 
qd_val = M(:,sz+1:sz*2);
tau_val =M(:,sz*2+1:end);
% AJ = full(M(:,1:sz*2)'); 
% BJ = full(M(:,sz*2+1:end)*dt);

qd_x=[zeros([length(q) length(q)]) eye(length(q))];
% fx = eye(x_size,x_size) + [qd_x;AJ']*dt;
AJ_2 = [q_val';qd_val']; 
fx = full(eye(x_size,x_size) + [qd_x;AJ_2']*dt);


% fu = [zeros(length(ui_un)+1,length(ui_un)); BJ(:,1:length(ui_un))]; 
fu = full([zeros(length(ui_un)+1,length(ui_un)); tau_val(:,1:length(ui_un))*dt]);


lx = full(params.L.Lx(xi)); 
lu = full(params.L.Lu(ui(1:end-1))); 

Hx = (lx + lambda'*fx)'; %Hx = lx + lambda'*fx 
Hu = (lu + lambda'*fu)'; %Hu = lu + lambda'*fu

% tic
% qdd = params.F.qddABA(q,qd,ui_un);
% q_val= params.F.q(q,qd,qdd);
% qd_val =params.F.qd(q,qd); 
% tau_val =params.F.tau(q);
% toc

mu_val = params.F.InvM_Mult(q,lambda2'); 
%%

if params.modRNEA == 0 % so it defaults to modRNEA
%Mod RNEA using ID relation
[qq, qd_qd, q_qd, d] =params.F.all_second(q,qd,qdd,mu_val,q_val,qd_val,tau_val);
% [a b c d2] = params.G.all_second(q,qd,ui_un,lambda2');
q_tau = d(:,1:end-1); q_qd = q_qd';
t3 =  [qq q_qd';q_qd qd_qd]*dt; t3 = full(0.5*(t3+t3')); 
else
%mod RNEA using casadi

% [H_mod, J_mod]=params.F.mod_second(q,qd,qdd,mu_val,q_val,qd_val,tau_val);
[qq1, qd_qd1, q_qd1, q_tau] = params.F.mod_second_all(q,qd,qdd,mu_val,q_val,qd_val,tau_val);

% A = - 1/2*H_mod(1:Nb,1:Nb) - J_mod(1:Nb,:);
% qq = (A+A');
% qd_qd = (-H_mod(Nb+1:end, Nb+1:end));
% q_qd = (-H_mod(1:Nb, Nb+1:end) - J_mod(Nb+1:2*Nb,:)')';
% q_tau = (-J_mod(2*Nb+1:end-1,:)');
% t3 =  [qq q_qd';q_qd qd_qd]*dt; t3 = full(0.5*(t3+t3')); 

t3 = [qq1 q_qd1;q_qd1' qd_qd1]*dt; t3 = full(0.5*(t3+t3')); 
q_tau = full(q_tau(:,1:end-1));
1==1;
end


1==1;

% 1==1;
% [qqC,qd_qdC,q_qdC,q_tauC]=params.G.all_second(q,qd,ui_un,lambda2'); 
% M11 = reshape(full(qq),sz); 
% qqB = Tens3byVec(M11,lambda2','pre');

% M12 = reshape(full(q_qd),sz);
% qqdB = Tens3byVec(M12,lambda2','pre');

% M13 = reshape(full(qd_qd),sz);
% qdqdB = Tens3byVec(M13,lambda2','pre');

% M14 = permute(M12,[1 3 2]); % This is a transpose
% qqdB2 = Tens3byVec(M14,lambda2','pre'); 

% M11 = permute(M11,[2 3 1]);
% M12 = permute(M12,[2 3 1]);
% M13 = permute(M13,[2 3 1]);
% M14 = permute(M14,[2 3 1]);

% hxx = [M11 M12;M14 M13];
% Fxx = cat(3,zeros(size(hxx)),hxx)*dt; Fxx = permute(Fxx,[3 1 2]);
% t2 = Tens3byVec(Fxx,lambda','pre');% t2 = 0.5*(t2+t2'); 

% t3 =  [qq q_qd';q_qd qd_qd]*dt; t3 = full(0.5*(t3+t3')); 

% M11 = permute(reshape(full(qq),sz),[2 3 1]); %M11 = permute(M11,[2 3 1]);
% M12 = permute(reshape(full(q_qd),sz),[2 3 1]);% M12 = permute(M12,[2 3 1]);
% M13 = permute(reshape(full(qd_qd),sz),[2 3 1]);% M13 = permute(M13,[2 3 1]);
% M14 = permute(M12,[1 3 2]); % This is a transpose
% M14 = permute(M14,[2 3 1]);
% 
% qq = Tens3byVec(M11,lambda2','pre');
% q_qd = Tens3byVec(M12,lambda2','pre');
% qd_qd = Tens3byVec(M13,lambda2','pre');
% q_qdO = Tens3byVec(M14,lambda2','pre'); 
% 
% hxx = [M11 M12;M14 M13];
% Fxx = cat(3,zeros(size(hxx)),hxx)*dt; Fxx = permute(Fxx,[3 1 2]);
% t2 = Tens3byVec(Fxx,lambda','pre'); t2 = 0.5*(t2+t2'); 
% 
Fxx = 1; %Does not matter for now. Come back
% t2 = [qq q_qd;q_qdO qd_qd]; t2 =full(0.5*(t2+ t2')); 

lxx = full(params.L.Lxx(xi));

Hxx = lxx +t3; 

% 1==1;
% M15 = reshape(full(q_tau),sztau); %qtau= Tens3byVec(M15,lambda2','pre');
% hux = [zeros(size(M15)); M15]; hux= permute(hux,[2 3 1]);
% fux = [zeros(size(hux));hux]*dt; 
% fux = full(q_tau);

% Hux = full(params.L.Lux(xi,ui)) + Tens3byVec(fux,lambda','pre');
Hux =  [full(q_tau*dt); zeros(size(q_tau))]; 
Huu = full(params.L.Luu(xi));

Qxx = Hxx + fx'*Vxxi*fx; Qxx = .5 *(Qxx + Qxx'); 
Qux = Hux' + fu'*Vxxi*fx;
Quu = Huu + fu'*Vxxi*fu; Quu = 0.5*(Quu+Quu');  

HamFxx = Hxx - lxx;
1==1;

end