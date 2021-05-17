function [Hx, Hu,HamFxx,Qxx, Quu, Qux,Fxx] = CasadiQinfoUpdated(xi,ui, Vxi, Vxxi,params)

% xi = zeros(size(xi)); 
% ui = zeros(size(ui)); 
% Vxi = zeros(size(Vxi)); Vxi(1) = 0; Vxi(2) = 0;
% Vxxi = zeros(size(Vxxi));


% Some definitions    
x_size = length(xi);
q =xi(1:length(xi)/2); 
qd = xi(length(xi)/2+1:end);

lambda = Vxi;
lambda1 = Vxi(1:length(xi)/2); %for simplification sake
lambda2 = Vxi(length(xi)/2+1:end);
sz = params.sz;
sztau = sz; sztau(end)= sz(end)-1;
dt=params.dt;

M = params.F.all_first(q,qd,ui);

AJ =full(M(:,1:sz*2)'); %
BJ = full(M(:,sz*2+1:end)*dt);



qd_x=[zeros([length(q) length(q)]) eye(length(q))];
% A1 =  params.F.q(q,qd,ui_un);  %
% A2 =  params.F.qd(q,qd);  %
% A = [full(A1)';full(A2)']; %
fx = eye(x_size,x_size) + [qd_x;AJ']*dt;
% B = full(params.F.tau(q))*dt;%
fu = [zeros(length(ui)+1,length(ui)); BJ(:,1:length(ui))]; 

lx = full(params.L.Lx(xi)); 
lu = full(params.L.Lu(ui)); 

Hx = (lx + lambda'*fx)'; %Hx = lx + lambda'*fx 
Hu = (lu + lambda'*fu)'; %Hu = lu + lambda'*fu

[qq,q_qd,q_tau,qd_qd]=params.F.all_second(q,qd,ui); 
M11 = reshape(full(qq),sz); 
qqB = Tens3byVec(M11,lambda2','pre');

M12 = reshape(full(q_qd),sz);
qqdB = Tens3byVec(M12,lambda2','pre');

M13 = reshape(full(qd_qd),sz);
qdqdB = Tens3byVec(M13,lambda2','pre');

M14 = permute(M12,[1 3 2]); % This is a transpose
qqdB2 = Tens3byVec(M14,lambda2','pre'); 

% M11 = permute(M11,[2 3 1]);
% M12 = permute(M12,[2 3 1]);
% M13 = permute(M13,[2 3 1]);
% M14 = permute(M14,[2 3 1]);

% hxx = [M11 M12;M14 M13];
% Fxx = cat(3,zeros(size(hxx)),hxx)*dt; Fxx = permute(Fxx,[3 1 2]);
% t2 = Tens3byVec(Fxx,lambda','pre'); t2 = 0.5*(t2+t2'); 
Fxx=1; %Not needed
t3 =  [qqB qqdB ;qqdB2 qdqdB]*dt; t3 = 0.5*(t3+t3'); 

lxx = full(params.L.Lxx(xi));
Hxx = lxx +t3; 

1==1;
M15 = reshape(full(q_tau),sztau); %qtau= Tens3byVec(M15,lambda2','pre');
hux = permute(M15,[2 3 1]); hux = cat(3,zeros(size(hux)),hux);
% hux = [M15; zeros(size(M15))]; hux= permute(hux,[2 3 1]);
% fux = [zeros(size(hux));hux]*dt; 
% fux = full(q_tau);

a = permute(hux,[3 2 1]) .* dt; fuxP = cat(3,a,zeros(size(a)));


Hux = full(params.L.Lux(xi,ui)) + Tens3byVec(fuxP,lambda','pre');
% Hux = [zeros(size(q_tau)); full(q_tau)];
Huu = full(params.L.Luu(xi));

Qxx = Hxx + fx'*Vxxi*fx; Qxx = .5 *(Qxx + Qxx'); 
Qux = Hux + fu'*Vxxi*fx;
Quu = Huu + fu'*Vxxi*fu; Quu = 0.5*(Quu+Quu');  

HamFxx = Hxx - lxx;
1==1;

end