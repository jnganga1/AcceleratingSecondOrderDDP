function [Hx, Hu,Qxx, Quu, Qux] = CasadiQinfoUpdated_OldMethod(xi,ui, Vxi, Vxxi,params)

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

fx = eye(x_size,x_size) + [qd_x;AJ']*dt;

fu = [zeros(length(ui)+1,length(ui)); BJ(:,1:length(ui))]; 

lx = full(params.L.Lx(xi)); 
lu = full(params.L.Lu(ui)); 

[qq,q_qd,q_tau,qd_qd]=params.F.all_second(q,qd,ui); 

M11 = reshape(full(qq),sz); 
qqB = Tens3byVec(M11,lambda2','pre');

M12 = reshape(full(q_qd),sz);
qqdB = Tens3byVec(M12,lambda2','pre');

M13 = reshape(full(qd_qd),sz);
qdqdB = Tens3byVec(M13,lambda2','pre');

M14 = permute(M12,[1 3 2]); % This is a transpose
qqdB2 = Tens3byVec(M14,lambda2','pre'); 

t3 =  [qqB qqdB ;qqdB2 qdqdB]*dt;

lxx = full(params.L.Lxx(xi));
Hxx = lxx +t3; 

M15 = reshape(full(q_tau),sztau); %qtau= Tens3byVec(M15,lambda2','pre');
hux = permute(M15,[2 3 1]); hux = cat(3,zeros(size(hux)),hux);
a = permute(hux,[3 2 1]) .* dt; fuxP = cat(3,a,zeros(size(a)));

Hx = (lx + lambda'*fx)'; %Hx = lx + lambda'*fx 
Hu = (lu + lambda'*fu)'; %Hu = lu + lambda'*fu

Hux = full(params.L.Lux(xi,ui)) + Tens3byVec(fuxP,lambda','pre');
Huu = full(params.L.Luu(xi));

Qxx    = Hxx + fx'*Vxxi*fx; Qxx = .5 *(Qxx + Qxx'); 
Qux    = Hux + fu'*Vxxi*fx;
Quu    = Huu + fu'*Vxxi*fu; Quu = 0.5*(Quu+Quu');  

end