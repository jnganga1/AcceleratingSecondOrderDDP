function [Hx, Hu,HamFxx,Qxx, Quu, Qux,Fxx] = QinfoCasadi(xi,ui, Vxi, Vxxi,params)

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

qd_x=[zeros([length(q) length(q)]) eye(length(q))];
A1 =  params.F.q(q,qd,ui);  
A2 =  params.F.qd(q,qd);  
A = [full(A1)';full(A2)']; 
fx = eye(x_size,x_size) + [qd_x;A']*params.dt;
B = full(params.F.tau(q))*params.dt;
fu = [zeros(length(ui)+1,length(ui)); B(:,1:length(ui))]; 

lx = full(params.L.Lx(xi)); 
lu = full(params.L.Lu(ui)); 

Hx = (lx + lambda'*fx)'; %Hx = lx + lambda'*fx 
Hu = (lu + lambda'*fu)'; %Hu = lu + lambda'*fu

[qq,qqd,qtau,qdqd]=params.F.all_second(q,qd,ui); 
M11 = permute(reshape(full(qq),sz),[2 3 1]);%; M11 = permute(M11,[2 3 1]);
M12 = permute(reshape(full(qqd),sz),[2 3 1]);% M12 = permute(M12,[2 3 1]);
M13 = permute(reshape(full(qdqd),sz),[2 3 1]);% M13 = permute(M13,[2 3 1]);
M14 = permute(M12,[1 3 2]); % This is a transpose
M14 = permute(M14,[2 3 1]);

qqB = Tens3byVec(M11,lambda2','pre');
qqdB = Tens3byVec(M12,lambda2','pre');
qdqdB = Tens3byVec(M13,lambda2','pre');
qqdB2 = Tens3byVec(M14,lambda2','pre'); 

% hxx = [M11 M12;M14 M13];
% Fxx = cat(3,zeros(size(hxx)),hxx)*dt; Fxx = permute(Fxx,[3 1 2]);
% t2 = Tens3byVec(Fxx,lambda','pre'); t2 = 0.5*(t2+t2'); 

t2 = [qqB qqdB ;qqdB2 qdqdB]; t2 =0.5*(t2+ t2'); 

Fxx =1;
lxx = full(params.L.Lxx(xi));
Hxx = lxx +t2; 

1==1;
M15 = reshape(full(qtau),sztau); %qtau= Tens3byVec(M15,lambda2','pre');
hux = [zeros(size(M15)); M15]; hux= permute(hux,[2 3 1]);
fux = [zeros(size(hux));hux]*dt; 

Hux = full(params.L.Lux(xi,ui)) + Tens3byVec(fux,lambda','pre');
Huu = full(params.L.Luu(xi));

Qxx = Hxx + fx'*Vxxi*fx; Qxx = .5 *(Qxx + Qxx'); 
Qux = Hux + fu'*Vxxi*fx;
Quu = Huu + fu'*Vxxi*fu; Quu = 0.5*(Quu+Quu');  

HamFxx = Hxx - lxx;
1==1;

end