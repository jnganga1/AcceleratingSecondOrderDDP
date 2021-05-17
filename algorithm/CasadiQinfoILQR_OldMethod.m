function [Hx, Hu,Qxx, Quu, Qux] = CasadiQinfoILQR(xi,ui, Vxi, Vxxi,params)

% Some definitions    
x_size = length(xi);
q =xi(1:length(xi)/2); 
qd = xi(length(xi)/2+1:end);
ui_un = ui;
ui = [ui;0];
sz = params.sz; 
dt = params.dt;

lambda = Vxi;
% lambda1 = Vxi(1:length(xi)/2); %for simplification sake
% lambda2 = Vxi(length(xi)/2+1:end);
% sz = params.sz;
% sztau = sz; sztau(end)= sz(end)-1;
% dt=params.dt;

%
M = params.F.all_first(q,qd,ui_un);

AJ =full(M(:,1:sz*2)'); %
BJ = full(M(:,sz*2+1:end)*dt);



qd_x=[zeros([length(q) length(q)]) eye(length(q))];
% A1 =  params.F.q(q,qd,ui_un);  %
% A2 =  params.F.qd(q,qd);  %
% A = [full(A1)';full(A2)']; %
fx = eye(x_size,x_size) + [qd_x;AJ']*dt;
% B = full(params.F.tau(q))*dt;%
fu = [zeros(length(ui_un)+1,length(ui_un)); BJ(:,1:length(ui_un))]; 


lx = full(params.L.Lx(xi)); 
lu = full(params.L.Lu(ui(1:end-1))); 

Hx = (lx + lambda'*fx)'; %Hx = lx + lambda'*fx 
Hu = (lu + lambda'*fu)'; %Hu = lu + lambda'*fu

lxx = full(params.L.Lxx(xi));
Hxx = lxx;

Hux = full(params.L.Lux(xi,ui_un)); %This is just zero
Huu = full(params.L.Luu(xi));

Qxx = Hxx + fx'*Vxxi*fx; Qxx = .5 *(Qxx + Qxx'); % good
Qux = Hux + fu'*Vxxi*fx; % good
Quu = Huu + fu'*Vxxi*fu; Quu = 0.5*(Quu+Quu');  

end