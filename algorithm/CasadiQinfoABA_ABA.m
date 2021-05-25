function [Hx, Hu,HamFxx,Qxx, Quu, Qux,Fxx] = CasadiQinfoABA_ABA(xi,ui, Vxi, Vxxi,params)

ver = params.ver;
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
ui_un = [ui;0];

if ver == 1
%     [qq,qd_qd,q_qd,q_tau]=params.F.all_second(q,qd,ui,lambda2'); 
       
    [M,t3,q_tau] = params.F.All_second(q,qd,ui,lambda2'); 
    
%     M =params.F.all_first(q,qd,ui);
    AJ = full(M(:,1:sz*2)'); 
    BJ = full(M(:,sz*2+1:end)*dt);
    qd_x=[zeros([length(q) length(q)]) eye(length(q))];
    fx = eye(x_size,x_size) + [qd_x;AJ']*dt;
    % B = full(params.F.tau(q))*dt; %
    fu = [zeros(length(ui)+1,length(ui)); BJ(:,1:length(ui))];

else
    [M,t3,q_tau] = params.F.All_second_V2(q,qd,ui,lambda2'); 

%         M =params.F.all_first(q,qd,ui);
    AJ = full(M(:,1:sz*2)'); 
    BJ = full(M(:,sz*2+1:end)*dt);
    qd_x=[zeros([length(q) length(q)]) eye(length(q))];
    fx = eye(x_size,x_size) + [qd_x;AJ']*dt;
    % B = full(params.F.tau(q))*dt; %
    fu = [zeros(length(ui)+1,length(ui)); BJ(:,1:length(ui))];
end


lx = full(params.L.Lx(xi)); 
lu = full(params.L.Lu(ui)); 

Hx = (lx + lambda'*fx)'; %Hx = lx + lambda'*fx 
Hu = (lu + lambda'*fu)'; %Hu = lu + lambda'*fu

% t3 =  [qq q_qd';q_qd qd_qd]*dt;
t3 = full(0.5*dt*(t3+t3')); 

Fxx = 1; %Does not matter for now. Come back
lxx = full(params.L.Lxx(xi));
Hxx = lxx +t3; 

Hux =  [full(q_tau*dt); zeros(size(q_tau))];
Huu = full(params.L.Luu(xi));

Qxx = Hxx + fx'*Vxxi*fx; Qxx = .5 *(Qxx + Qxx'); 
Qux = Hux' + fu'*Vxxi*fx;
Quu = Huu + fu'*Vxxi*fu; Quu = 0.5*(Quu+Quu');  

HamFxx = Hxx - lxx;
1==1;

end