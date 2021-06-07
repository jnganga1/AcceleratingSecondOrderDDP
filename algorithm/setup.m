%setup function - don't run alone 
%assumes variable Nb x0_diff exist in the workspace

addpath([pwd '/algorithm']);
addpath([pwd '/support']);
addpath([pwd '/spatial_v2']); %using a local copy of spatialV2
addpath(genpath([pwd '/spatial_v2_casadi']));

import casadi.*
warning('off','all')

% saver = 0;

%Model Setup
% Nb = 3;
%Symbolics
q_MX =  MX.sym('q_MX',[Nb 1]);
qd_MX = MX.sym('qd_MX',[Nb 1]);
x_MX = [q_MX ; qd_MX];   % State    

u = MX.sym('u',[Nb-1 1]); %Last link not controlled

u_size = Nb -1; x_size = length(x_MX);
params.u_size = u_size;
params.x_size =x_size;
% params.modRNEA =modRNEA;

sz = (Nb) * [1 1 1]; %Size of reshape of second derivs 
params.sz = sz;     
fprintf('Nb = %d \n',Nb);
 
% rbtNmber= 1; %1==acrobot, 2 ==  cartpole

% Generate the model
% branching_factor = 1; %1=chain, 2=binary tree
% skewing = pi/exp(1);  %0->2D mechanism, anything else -> 3D mechanism
% model = autoTree( Nb, branching_factor, skewing);

n = Nb;
model= robotModel(1/n,1/n,Nb,rbtNmber);


% N = 150;
dt= 0.0025;
params.dt = dt; %Incase I use it in other scripts 

x0 = zeros(size(x_MX));
% x0(1) = -pi/2 + 0.5*pi/2 ;% First link
% x0(rbtNmber) = pi/2 + 0.150;%
x0(rbtNmber) = -pi/2; 
x0 = x0 + x0_diff;
% x0(rbtNmber) =x0(rbtNmber)+0.01;
% x0 = x0';

q_desired= zeros(size(q_MX));
qd_desired=zeros(size(qd_MX));
u_desired = zeros(length(u),1);
q_desired(rbtNmber) = pi/2;

R = diag(500*ones(1,u_size)); 
Q_pos = diag(900*ones(1,size(q_MX,1)));
Q_vel = diag(400*ones(1,size(q_MX,1))); 
Q = [Q_pos zeros(size(Q_pos));zeros(size(Q_vel)) Q_vel];
Q(1,1) = 1200;
x_desired= [q_desired;qd_desired];


l_cont = (x_MX- x_desired)'*Q*(x_MX-x_desired) +  (u)'*R*(u);
l = l_cont*dt;
vars = [x_MX; u];


Lfunc.L = Function('L', {x_MX,u}, {l} ); 
Lfunc.Lx = Function('Lx', {x_MX}, {jacobian(l,x_MX)});
Lfunc.Lu = Function('Lu', {u}, {jacobian(l,u)}); 
Lfunc.Lxx = Function('Lxx', {x_MX}, { hessian(l,x_MX)});
Lfunc.Lux = Function('Lux', {x_MX,u}, { jacobian(jacobian(l,u),x_MX) });
Lfunc.Luu =Function('Luu', {x_MX}, { hessian(l,u)});






%%
function robot = robotModel(mass_pole,length_pole,N_links,whichRobot) 
%1 == acrobot 
%2 == Cartpole 


if whichRobot == 1 
    %N-linked Acrobot 
    robot.NB = N_links; 
    robot.parent =[0:N_links-1];

    robot.gravity = [0 -9.81]';

   for idx =1:N_links
        robot.jtype{idx} ='r'; 
        robot.I{idx} = mcI(mass_pole,[length_pole 0],0); 
    %        robot.I{idx} = mcI(mass_pole,[length_pole/2 0],1/12*mass_pole*length_pole^2); 
        if idx == 1  
           robot.Xtree{idx} = eye(3); %1
    %            robot.Xtree{idx} = plnr(deg2rad(-180),[0 0]); %2
        else
           robot.Xtree{idx} = plnr(0,[length_pole 0]);
        end      
   end
else
    %N-linked Cartpole
        
%         NB = N_links+1;
        robot.NB = N_links;
        robot.parent = [0:N_links-1];
        robot.gravity = [0 -9.81]';

%         NB = N_links+1;
%         robot.NB = NB;
%         robot.parent = [0:N_links];
%         robot.gravity = [0 -9.81]';


        robot.jtype{1} = 'px'; % The first joint is a prismatic in the +X0 direction
        robot.jtype{2} = 'r';

        robot.Xtree{1} = eye(3); % The transform from {0} to just before joint 1 is identity
        robot.Xtree{2} = eye(3);

        mass_cart = 1;
        robot.I{1} = mcI( mass_cart, [0 0], 0 );
        robot.I{2} = mcI( mass_pole, [length_pole 0], 0 );    

        for i=3:N_links
            robot.jtype{i} = 'r'; %Revolute in +z1; 
            robot.Xtree{i} = plnr(0,[length_pole 0]); 
    %           robot.Xtree{i+1} = xlt([length_pole 0 0]);
            robot.I{i}= mcI( mass_pole, [length_pole 0], 0 );

        end
end   
end