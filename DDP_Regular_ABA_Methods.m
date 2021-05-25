function Out = DDP_Regular_ABA_Methods(Nb,ver,iLQR,N,rbtNmber,x0_diff) %'end' is 351 

% ver=1;
% Nb=4;
% profile on -history 



%

addpath([pwd '/algorithm']);
addpath([pwd '/support']);
addpath([pwd '/spatial_v2']); %using a local copy of spatialV2
addpath(genpath([pwd '/spatial_v2_casadi']));
addpath(genpath([pwd])); %adds everything. Supersedes above adds

import casadi.*
warning('off','all')

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

sz = (Nb) * [1 1 1]; %Size of reshape of second derivs 
params.sz = sz;     
fprintf('Nb = %d \n',Nb);
 
% rbtNmber= 1;
% Generate the model
branching_factor = 1; %1=chain, 2=binary tree
skewing = pi/exp(1);  %0->2D mechanism, anything else -> 3D mechanism
% model = autoTree( Nb, branching_factor, skewing);
n = Nb;
model= robotModel(1/n,1/n,Nb,rbtNmber);

second_order =1;
if iLQR
    second_order = 0;
end
fprintf('\tComputing Direct Functions\n');
ABA_Deriv_Funcs = GetDerivativeFunctionsDirect(model,second_order);
% fprintf('\tComputing Alt Functions\n');
% ABA_Deriv_via_RNEA_Funcs    = GetDerivativeFunctionsAlt(model,second_order);
% fprintf('\tComputing Old Functions\n');
% funcsOld = GetDerivativeFunctionsOld(model,second_order);

%fnVar:: 1 = ABA_Deriv_Funcs
%fnVar:: 2 = ABA_Deriv_via_RNEA_Funcs
%fnVar:: 3 = funcsOld
%fnVar:: 4 = mod_RNEA

%Problem Setup 
% N=1200;%Horizon
% dt = 0.03; %state Initialization
% dt = 0.025*10^-1; %Cartpole

% N = 150;
dt= 0.0025;
% dt =0.03;

params.N = N;
params.dt = dt; %Incase I use it in other scripts 

x0 = zeros(size(x_MX));
% x0(1) = -pi/2 + 0.5*pi/2 ;% First link
x0(rbtNmber) = -pi/2;
x0 = x0 + x0_diff;
% % x0(rbtNmber) = pi/2 + x0_diff;
% x0(rbtNmber+1) =pi/2 + x0_diff; 
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


%Lqr Mats
x_size = length(x_desired);
qd_x=[zeros([length(q_MX) length(q_MX)]) eye(length(q_MX))];
A1 =  ABA_Deriv_Funcs.q(q_desired,qd_desired,u_desired);  
A2 =  ABA_Deriv_Funcs.qd(q_desired,qd_desired);  
A = [full(A1)';full(A2)']; 
A = eye(x_size,x_size) + [qd_x;A']*dt;
B = full(ABA_Deriv_Funcs.tau(q_desired))*dt;
B = [zeros(length(u)+1,length(u)); B(:,1:length(u))]; 

contrMat = ctrb(A,B); 
if ~(rank(A) == rank(contrMat))
    warning('Pair not controllable');   
    1==1;
end

[P,~,G] = dare(A,B,Q,R);
% P(1,1)=P(1,1)*70;
LF   = 1/2*(x_MX-x_desired)'*P*(x_MX-x_desired);  


Lfunc.L = Function('L', {x_MX,u}, {l} ); 
Lfunc.Lx = Function('Lx', {x_MX}, {jacobian(l,x_MX)});
Lfunc.Lu = Function('Lu', {u}, {jacobian(l,u)}); 
Lfunc.Lxx = Function('Lxx', {x_MX}, { hessian(l,x_MX)});
Lfunc.Lux = Function('Lux', {x_MX,u}, { jacobian(jacobian(l,u),x_MX) });
Lfunc.Luu =Function('Luu', {x_MX}, { hessian(l,u)});
Lfunc.LF = Function('LF', {x_MX},{LF});
Lfunc.LFderivs = Function('LFderivs', {x_MX},{jacobian(LF,x_MX),hessian(LF,x_MX)});

params.ver = ver; 
params.F = ABA_Deriv_Funcs; 
params.L = Lfunc;
1==1;
%%
%}
%%
% Parameters
params.iLQR = iLQR;
params.N = N;

params.beta = .5;
params.gamma = .01;
params.Debug = 0;

%PlaceHolder
callback_params.plot_V_prediction = 0;
callback_params.pause = 0;
callback_params.UserAdvanceIteration = 0;


% Initialize
Vbar = 0;
ubar = zeros(params.u_size, params.N); %1 *150
xbar = zeros(params.x_size, params.N); % 3*150
du = zeros(params.u_size, params.N); %1 *150
K = zeros(params.u_size, params.x_size, params.N);

xbar(:,1) = x0; %Store states
xi = x0;
% xi = zeros(size(x0));xi(2)=pi/2*.70;

% Make initial Guess
for i = 1:(params.N-1)
    
%     ui = -G*(xi-x_desired);% + .05*rand(Nb-1,1);
%     ui = 0;
%     + noiselvl(1:length(u),i)*.75;
     ui = -.0025*(xi(ceil(length(xi)/2)+1:end-1)); 


    ubar(:,i) = ui;
    xi = CasadiDyn(params.F.qddABA,xi,ui,dt);
    xbar(:,i+1) = xi;
end

%{
close all;
for i = 1:(Nb)
    figure; 
    hold on 
    h=plot(xbar(i,:),'DisplayName','state');h.LineWidth = 2;
    if  i==2
       b = plot(pi/2*ones(1,length(xbar)),'k--','DisplayName','pi/2 line');
       b.LineWidth =2;
    else    
        b=plot(zeros(1,length(xbar)),'k--','DisplayName','zero line');
        b.LineWidth =2;
    end
    legend
    hold off 
end

fprintf('\nPaused at line 239\n') 
%}
1==1;


%}

%%
[Vbar, xbar, ubar] = ForwardPass(xbar, ubar, du, K,1, params);

iter = 0;
regularization = 0; % for full second order DDP 
% profile clear
% profile on 
Vstore = [Vbar];
iLQR_store = [];
du_store = [];
DDPstart = tic;
idx2 = 1;
du_prev = zeros(1,N); 
norm_store = []; iLQR_du_store = [];DDP_du_store =[]; reg = 0;
change_store = []; change_holder = 0;
rho_store = [];FwdTrkcnter = [];BckTrkcnter=[];
Tracker.BckSuccess = 0; Tracker.FwdSuccess = 0;
Tracker.BckNonSuccess = 0; Tracker.FwdNonSuccess = 0;

% xf_dev = []; 
iterTimerTracker=[];
while 1 == 1
    iter= iter+1
    iterStart = tic;
    params.iter =iter;
    
    % Backward Pass 
    if params.iLQR
        [dV, Vx, Vxx, du, K] = BackwardPass_ABA(xbar, ubar, params);
        warning('Doing iLQR');
        BckTrkcnt=0;
    else % It is full second order DDP, so we must be careful
        success = 0;BckTrkcnt = 0;
        while success == 0
            if params.Debug
                fprintf('\t reg=%.3e\n',regularization);
            end
            bckTime = tic;
            [dV, Vx, Vxx, du, K, success] = BackwardPass_ABA(xbar, ubar, params,regularization);
            bckEndTime = toc(bckTime);
            if success == 0
                regularization = max(regularization*4, 1e-3);
                BckTrkcnt = BckTrkcnt + 1;
                Tracker.BckNonSuccess = Tracker.BckNonSuccess + bckEndTime;
            end
        end
        Tracker.BckSuccess = Tracker.BckSuccess + bckEndTime;
        regularization = regularization / 4;
        if regularization < 1e-6
            regularization = 0;
        end
    end
%     if iter ==5 
%         1==1;
% %         qdd_func =  ABA_Deriv_Funcs.qddABA;
%         save('DDP_gains.mat','du','K','xbar','ubar','params')
%         save('iLQR_gains.mat','du','K','xbar','ubar','params')
%     end
%     
    Vprev = Vbar;
    xprev =xbar; 
  
    %% Forward pass : Stepsize selection via backtracking line search
    [Vbar, xbar, ubar,rho,ThreeNumbers] = ForwardIteration(xbar, ubar, Vprev, du, K, dV , params);
    Tracker.FwdSuccess =  Tracker.FwdSuccess + ThreeNumbers.FwdSuc;  
    Tracker.FwdNonSuccess = Tracker.FwdNonSuccess + ThreeNumbers.FwdNonSuc;
    
%     xf_dev(end+1) = norm(xbar(1:Nb*2,end)- xprev(1:Nb*2,end));
    Vstore(end+1) = Vbar;
    Change = Vprev - Vbar 
    rho_store(end+1)=rho;
    FwdTrkcnter(end+1) = ThreeNumbers.FwdFailTrkCnt;
    BckTrkcnter(end+1) = BckTrkcnt;
%     change_holder = change_holder + Change; 
%     if ~mod(iter,5)
%         if change_holder/5 > 15 || change_holder/5 < 1
%             params.iLQR = 0;
%         else
%             params.iLQR = 1;
%         end 
%         change_store(end+1,:) = change_holder/5;
% 
%         change_holder = 0;
%     end
    

      iterTimerTracker(end+1)=toc(iterStart);
%     iter = iter+1;
    if Change < 1e-9
        fprintf('CONVERGED: Change: %f',Change)
        break
    end
end
Out.Vstore = Vstore;
Out.Time = toc(DDPstart);
Out.iterTimerTracker = iterTimerTracker; 
Out.Iters = iter;
1==1;
% Out.xf_dev = xf_dev;
% fprintf('Total time for DDP w/ TrustRegion was: %f\n seconds',outTime)



%% Plots 
%Undo The next %{ for plots
%{
time = 0:dt:dt*(params.N-1); eTime = ceil(dt*(params.N-1)); 
close all 
for i=1:length(q_MX)
    figure; 
    hold on 
    DP = plot(1:params.N,xbar(i,:),'b');
    if i==2
        ZP=plot(pi/2*ones (1,length(xbar)),'k--');
    else
        ZP=plot(pi/2*ones(1,length(xbar)),'k--');
    end
    DP.LineWidth =2;ZP.LineWidth =2;    
    hold off    
    xlabel('Time'); ylabel(['q',num2str(i)]);
    title(['State, q',num2str(i),', over time']);
    legend('State','Zero Line','interpreter','latex');
end
for i=1:length(qd_MX)
    figure; 
    hold on 
    k = i + length(q_MX);
    DP = plot(1:params.N,xbar(k,:),'b');
    ZP=plot(zeros(1,length(xbar)),'k--');
    DP.LineWidth =2;ZP.LineWidth =2;    
    hold off    
    xlabel('Time'); ylabel(['qd',num2str(i)]);
    title(['State, qd',num2str(i),', over time']);
    legend('State','Zero Line','interpreter','latex');
end

for i = 1:length(u)
    figure; 
    hold on
    ZP = plot(zeros(1,length(ubar(i,:))),'k--');
    DP = plot(1:params.N,ubar(i,:),'b');
    DP.LineWidth = 2; ZP.LineWidth = 2;
    hold off
    xlabel('Time'); ylabel(['U',num2str(i)]); title(['Control, U',num2str(i),', over time'])
    legend('Zero Line','Theta over Time');
end 
figure;  
hold on
DP = plot(Vstore,'b');
hold off
xlabel('Time'); ylabel('U');title('Value Over Time'); 

figure;
% sp(1) = subplot(2,2,1);
hold on;
set(gca, 'YScale','log');
CD = semilogy(Vstore-Vstore(end),'b-o');
CD.LineWidth = 0.5;
grid on
hold off 
xlabel('Iterations');
ylabel('Suboptimality');
title('Convergence Decay');

figure; 
plot(rho_store(1:end-1))
xlabel('Iters') 
ylabel('rho ratio')
text(iter-1, rho_store(end-1),'No Last pt \rightarrow ','HorizontalAlignment','right')
title('Rho ratio over time')

figure;
a=subplot(2,1,1);
c=plot(FwdTrkcnter);
xlabel('Iters') 
ylabel('Value')
title('''Fwd Line search failures''') 
c.LineWidth = 2;

b=subplot(2,1,2);
d=plot(BckTrkcnter); 
xlabel('Iters') 
title('''Rejected Bckwd Pass''')
ylabel('Value')
d.LineWidth =2;


vals = [Tracker.BckSuccess Tracker.BckNonSuccess Tracker.FwdSuccess Tracker.FwdNonSuccess];
vals = [vals; nan(1,4)];
vals = vals./sum(vals,2);
% vals(1,:) = vals(1,:)/sum(vals(1,:)); vals(2,:) = vals(2,:)/sum(vals(2,:));
figure; 
xvals = {'Regular DDP'; ''};
h=bar(vals,'stacked'); 
set(h, {'DisplayName'}, {'BckSuccess','BckNonSuccess','FwdSuccess','FwdNonSuccess'}')
set(gca, 'xticklabel',xvals)
legend('location','best')
title('Relative Time Spent')

%}
end


%{
fprintf('Final V: %.3f',Vbar) 
fprintf('Iters: %i',iter)


p = ['RegularDDP_data.txt'];
fileID = fopen(p,'a');
mm = ['Link',num2str(Number_Links)];
here = [mm,'::, Total Time: %f,',...
    ' Total Iters: %i, Time per Iter: %f, Final V Value: %f\n'];
hereVals = [outTime,iter,outTime/iter,Vbar];
fprintf(fileID,here,hereVals);
    
fclose(fileID);

FolderName= [searchPath,'/figures/'];
mkdir(FolderName);
save([FolderName,'rhoratio.mat'], 'rho_store')
save([FolderName,'RejectedIter.mat'],'BckTrkcnter','FwdTrkcnter')
save([FolderName,'TrackerInfo.mat'],'Tracker')
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName = ['figure',num2str(iFig)];
  savefig(FigHandle, [FolderName, FigName, '.fig']);
end

% end
%}

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