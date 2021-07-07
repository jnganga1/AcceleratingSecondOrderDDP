%%
% simulation of Ornstein–Uhlenbeck adapted from 
% https://math.stackexchange.com/questions/1287634/implementing-ornstein-uhlenbeck-in-matlab
% 
% clear; close all;
% 
% 
% N =250; 
% 
% th = 1;
% mu = 1.2;
% sig = 0.3;
% dt = 1e-2;
% %t = 0:dt:2;            % Time vector
% t = linspace(0,10,N); %linspace(0,2,N); 
% dt = t(2) - t(1);
% fwdU = zeros(1,length(t)); % Allocate output vector, set initial condition
% rng(1);                 % Set random seed
% for i = 1:length(t)-1
%     fwdU(i+1) = fwdU(i)+th*(mu-fwdU(i))*dt+sig*sqrt(dt)*randn;
% end


%%
clear all; 
rng('default');
num_sims =40;
% rndIdx = randi([1, length(fwdU)],[1 num_sims]); %Choose five indices of openloop ctrl


cd 'DDP_Regular_ABA_Methods - LBRModel'

addpath(genpath('../casadi-windows-matlabR2016a-v3.5.3'));
addpath(genpath('../casadi-linux-matlabR2014b-v3.5.5'));
addpath(genpath('../spatial_v2'));
addpath(genpath('../spatial_v2_casadi')); 
addpath(genpath([pwd])); %adds everything. Supersedes above adds


N = 500;
Nb = 7;



Output = zeros(4,num_sims);
Output_iLQR =Output;

for i=1:num_sims
%    rndIdx = randi(length(fwdU), 1); %random index 
%    cntrl = fwdU; %Feedforward Gain
   for idx=1:Nb*2
       cntrl(idx,:) = NoiseGen(N); 
   end
    Out_ddp=DDP_Regular_ABA_Methods_LBR(0,N,cntrl)
%    Out_ddp = DDP_Regular_OldMethod(Nb,0,N,1,0,cntrl);
   Output(:,i) = [Out_ddp.Vstore(end),sum(Out_ddp.Vstore),Out_ddp.Time,Out_ddp.Iters];
%    Out_ilqr = DDP_Regular_OldMethod(Nb,1,N,1,0,cntrl);
    Out_ilqr = DDP_Regular_ABA_Methods_LBR(1,N,cntrl)
   Output_iLQR(:,i) = [Out_ilqr.Vstore(end),sum(Out_ilqr.Vstore),Out_ilqr.Time,Out_ilqr.Iters];
end

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultAxesFontSize', 18);

save lbr_randomization.mat
care = 2;
figure; %subplot(2,2,1)
bar([Output(care,:);Output_iLQR(care,:)]');
legend('DDP','iLQR');
xlabel('Simulations'); ylabel('Total Cost'); 

%%
figure;
care =1;
% figure;
subplot(2,2,[3 4])
bar([Output(care,:);Output_iLQR(care,:)]');
legend('DDP','iLQR');
xlabel('Simulations'); ylabel('Cost at Last Iteration'); 

care =4;
% figure;
subplot(2,2,[1 2])
bar([Output(care,:);Output_iLQR(care,:)]');
legend('DDP','iLQR');
xlabel('Simulations'); ylabel('Number of Iterations'); 



%pdf 
figure;
care =4;
% figure;
subplot(2,2,[1 2])
bar([Output(care,:);Output_iLQR(care,:)]');
legend('DDP','iLQR');
xlabel('Simulations'); ylabel('Number of Iterations'); 


care =1;
% figure;
subplot(2,2,[3 4])
histogram(log10(Output(care,:)),'Normalization','pdf','DisplayName','DDP','BinWidth',0.2); hold on 
histogram(log10(Output_iLQR(care,:)),'Normalization','pdf','DisplayName','iLQR','BinWidth',0.2);
legend
ylabel('PDF'); xlabel('Log Cost at Last Iteration'); 


%% Save pdfs of current figure
%{
h = gcf; 
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'Sim_lbr','-dpdf','-r0')
%}

%%

% plot(NoiseGen(250))
function fwdU = NoiseGen(nmbr)
    rng('shuffle');
    th = 1;
%     mu = 1.2;
    mu = 5;
    sig = 1.2;%sig = 0.6;
    dt = 1e-2;
    %t = 0:dt:2;            % Time vector
    t = linspace(randn,10,nmbr); %linspace(0,2,N); 
    dt = t(2) - t(1);
    fwdU = zeros(1,length(t)); % Allocate output vector, set initial condition
    fwdU(1) = randn;
    rng(1);                 % Set random seed
    for i = 1:length(t)-1
        fwdU(i+1) = fwdU(i)+th*(mu-fwdU(i))*dt+sig*sqrt(dt)*randn;
    end
end