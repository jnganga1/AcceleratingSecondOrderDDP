% Script to compare effects of varying Horizon length
% Ensure that initial control is as desired 

%rbtNmbr: 1==Acrobot; 2==Cartpole (they use same L func) 
%DDP_Regular_ABA_Methods(iNb,ver,iLQR,150,rbtNmber,x0_diff);
    %ver: 1== _____  
%DDP_RegularRNEA(Nb,iLQR,N,carp,rbtNmber,x0_diff); 
   %carp: 1==carpentier; 0==Modified 
%DDP_Regular_ABA_Methods(Nb,ver,iLQR,N,rbtNmber,x0_diff)
% DDP_Regular_OldMethod(Nb,iLQR,N,rbtNmber,x0_diff)

%Ensure x0_diff does what you want :: currently it's doing this
%x0_diff changes the following
% x0 = zeros(2*n,1);n 
% x0(rbtNmber) = pi/2 + x0_diff

%Select robot here - ensures all code runs the same thing
rbtNmber =1;
addpath(genpath([pwd])); %adds everything. Supersedes above adds


%% Make sure scripts run want you want!!! check 
%% DDP_Regular_ABA_Methods(); DDP_RegularRNEA(); DDP_Regular_OldMethod();
%% Comparing the effects of varying Horizon Length

x0_diff =0.15;

VstoreA={};
VStoreB={};
VStoreC={};
VStoreD={};
VStoreE = {};

N = [100 150 200 250 300 400];
nbd = [4 8 10]
for  iNb = nbd
    %Didn't bother to avg repititions
    % 
    Timer_ABA = []; Timer_ABA_iLQR=[];
    Iters_ABA = []; Iters_ABA_iLQR=[]; 

    Timer_RNEA =[]; Timer_RNEA_iLQR=[];
    Iters_RNEA =[]; Iters_RNEA_iLQR =[];
    TimerRNEA_Carp=[]; ItersRNEA_Carp=[];

    Old_Time = []; Old_Iters =[];
    Old_Time_iLQR=[];Old_Iters_iLQR =[];


    Vstore_ABA_iLQR = [];
    Vstore_ABA = []; 
    Vstore_old =[];
    idx = 1; 
    for iN = N

        Out= DDP_Regular_ABA_Methods(iNb,2,0,iN,rbtNmber,x0_diff);
        OutiLQR= DDP_Regular_ABA_Methods(iNb,2,1,iN,rbtNmber,x0_diff);
        Vstore_ABA_iLQR(end+1) = OutiLQR.Vstore(end);
        Vstore_ABA(end+1) = Out.Vstore(end);
        Timer_ABA(end+1) = Out.Time; Iters_ABA(end+1)=Out.Iters;
        Timer_ABA_iLQR(end+1)=OutiLQR.Time; Iters_ABA_iLQR(end+1)=OutiLQR.Iters;

        Out= DDP_Regular_OldMethod(iNb,0,iN,rbtNmber,x0_diff); %DDP
        Old_Time(end+1)=Out.Time; Old_Iters(end+1)=Out.Iters; 
        Vstore_old(end+1) = Out.Vstore(end);

    end
    fname= ['n',num2str(iNb),'data']
    save(fname)

% 
% figure;
% plot3(N,Old_Time,Vstore_old,'DisplayName','DDP via Tensor Contraction'); hold on 
% plot3(N,Timer_ABA_iLQR,Vstore_ABA, 'DisplayName','iLQR via ABA')
% plot3(N,Timer_ABA,Vstore_ABA_iLQR, 'DisplayName', 'DDP via ABA')
% ylabel('Time'); xlabel('N'); 
% set(gca,'YScale','log')
% legend
% grid on; grid minor; title('$n=6$')
% 

%% Saves fig as pdf without improperly cropping them 
%Remove next { to activate code  
%
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'TimeVaryn10N','-dpdf','-r0')
%}
%%
%%
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultLineLineWidth', 5);


figure;
plot(N,Old_Time,'DisplayName','DDP via Tensor Contraction'); hold on 
plot(N,Timer_ABA_iLQR, 'DisplayName','iLQR via ABA')
plot(N,Timer_ABA, 'DisplayName', 'DDP via ABA')
ylabel('Time'); xlabel('N'); 
set(gca,'YScale','log')
legend
grid on; grid minor; title('$n=6$')

figure;
plot(N,Old_Iters,'DisplayName','DDP via Tensor Contraction'); hold on
plot(N,Iters_ABA_iLQR, 'DisplayName','iLQR via ABA')
plot(N,Iters_ABA, 'DisplayName', 'DDP via ABA')
ylabel('Iterations'); xlabel('N'); 
set(gca,'YScale','log')
legend
grid on; grid minor;
sTitle = ['n=',num2str(iNb)]
title(sTitle)
%%

end
%%
%{
figure;
hold on;
set(gca, 'YScale','log');

DE1 = semilogy(VStoreA.ABA-VStoreA.ABA(end),'b-o','DisplayName',...
    strcat('DDP via ABA, N=',string(VStoreA.N))); DE1.LineWidth = 2.5; 
DE2 = semilogy(VStoreB.ABA-VStoreB.ABA(end),'k-*','DisplayName',...
    strcat('DDP via ABA, N=',string(VStoreB.N))); DE2.LineWidth = 2.5;
DE3 = semilogy(VStoreC.ABA-VStoreC.ABA(end),'r-s','DisplayName',...
    strcat('DDP via ABA, N=',string(VStoreC.N))); DE3.LineWidth = 2.5;
DE4 = semilogy(VStoreD.ABA-VStoreD.ABA(end),'m-d','DisplayName',...
    strcat('DDP via ABA, N=',string(VStoreD.N))); DE4.LineWidth = 2.5;


DF1 = semilogy(VStoreA.ABA_iLQR-VStoreA.ABA_iLQR(end),'b-.o','DisplayName',...
    strcat('iLQR via ABA, N=',string(VStoreA.N))); DF1.LineWidth = 2.5;
DF2 = semilogy(VStoreB.ABA_iLQR-VStoreB.ABA_iLQR(end),'k-.*','DisplayName',...
    strcat('iLQR via ABA, N=',string(VStoreB.N))); DF2.LineWidth = 2.5;
DF3 = semilogy(VStoreC.ABA_iLQR-VStoreC.ABA_iLQR(end),'r-.s','DisplayName',...
    strcat('iLQR via ABA, N=',string(VStoreC.N))); DF3.LineWidth = 2.5;
DF4 = semilogy(VStoreD.ABA_iLQR-VStoreD.ABA_iLQR(end),'m-.d','DisplayName',...
    strcat('iLQR via ABA, N=',string(VStoreD.N))); DF4.LineWidth = 2.5;

grid on
hold off 
ld = legend; ld.FontSize = 20; legend
xlabel('Iterations');
ylabel('Suboptimality');
title(strcat('ABA: Suboptimality Decay vs Iterations for ',string(VStoreA.nbd),' bodies'));

figure;
hold on;
set(gca, 'YScale','log');
DG1 = semilogy(VStoreA.RNEA-VStoreA.RNEA(end),'b-o','DisplayName',...
    strcat('DDP via RNEA, N=',string(VStoreA.N))); DG1.LineWidth = 2.5; 
DG2 = semilogy(VStoreB.RNEA-VStoreB.RNEA(end),'k-*','DisplayName',...
    strcat('DDP via RNEA, N=',string(VStoreB.N))); DG2.LineWidth = 2.5;
DG3 = semilogy(VStoreC.RNEA-VStoreC.RNEA(end),'r-s','DisplayName',...
    strcat('DDP via RNEA, N=',string(VStoreC.N))); DG3.LineWidth = 2.5;
DG4 = semilogy(VStoreD.RNEA-VStoreD.RNEA(end),'m-d','DisplayName',...
    strcat('DDP via RNEA, N=',string(VStoreD.N))); DG4.LineWidth = 2.5;


DT1 = semilogy(VStoreA.RNEA_iLQR-VStoreA.RNEA_iLQR(end),'b-.o','DisplayName',...
    strcat('iLQR via RNEA, N=',string(VStoreA.N))); DT1.LineWidth = 2.5;
DT2 = semilogy(VStoreB.RNEA_iLQR-VStoreB.RNEA_iLQR(end),'k-.*','DisplayName',...
    strcat('iLQR via RNEA, N=',string(VStoreB.N))); DT2.LineWidth = 2.5;
DT3 = semilogy(VStoreC.RNEA_iLQR-VStoreC.RNEA_iLQR(end),'r-.s','DisplayName',...
    strcat('iLQR via RNEA, N=',string(VStoreC.N))); DT3.LineWidth = 2.5;
DT4 = semilogy(VStoreD.RNEA_iLQR-VStoreD.RNEA_iLQR(end),'m-.d','DisplayName',...
    strcat('iLQR via RNEA, N=',string(VStoreD.N))); DT4.LineWidth = 2.5;

grid on
hold off 
ld = legend; ld.FontSize = 20; legend
xlabel('Iterations');
ylabel('Suboptimality');
title(strcat('RNEA: Suboptimality Decay vs Iterations for ',string(VStoreA.nbd),' bodies'));


figure;
hold on;
set(gca, 'YScale','log');
DG1 = semilogy(VStoreA.Old_DDP-VStoreA.Old_DDP(end),'b-o','DisplayName',...
    strcat('DDP via Tensor Contraction, N=',string(VStoreA.N))); DG1.LineWidth = 2.5; 
DG2 = semilogy(VStoreB.Old_DDP-VStoreB.Old_DDP(end),'k-*','DisplayName',...
    strcat('DDP via Tensor Contraction, N=',string(VStoreB.N))); DG2.LineWidth = 2.5;
DG3 = semilogy(VStoreC.Old_DDP-VStoreC.Old_DDP(end),'r-s','DisplayName',...
    strcat('DDP via Tensor Contraction, N=',string(VStoreC.N))); DG3.LineWidth = 2.5;
DG4 = semilogy(VStoreD.Old_DDP-VStoreD.Old_DDP(end),'m-d','DisplayName',...
    strcat('DDP via Tensor Contraction, N=',string(VStoreD.N))); DG4.LineWidth = 2.5;


DT1 = semilogy(VStoreA.Old_iLQR-VStoreA.Old_iLQR(end),'b-.o','DisplayName',...
    strcat('iLQR via Old Method(ABA), N=',string(VStoreA.N))); DT1.LineWidth = 2.5;
DT2 = semilogy(VStoreB.Old_iLQR-VStoreB.Old_iLQR(end),'k-.*','DisplayName',...
    strcat('iLQR via Old Method(ABA), N=',string(VStoreB.N))); DT2.LineWidth = 2.5;
DT3 = semilogy(VStoreC.Old_iLQR-VStoreC.Old_iLQR(end),'r-.s','DisplayName',...
    strcat('iLQR via Old Method(ABA), N=',string(VStoreC.N))); DT3.LineWidth = 2.5;
DT4 = semilogy(VStoreD.Old_iLQR-VStoreD.Old_iLQR(end),'m-.d','DisplayName',...
    strcat('iLQR via Old Method(ABA), N=',string(VStoreD.N))); DT4.LineWidth = 2.5;

grid on
hold off 
ld = legend; ld.FontSize = 20; legend
xlabel('Iterations');
ylabel('Suboptimality');
title(strcat('Tensor Contraction: Suboptimality Decay vs Iterations for ',string(VStoreA.nbd),' bodies'));
%}