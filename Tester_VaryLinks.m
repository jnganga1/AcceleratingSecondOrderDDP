%% varying number of links
%looking at specifics
%This section -> varying number of links
%Compares Suboptimality vs Iterations 
%Compares EvaluationTime vs number of links

addpath(genpath([pwd])); %adds everything. Supersedes above adds

%Select robot here - ensures all code runs the same thing
rbtNmber =1;
x0_diff = 0.15;

N = 400;
% nbd = [2 3 4 5 7 9 12 16 22 26 30];
nbd = [2 4 6 10 15 19 20]; % 25 30];
% nbd=[20 25 30 35 40]%16 20];
% nbd = 15;

% nbd = 4
nmbRepts =2;


%ABA::  DDP_Regular_ABA_Methods(Nb,ver,iLQR,N,rbtNmber,x0_diff)
Timer_ABA = []; Timer_ABA_iLQR=[];
Iters_ABA = []; Iters_ABA_iLQR=[]; 
%ver == 2; Casadi does the derivs for you
for iNb = nbd
    fprintf('\ABA, Number of body %i \n',iNb)
    timer = zeros(1,nmbRepts); iterI= timer;
    iterD = timer; timer_iLQR =timer;
    for iRepts=1:nmbRepts
        Out= DDP_Regular_ABA_Methods(iNb,2,0,N,rbtNmber,x0_diff); %DDP
        OutiLQR= DDP_Regular_ABA_Methods(iNb,2,1,N,rbtNmber,x0_diff);%iLQR
        timer(iRepts)= Out.Time;
        timer_iLQR(iRepts)= OutiLQR.Time;
        iterI(iRepts) = OutiLQR.Iters;
        iterD(iRepts)= Out.Iters;
    end
    VstoreABA = Out.Vstore;
    VstoreABA_iLQR = OutiLQR.Vstore;
    
    Timer_ABA(:,end+1) = mean(timer);
    Timer_ABA_iLQR(:,end+1) = mean(timer_iLQR);
    
    Iters_ABA(:,end+1) = mean(iterD);
    Iters_ABA_iLQR(end+1) = mean(iterI);
end
%}
%%
1==1;

Timer_RNEA =[]; Timer_RNEA_iLQR=[];
Iters_RNEA =[]; Iters_RNEA_iLQR =[];
TimerRNEA_Carp=[]; ItersRNEA_Carp=[];
%RNEA

modRnea=1;%modified RNEA
notModRnea=0;
%RNEA :: DDP_RegularRNEA(Nb,iLQR,N,modRNEA,rbtNmber,x0_diff)
for iNb = nbd
    fprintf('\RNEA, Number of body %i \n',iNb)
    timer = zeros(1,nmbRepts); iterD = timer;
    timer_iLQR =timer; iterI = timer;
    timerCarp = timer; iterCarp = timer; 
    for iRepts=1:nmbRepts
        Out= DDP_RegularRNEA(iNb,0,N,modRnea,rbtNmber,x0_diff);%relations of modRNEA
        OutCarp = DDP_RegularRNEA(iNb,0,N,notModRnea,rbtNmber,x0_diff); %$Casadi does derivs of modRNEA
        OutiLQR = DDP_RegularRNEA(iNb,1,N,modRnea,rbtNmber,x0_diff); %same for both mod/notmod Rnea
        timer(iRepts)= Out.Time;
        timerCarp(iRepts)=OutCarp.Time;
        timer_iLQR(iRepts)= OutiLQR.Time;
        
        iterD(iRepts)= Out.Iters;
        iterCarp(iRepts) =OutCarp.Iters;
        iterI(iRepts)= OutiLQR.Iters;
    end
    VstoreRNEA  = Out.Vstore;
    VstoreRNEACarp  = OutCarp.Vstore;
    VstoreRNEA_iLQR = OutiLQR.Vstore;
    
    Timer_RNEA(:,end+1) = mean(timer);%modRNEA
    TimerRNEA_Carp(:,end+1)= mean(timerCarp);
    Timer_RNEA_iLQR(:,end+1) = mean(timer_iLQR);
    
    Iters_RNEA(:,end+1) = mean(iterD);
    Iters_RNEA_iLQR(end+1) = mean(iterI);
    ItersRNEA_Carp(end+1)= mean(iterCarp);
end

%%
Old_Time = []; Old_Iters =[];
Old_Time_iLQR=[];Old_Iters_iLQR =[];
%
%OldMethod ::  DDP_Regular_OldMethod(Nb,iLQR,N,rbtNmber,x0_diff) 
for iNb =nbd
    fprintf('\tOldWay, Number of body %i \n',iNb)
    timer = zeros(1,nmbRepts); iterD=zeros(1,nmbRepts);
    timer_iLQR = timer; iterI= iterD; 
    for iRepts=1:nmbRepts
        Out= DDP_Regular_OldMethod(iNb,0,N,rbtNmber,x0_diff); %DDP
%         OutiLQR =DDP_Regular_OldMethod(iNb,1,N,rbtNmber,x0_diff);  %iLQR
        timer(iRepts)= Out.Time;
%         timer_iLQR(iRepts)= OutiLQR.Time;
        iterD(iRepts)= Out.Iters; 
%         iterI(iRepts)= OutiLQR.Iters;
    end
    VstoreOld_DDP = Out.Vstore; 
%     VstoreOld_iLQR = OutiLQR.Vstore;
    
    Old_Time(end+1) = mean(timer);
    Old_Time_iLQR(end+1) = mean(timer_iLQR);
    
    Old_Iters(end+1) = mean(iterD);
    Old_Iters_iLQR(end+1) = mean(iterI);
end
%}
%%
1==1;
% Some plots

nbd = [2 4 6 10 15 19 20];
figure;
hold on;
set(gca, 'YScale','log');

DE = semilogy(VstoreABA-VstoreABA(end),'b-o','DisplayName','DDP via ABA'); DE.LineWidth = 2.5;
CD = semilogy(VstoreRNEA-VstoreRNEA(end),'k-o','DisplayName','DDP via modified RNEA'); CD.LineWidth = 2.5;
DC = semilogy(VstoreRNEACarp-VstoreRNEACarp(end),'m-o','DisplayName','DDP via RNEA (Carpentier)'); DC.LineWidth = 2.5;
EF = semilogy(VstoreOld_DDP-VstoreOld_DDP(end),'r-o','DisplayName','DDP via Tensor Contraction'); EF.LineWidth = 2.5;

SF =semilogy(VstoreABA_iLQR-VstoreABA_iLQR(end),'b-.s','DisplayName','iLQR via ABA'); SF.LineWidth = 3.5;
AB = semilogy(VstoreRNEA_iLQR-VstoreRNEA_iLQR(end),'k-.s','DisplayName','iLQR via RNEA'); AB.LineWidth = 2.5;
% AB.Color =[.4660 .6740 .1880];
% FG = semilogy(VstoreOld_iLQR-VstoreOld_iLQR(end),'r-.s','DisplayName','iLQR via Old Method (ABA)'); FG.LineWidth = 2.5;

grid on
hold off 
set(gca,'Fontsize',15)
ld = legend; ld.FontSize = 20; legend
xlabel('Iterations');
ylabel('Suboptimality');

if rbtNmber==1
    title(strcat('Acrobot: Suboptimality Decay vs Iterations for ',{' '},string(nbd(end)),' Links'));
else
    title(strcat('Cartpole: Suboptimality Decay vs Iterations for ',{' '},string(nbd(end)),' Links'));
end 

figure; 
%DDP
hold on
A0=loglog(nbd,Timer_ABA,'b-o'); hold on
A1=loglog(nbd,Timer_RNEA,'k-o'); %Mod RNEA
A1_2=loglog(nbd,TimerRNEA_Carp,'m-o');%Not Mod RNEA
A2=loglog(nbd,Old_Time,'r-o');
%iLQR
A3=loglog(nbd,Timer_ABA_iLQR,'b-.s');
A4=loglog(nbd,Timer_RNEA_iLQR,'k-.s');
% A5=loglog(nbd,Old_Time_iLQR,'r-.s');

if rbtNmber==1
    title('Acrobot: DDP Evaluation Time')
else
    title('Cartpole: DDP Evaluation Time')
end

%DDP
A0.DisplayName = 'DDP via ABA'; A0.LineWidth=2.5;A0.Color='b';
A1.DisplayName = 'DDP via RNEA'; A1.LineWidth=2.5;A1.Color='k';
A1_2.DisplayName = 'DDP via Modified RNEA'; A1_2.LineWidth=2.5;A1_2.Color='m';
A2.DisplayName = 'DDP via Tensor Contraction'; A2.LineWidth=2.5;A2.Color='r';

%iLQR
A3.DisplayName = 'iLQR via ABA'; A3.LineWidth=2.5;%A3.Color=[.4660 .6740 .1880];
A4.DisplayName = 'iLQR via RNEA'; A4.LineWidth=2.5;
% A5.DisplayName = 'iLQR via Old Method (ABA)'; A5.LineWidth=2.5;

lgd = legend; lgd.FontSize=20; legend; 
xlabel('Dofs'); ylabel('Time (s)');
set(gca,'Fontsize',15)
% set(gca, 'YScale', 'log')
xlim([nbd(1) nbd(end)]);
xticks(nbd);
grid on
set(gca,'Yscale','linear')
%%

figure; 
bp=boxplot([Timer_RNEA',Timer_ABA',TimerRNEA_Carp',Old_Time',...
    Timer_ABA_iLQR',Timer_RNEA_iLQR'],'labels',{A0.DisplayName,...
    A1.DisplayName,A1_2.DisplayName, A2.DisplayName,...
    A3.DisplayName,A4.DisplayName} );
xtickangle(20)
set(gca, 'YScale', 'log')
set(gca,'FontSize',12);
set(bp,'LineWidth', 2);
title('Evaluation Time for $n=2$ to $n=20$ Acrobot','Interpreter','latex')
ylabel('Time')



figure; 
bp=boxplot([Iters_RNEA',Iters_ABA',ItersRNEA_Carp',Old_Iters',...
    Iters_ABA_iLQR',Iters_RNEA_iLQR'],'labels',{A0.DisplayName,...
    A1.DisplayName,A1_2.DisplayName, A2.DisplayName,...
    A3.DisplayName,A4.DisplayName} );
xtickangle(20)
set(gca, 'YScale', 'log')
set(gca,'FontSize',12);
set(bp,'LineWidth', 2);
% title('Evaluation Time for $n=2$ to $n=20$ Acrobot','Interpreter','latex')
ylabel('Iters')




%%
h=gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'LastRun_N7BoxPlotTime','-dpdf','-r0')

%%
f

