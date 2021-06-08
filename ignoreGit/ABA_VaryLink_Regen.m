%% varying number of links
%looking at specifics
%This section -> varying number of links
%Compares Suboptimality vs Iterations 
%Compares EvaluationTime vs number of links

%Select robot here - ensures all code runs the same thing
rbtNmber =1;
x0_diff = 0.15;

N = 400;
% nbd = [2 3 4 5 7 9 12 16 22 26 30];
nbd = [2 4 6 10 15 19 20 22]; % 25 30];
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
save('ABA_varyLink_regen')
