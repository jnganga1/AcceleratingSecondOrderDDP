

%Partials





% nbd = [2 4 6 10 15 19 20]; % 25 30];
nbd = [5 6 7] ;
nmbRepts=2;

%ABA::  DDP_Regular_ABA_Methods(Nb,ver,iLQR,N,rbtNmber,x0_diff)
Timer_ABA_v1 = []; Timer_ABA_v2 = []; Timer_ABA_iLQR=[];
Iters_ABA_v1 = []; Iters_ABA_v2 = []; Iters_ABA_iLQR=[]; 

Timer_MODRNEA =[]; Timer_NOTMODRNEA =[]; Timer_RNEA_ILQR =[];
iters_MODRNEA =[]; iters_NOTMODRNEA =[]; iters_RNEA_ILQR =[];

Timer_Old = []; 
Iters_Old = [];

x0_diff = 0.15;rbtNmber=1;N = 400;
%ver == 2; Casadi does the derivs for you
for iNb = nbd
    fprintf('Number of body %i \n',iNb)
    timer_v1_ABA = zeros(1,nmbRepts); iter_v1_ABA= zeros(1,nmbRepts);
    timer_v2_ABA = zeros(1,nmbRepts); iter_v2_ABA = zeros(1,nmbRepts);
    Timer_iLQR = zeros(1,nmbRepts);   iter_iLQR = zeros(1,nmbRepts); 

    timer_notModRNEA = zeros(1,nmbRepts); iters_notModRNEA =  zeros(1,nmbRepts);
    timer_modRNEA = zeros(1,nmbRepts); iter_modRNEA = zeros(1,nmbRepts);
    timer_RNEA_ilqr = zeros(1,nmbRepts); Iters_RNEA_iLQR = zeros(1,nmbRepts);
    
    timer_tens  = zeros(1,nmbRepts);
    iters_tens   = zeros(1,nmbRepts);
    for iRepts=1:nmbRepts
        Out_v1= DDP_Regular_ABA_Methods(iNb,1,0,N,rbtNmber,x0_diff); %v1 ABA
        Out_v2= DDP_Regular_ABA_Methods(iNb,2,0,N,rbtNmber,x0_diff); %v2 ABA
        Out_ABA_iLQR= DDP_Regular_ABA_Methods(iNb,2,1,N,rbtNmber,x0_diff);%ABA iLQR
        
        timer_v1_ABA(iRepts)= Out_v1.Time;
        timer_v2_ABA(iRepts)= Out_v2.Time;
        Timer_iLQR(iRepts)= Out_ABA_iLQR.Time;
        
        iter_v1_ABA(iRepts) = Out_v1.Iters;
        iter_v2_ABA(iRepts) = Out_v2.Iters;
        iter_iLQR(iRepts)= Out_ABA_iLQR.Iters;
        
        modRnea=1; notModRnea = 0; 
        Out_notMod = DDP_RegularRNEA(iNb,0,N,notModRnea,rbtNmber,x0_diff); %$Casadi does derivs of modRNEA
        Out_modRNEA= DDP_RegularRNEA(iNb,0,N,modRnea,rbtNmber,x0_diff);%relations of modRNEA
        Out_Rnea_iLQR = DDP_RegularRNEA(iNb,1,N,modRnea,rbtNmber,x0_diff); %same for both mod/notmod Rnea
        
        timer_modRNEA(iRepts)= Out_modRNEA.Time;
        timer_notModRNEA(iRepts)=Out_notMod.Time;
        timer_RNEA_ilqr(iRepts)= Out_Rnea_iLQR.Time;
        
        iter_modRNEA(iRepts)= Out_modRNEA.Iters;
        iters_notModRNEA(iRepts) =Out_notMod.Iters;
        Iters_RNEA_iLQR(iRepts)= Out_Rnea_iLQR.Iters;
        
        Out_tens= DDP_Regular_OldMethod(iNb,0,N,rbtNmber,x0_diff);%DDP Tensor
        
        timer_tens(iRepts) = Out_tens.Time; 
        iters_tens(iRepts) = Out_tens.Iters;
    end
        
    Timer_ABA_v1(end+1) = mean(timer_v1_ABA);
    Timer_ABA_v2(end+1) = mean(timer_v2_ABA);
    Timer_ABA_iLQR(end+1) = mean(Timer_iLQR);
    
    Iters_ABA_v1(end+1) = mean(iter_v1_ABA);
    Iters_ABA_v2(end+1) = mean(iter_v2_ABA);
    Iters_ABA_iLQR(end+1) = mean(iter_iLQR);
    
    
    Timer_MODRNEA(end+1) = mean(timer_modRNEA);%modRNEA
    Timer_NOTMODRNEA(end+1)= mean(timer_notModRNEA);
    Timer_RNEA_ILQR(end+1) = mean(timer_RNEA_ilqr);
    
    iters_MODRNEA(end+1) = mean(iter_modRNEA);
    iters_NOTMODRNEA(end+1) = mean(iters_notModRNEA);
    iters_RNEA_ILQR(end+1)= mean(Iters_RNEA_iLQR);
    
    Timer_Old(end+1) = mean(timer_tens);
    Iters_Old(end+1) = mean(iters_tens);

end


%%

%Did it this way since struct sometimes crash
save('RepeatABARNEA.mat',...
    'Timer_ABA_v1','Timer_ABA_v2','Timer_ABA_iLQR',...
    'Timer_MODRNEA','Timer_NOTMODRNEA','Timer_ABA_iLQR',...
    'Iters_ABA_v1','Iters_ABA_v2','Iters_ABA_iLQR',...
    'iters_MODRNEA','iters_NOTMODRNEA','iters_RNEA_ILQR');



figure; 
plot(nbd,Timer_ABA_v1,'DisplayName','ABA v1'); hold on 
plot(nbd,Timer_ABA_v2,'DisplayName','ABA v2 - casadi grouped') 
plot(nbd,Timer_ABA_iLQR,'DisplayName','ABA iLQR')

plot(nbd,Timer_MODRNEA,'DisplayName','ModRNEA');
plot(nbd,Timer_NOTMODRNEA, 'DisplayName','NotModRNEA');
plot(nbd,Timer_RNEA_ILQR,'DisplayName','RNEA iLQR');
legend


%Randomization 