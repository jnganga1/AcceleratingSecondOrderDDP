
Nb = 7; 

%Select robot here - ensures all code runs the same thing
rbtNmber = 1; %pendubot model
x0_diff = 0.15;

N = 400;

Nb = [2 4 6 7];

numrepts = length(Nb); 


RNEA_Timer = zeros(1,numrepts); 
Tensor_Timer = zeros(1,numrepts); 
%global TensorBackTime TensorBackIters
% global RNEAbackTime RNEAbackIters

% RNEAbackTime = 0;
% RNEAbackIters = 0;
% TensorBackTime = 0;
% TensorBackIters = 0;
i=0;
ilqr =0;
for idx = Nb
    
    
     i=i+1;
    OutCarp = DDP_RegularRNEA(idx,ilqr,N,0,rbtNmber,x0_diff); %RNEA_DDP: Casadi does derivs of modRNEA
    RNEA_Timer(i) = OutCarp.Time; 

    Out= DDP_Regular_OldMethod(idx,ilqr,N,rbtNmber,x0_diff);  %Tensor DDP
    Tensor_Timer(i) = Out.Time;
    
    OutABA = DDP_Regular_ABA_Methods(idx,1,ilqr,N,rbtNmber,x0_diff) 
    OutABA.Time

end

fprintf('\nDDP via RNEA took: %.2f Seconds\n',mean(RNEA_Timer));
%fprintf('\n Qinfo took: %.2f Seconds\n',RNEAbackTime);



fprintf('\nDDP via TensorContraction took: %.2f Seconds\n',mean(Tensor_Timer));
%fprintf('\n Qinfo took: %.2f Seconds\n',TensorBackTime);

figure; 
plot(Nb,RNEA_Timer,'DisplayName','RNEA');hold on 
plot(Nb,Tensor_Timer,'DisplayName','Tensor'); 
legend

%% For Pat:: ModRNEA vs RNEA
x0_diff = 0.15; N = 400; rbtNmber=1;
ilqr =0;
modRNEA = 1; notModRNEA =0;
RNEA_TimeD =[]; notModRNEA_TimeD=[];
nb= [2 4 6 8 10]
for i=nb
    RNEA_Time=[];notModRNEA_Time=[];
  for idx=1:2
     Out= DDP_RegularRNEA(i,ilqr,N,modRNEA,rbtNmber,x0_diff);%relations of modRNEA
     RNEA_Time(end+1) = Out.Time; 
     OutNotRNEA= DDP_RegularRNEA(i,ilqr,N,notModRNEA,rbtNmber,x0_diff); 
     notModRNEA_Time(end+1) = OutNotRNEA.Time;  
  end
     RNEA_TimeD(end+1) = mean(RNEA_Time);
     notModRNEA_TimeD(end+1) = mean(notModRNEA_Time);
end
figure; 
plot(nb,RNEA_TimeD,'DisplayName','modRNEA');hold on; plot(nb,notModRNEA_TimeD,'DisplayName','notModRNEA')
legend
xlabel('Dofs');ylabel('4 iterations of DDP');

