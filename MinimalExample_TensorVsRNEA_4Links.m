
Nb = 2; 

%Select robot here - ensures all code runs the same thing
rbtNmber = 1; %pendubot model
x0_diff = 0.15;

N = 400;

numrepts = 1; 

RNEA_Timer = zeros(1,numrepts); 
Tensor_Timer = zeros(1,numrepts); 
global TensorBackTime TensorBackIters
global RNEAbackTime RNEAbackIters

RNEAbackTime = 0;
RNEAbackIters = 0;
TensorBackTime = 0;
TensorBackIters = 0;

for idx = 1:numrepts

    OutCarp = DDP_RegularRNEA(Nb,0,N,1,rbtNmber,x0_diff); %RNEA_DDP: Casadi does derivs of modRNEA
    RNEA_Timer(idx) = OutCarp.Time; 

    Out= DDP_Regular_OldMethod(Nb,0,N,rbtNmber,x0_diff);  %Tensor DDP
    Tensor_Timer(idx) = Out.Time;

end

fprintf('\nDDP via RNEA took: %.2f Seconds\n',mean(RNEA_Timer));
fprintf('\n Qinfo took: %.2f Seconds\n',RNEAbackTime);



fprintf('\nDDP via TensorContraction took: %.2f Seconds\n',mean(Tensor_Timer));
fprintf('\n Qinfo took: %.2f Seconds\n',TensorBackTime);