% Script for varying the initial condition 
% Ensure that guess control is as desired 

%Code is modular 
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

%% Make sure scripts run want you want!!! check 
%% DDP_Regular_ABA_Methods(); DDP_RegularRNEA(); DDP_Regular_OldMethod();
%% Comparing Varying the initial condition 

a = -1; 
b =1; 
N=400;

rbtNmber =1;
nbd =7; 
numRepts = 20;

%row 1 is timer row 2 is iters
ABA_store =zeros(2,numRepts);ABA_iLQR_store=ABA_store; 
RNEA_store =zeros(2,numRepts);RNEA_iLQR_store=RNEA_store; 
modRNEA_store=zeros(2,numRepts);
Old_store = zeros(2,numRepts);

for idx =1:numRepts
    x0_diff = (b-a).*rand(nbd*2,1) + a; %x0_diff=x0_diff;
    fprintf('\nInitial x0 = [pi/2 + %.3f;0;0...]\n',x0_diff);
 
    
%     x0_diff = [pi/2*ones(1,7) zeros(1,7)]; x0_diff(1) = x0_diff(1) +pi/2;
    
    1==1;
    ABA= DDP_Regular_ABA_Methods(nbd,1,0,N,rbtNmber,x0_diff);
    ABA_store(:,idx) = [ABA.Time;ABA.Iters];
    
    ABAiLQR= DDP_Regular_ABA_Methods(nbd,1,1,N,rbtNmber,x0_diff);
    ABA_iLQR_store(:,idx) = [ABAiLQR.Time; ABAiLQR.Iters];

    
    RNEA= DDP_RegularRNEA(nbd,0,N,1,rbtNmber,x0_diff);
    RNEA_store(:,idx)=[RNEA.Time; RNEA.Iters];
    
    modRNEA= DDP_RegularRNEA(nbd,0,N,0,rbtNmber,x0_diff);
    modRNEA_store(:,idx)=[modRNEA.Time; modRNEA.Iters];

    RNEAiLQR = DDP_RegularRNEA(nbd,1,N,0,rbtNmber,x0_diff);
    RNEA_iLQR_store(:,idx) = [RNEAiLQR.Time; RNEAiLQR.Iters];
    

    Old= DDP_Regular_OldMethod(nbd,0,N,rbtNmber,x0_diff); %DDP
    Old_store(:,idx) = [Old.Time;Old.Iters];
 
end

save('VaryInitConditions')

labels_iLQR ={'iLQR via ABA', 'iLQR via RNEA'};
labels_DDP ={'DDP via ABA', 'DDP via RNEA', 'DDP via Modified RNEA',...
    'DDP via Tensor Contraction'};
labels = [labels_iLQR labels_DDP]; 

figure;
row = 1;%Time
Tme = boxplot([RNEA_iLQR_store(row,:)',ABA_iLQR_store(row,:)',...
    ABA_store(row,:)',RNEA_store(row,:)',modRNEA_store(row,:)',Old_store(row,:)'],...
    'labels', labels); 

xtickangle(20);
set(gca, 'YScale', 'log'); 
set(gca,'FontSize',12);
set(Tme,'LineWidth', 2);
L = gca;
L.XAxis.TickLabelInterpreter = 'latex';
ylabel('Logarithm of Time(s)','Interpreter','latex')


figure; 
row =2;%Iters
Itr = boxplot([RNEA_iLQR_store(row,:)',ABA_iLQR_store(row,:)',...
    ABA_store(row,:)',RNEA_store(row,:)',modRNEA_store(row,:)',Old_store(row,:)'],...
    'labels',labels); 
xtickangle(20);
set(gca,'FontSize',12);
set(Itr,'LineWidth', 2);
L = gca;
L.XAxis.TickLabelInterpreter = 'latex';
ylabel('Iterations','Interpreter','latex')

% Adjust these values to match the plots
txt1 = {'First-order'};
t1= text(1,22.5,txt1);
t1.Interpreter ='latex';
t1.FontSize = 12;

txt2 = {'Second-order Information'};
t2= text(2.5,7.5,txt2);
t2.Interpreter ='latex';
t2.FontSize = 12;

hold on
%FaceAlpha = 0.5 is light. 1 is dark
h = findobj(gca,'Tag','Box');
colors = {'r','k','m','b','k','b'}
FaceAl = [1;1;1;1;0.5;0.5]
% FaceAl = [0.5;0.5;1;1;1;1]
for i=1:length(h)
    bb(i) = patch(get(h(i),'XData'),get(h(i),'YData'),colors{i},'FaceAlpha',FaceAl(i));
end

b2b = plot(2.5*ones(1,25),0:24);
b2b.LineWidth = 2;
b2b.Color = 'k';

b3b = plot(2.5*ones(1,67),4:70) 
b3b.Color = 'k'
b3b.LineWidth = 2;

hold off