clear 
addpath([pwd '/spatial_v2_casadi']);
addpath(genpath([pwd])); 

import casadi.*

% Loop through kinematics chains of different sizes

Nb_lower = 2;         % Lower bound on sweep over # bodies
Nb_upper = 40;        % Upper bound on sweep over # bodies
sweep_points = 20;    % # of points in the sweep
second_order = 1;     % 0=1st order partials only, 1=2nd order too
direction    = 0;     % +1=forward mode AD, -1=reverse mode AD, 0=CasADi choice

% Model parameters
branching_factor = 1; %1=chain, 2=binary tree
skewing = pi/exp(1);  %0=2D mechanism, anything else=3D mechanism
    
ii = 0;
Ts = [];


Nb_lst = unique(round(logspace(log10(Nb_lower),log10(Nb_upper),sweep_points)));
%Nb_lst = 5;
for Nb = Nb_lst
    fprintf('Nb = %d \n',Nb);
    
    % Generate the model
    model = autoTree( Nb, branching_factor, skewing);
    
    for i = 1:Nb
        Inertia = rand(3);
        model.I{i}(1:3,1:3) = Inertia'*Inertia;
    end
    
    fprintf('\tComputing Direct Functions (Derivs of ABA)\n');
    ABA_Deriv_Funcs = GetDerivativeFunctionsDirect(model,second_order,direction);
    fprintf('\tComputing Alt Functions (via Derivs of RNEA)\n');
    ABA_Deriv_via_RNEA_Funcs    = GetDerivativeFunctionsAlt_mod(model,second_order,direction);
    fprintf('\tComputing Old Functions\n');
    funcsOld = GetDerivativeFunctionsOld(model,1);

    % Values for numericallly evaluating the resulting functions
    q_val = rand(Nb,1);
    qd_val = rand(Nb,1);
    tau_val = rand(Nb-1,1);
    qdd_val = FDab(model, q_val, qd_val, [tau_val;0]);
    eta_val = rand(Nb,1); % qdot part of Vx in DDP
    
    sz = Nb * [1 1 1]; %Size of reshape of second derivs 
    sztau = sz; sztau(end)= sz(end)-1;

    % Optional, check all of the methods for computing the partials to make sure they agree
    CheckPartials = 1;   
    if CheckPartials
        disp('==========================')
        disp('First order partials of q')
        disp('==========================')

        Jq_1 = full(ABA_Deriv_Funcs.q(q_val,qd_val,tau_val));
        Jq_2 = full(ABA_Deriv_via_RNEA_Funcs.q(q_val,qd_val,qdd_val));
        err = norm(Jq_1-Jq_2)
        assert( err < 1e-6)

        disp('==========================')
        disp('First order partials of qd')
        disp('==========================')

        Jqd1 = full(ABA_Deriv_Funcs.qd(q_val,qd_val));
        Jqd2 =full(ABA_Deriv_via_RNEA_Funcs.qd(q_val,qd_val));
        err = norm(Jqd1-Jqd2)
        assert( err < 1e-6)

        disp('==========================')
        disp('First order partials of tau')
        disp('==========================')

        Jtau1 = full(ABA_Deriv_Funcs.tau(q_val));
        Jtau2 = full(   ABA_Deriv_via_RNEA_Funcs.tau(q_val));
        err = norm(Jtau1-Jtau2(:,1:end-1))
        assert( err < 1e-6)
        
        
%         kk =ABA_Deriv_via_RNEA_Funcs.all_first(q_val,qd_val,tau_val);

        disp('==========================')
        disp('Second order partials of q')
        disp('==========================')

        mu_val = ABA_Deriv_via_RNEA_Funcs.InvM_Mult(q_val,eta_val);
        
        [H_mod, J_mod] = ABA_Deriv_via_RNEA_Funcs.mod_second(q_val, qd_val, qdd_val, mu_val, Jq_1, Jqd1, Jtau2);
        
        %H_mod = full( ABA_Deriv_via_RNEA_Funcs.mod_hessian(q_val, qd_val, qdd_val, mu_val) );
        %J_mod = full( ABA_Deriv_via_RNEA_Funcs.mod_jac(q_val, qd_val, qdd_val, mu_val, Jq_1, Jqd1, Jtau1) );
        
        A = - 1/2*H_mod(1:Nb,1:Nb) - J_mod(1:Nb,:);
        H_qq_mod = full(A+A');
        H_qd_qd_mod = full(-H_mod(Nb+1:end, Nb+1:end));
        H_q_qd_mod = full(-H_mod(1:Nb, Nb+1:end) - J_mod(Nb+1:2*Nb,:)');
        H_q_tau = full(-J_mod(2*Nb+1:end,:)');
        
        

        a1 = full( ABA_Deriv_Funcs.q_q(q_val ,qd_val, tau_val, eta_val) );    
        a2 = full(    ABA_Deriv_via_RNEA_Funcs.q_q(q_val ,qd_val, qdd_val, mu_val, Jq_1) );
        err = norm(a1-a2)
        
        e_qq= norm(H_qq_mod-a1)
        
        assert( err < 1e-6)
        assert( e_qq < 1e-6)

        disp('===========================')
        disp('Second order partials of qd')
        disp('===========================')

        a1 = full( ABA_Deriv_Funcs.qd_qd(q_val, qd_val,eta_val));
        a2 =  full( ABA_Deriv_via_RNEA_Funcs.qd_qd(q_val, qd_val, mu_val) );
        err = norm(a1-a2)
        assert( err < 1e-6)
        
        
        
        e_qd_qd = norm(H_qd_qd_mod-a2)
        assert( e_qd_qd < 1e-6)


        disp('===================================')
        disp('Mixed Second order partials of q,qd')
        disp('===================================')

        a1 =  full( ABA_Deriv_Funcs.q_qd(q_val, qd_val,eta_val) );
        a2 = full(ABA_Deriv_via_RNEA_Funcs.q_qd(q_val,qd_val,mu_val,Jqd1));
        err = norm(a1-a2)
        assert( err < 1e-6)
        
        e_q_qd = norm(H_q_qd_mod - a1)
        
        assert(e_q_qd < 1e-6);

        disp('=====================================')
        disp('Mixed Second order partials of q,tau')
        disp('=====================================')

        a1 = full( ABA_Deriv_Funcs.q_tau(q_val, tau_val,eta_val) );
        a2 = full(    ABA_Deriv_via_RNEA_Funcs.q_tau(q_val, qdd_val,mu_val,Jtau2));
        err = norm(a1-a2(:,1:end-1))
        assert( err < 1e-6)
        
        
        
        
        
        e_q_tau = norm(H_q_tau(:,1:end-1)-a1)
        assert( e_q_tau < 1e-6)

        
        %return
    end


    fprintf('\tTesting Timing\n');
    % Timing
    t_direct_first = []; % timing for direct first-order derivatives of ABA
    t_direct_second=[];  % timing for direct second-order  "" ""
    t_direct_second_v2 = []; % Alternate strategy for the second order derivatives
    t_alt_first = [];    % timing for first-order derivatives of ABA via RNEA
    t_alt_second = [];   % timing for second-order derivatives of ABA via RNEA
    t_mod_second = [];   % Method using modified RNEA
    t_old_first = [];  %timing first-order derivs of tensor contraction
    t_old_second = []; %timing second-order derivs of tensor contraction
    
    for i = 1:10
        tic
        f1_all = ABA_Deriv_Funcs.all_first(q_val,qd_val,tau_val);
        t_direct_first(end+1) = toc;

        tic
        f2_all = ABA_Deriv_via_RNEA_Funcs.all_first(q_val,qd_val,qdd_val);
        t_alt_first(end+1) = toc;
        
        tic 
        f3_all = funcsOld.all_first(q_val,qd_val,tau_val);  
        t_old_first(end+1)=toc; 
        
        if second_order
            tic;
            [first_derivs,second_derivs,H1_q_tau] = ABA_Deriv_Funcs.All_second(q_val,qd_val, tau_val, eta_val);
            t_direct_second(end+1) = toc;

            tic;
            [first_derivs,second_derivs,H1_q_tau] = ABA_Deriv_Funcs.All_second_V2(q_val,qd_val, tau_val, eta_val);
            t_direct_second_v2(end+1) = toc;

            tic; 
            [first_derivs,second_derivs,H1_q_tau] = ABA_Deriv_via_RNEA_Funcs.All_second(q_val,qd_val,[tau_val;0], eta_val);%RNEA
            t_alt_second(end+1) = toc;
            
        
            
            tic
            [first_derivs,second_derivs,H1_q_tau] = ABA_Deriv_via_RNEA_Funcs.mod_all(q_val, qd_val,[tau_val;0],eta_val);%mod_RNEA
            t_mod_second(end+1) = toc;

            1==1;

            tic 
            [first_derivs, qq,qqd,qtau,qdqd]=funcsOld.all(q_val,qd_val,tau_val); 
            M11 = reshape(full(qq),sz); qq = Tens3byVec(M11,eta_val','pre');
            M12 = reshape(full(qqd),sz); qqd = Tens3byVec(M12,eta_val','pre');            
            M15 = reshape(full(qtau),sztau); hux = permute(M15,[2 3 1]); 
            qtau = Tens3byVec(hux,eta_val','pre');
            M14 = reshape(full(qdqd),sz); qdqd = Tens3byVec(M14,eta_val','pre');
            t_old_second(end+1)=toc;
            
        end



    end
    ii = ii+1;
	% Collecing the Time Spent for each algo
    Ts(1,ii) = mean(t_direct_first);
    Ts(2,ii) = mean(t_alt_first);
    if second_order
        Ts(3,ii) = mean(t_direct_second);
        Ts(4,ii) = mean(t_alt_second);
        Ts(5,ii) = mean(t_direct_second_v2);
        
        Ts(6,ii) = mean(t_mod_second);
%         Ts(7,ii) = mean(t_old_first); %It's not different from Ts(1)
        Ts(7,ii)= mean(t_old_second);
    end
    
end

%%
figure(1);
clf
%subplot(211)
leg = {};
h{1} = loglog(Nb_lst, Ts(1,:)); hold on;
h{2} = loglog(Nb_lst, Ts(2,:));
% h{3} = loglog(Nb_lst, Ts(7,:));  %fIND ME
if second_order
    h{3} = loglog(Nb_lst, Ts(3,:)); 
    h{4} = loglog(Nb_lst, Ts(4,:));
    h{5} = loglog(Nb_lst, Ts(5,:));
    h{6} = loglog(Nb_lst, Ts(6,:));
    h{7} = loglog(Nb_lst, Ts(7,:));
end


leg{1} = 'First-Order Partials of $\ddot{\mathbf{q}}$ via ABA (directly)';
leg{2} = 'First-Order Partials of $\ddot{\mathbf{q}}$ via RNEA';
% leg{8} 
if second_order
    leg{3} = 'Second-Order partials of ${\eta}^T\ddot{\mathbf{q}}$ via ABA (Direct)'; %find me
    leg{4} = 'Second-Order Partials of ${\eta}^T\ddot{\mathbf{q}}$ via RNEA';
    leg{5} = 'Second-Order Partials of ${\eta}^T\ddot{\mathbf{q}}$ via ABA (direct, alternate strategy)';
    leg{6} = 'Second-Order Partials of ${\eta}^T\ddot{\mathbf{q}}$ via mod RNEA'; 
    leg{7} = 'Second-Order Partials of ${\eta}^T\ddot{\mathbf{q}}$ via Tensor Contraction';
end

for i = 1:length(h)
    h{i}.LineWidth = 3;
    p{i} = polyfit(log10(h{i}.XData(h{i}.XData>10)),log10(h{i}.YData(h{i}.XData>10)),1);
    pi = p{i};
    
    h_fit{i}= loglog(Nb_lst, 10^pi(2) * Nb_lst.^pi(1));
    h_fit{i}.LineWidth = 2;
    h_fit{i}.Color = h{i}.Color;
    h_fit{i}.LineStyle = ":";
    
    leg{end+1} = ['$' sprintf('{\\rm log}_{10}(t) = %.3f \\, {\\rm log}_{10}(N_b) %.3f   )',pi(1), pi(2)) '$'];
end

grid on
l= legend([h{:} h_fit{:}], leg);%,);
l.Interpreter = 'Latex';
xlabel('Dofs');
ylabel('Time (s)');
g = gca;
set(g,'FontSize',20)
grid on;
xlim([Nb_lst(1) Nb_lst(end)]);

1==1;
figure; hold on; grid on;
for idx =1:length(h)
    b{idx}=plot(Nb_lst,Ts(idx,:));
    b{idx}.Color = h{idx}.Color;
    b{idx}.LineWidth=2;
end
xlabel('Dofs');
ylabel('Time (s)');
grid on;
lgd=legend([leg(1:length(h))]); lgd.Interpreter = 'Latex';
lgd.Location= 'best';
g = gca; set(g,'FontSize',20)
xlim([Nb_lst(1) Nb_lst(end)]);
% 
% subplot(212)
% semilogx(Nb_lst, Ts(1,:)./Ts(2,:))
% xlabel('DoFs');
% ylabel('Speedup of using Partials of RNEA')
% g = gca;
% set(g,'FontSize',20)
% grid on;

