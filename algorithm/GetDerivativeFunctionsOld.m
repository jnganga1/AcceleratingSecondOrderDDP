function funcs = GetDerivativeFunctionsOld(model,second_order)
    import casadi.*
    
    % Casadi symbolics
    q = MX.sym('q',[model.NB 1]);
    qd = MX.sym('qd',[model.NB 1]);
    tau = MX.sym('tau',[model.NB-1 1]);
    
     
    % Call articulated body algorithm to get expression for qdd
    qddABA = FDab_casadi(model,q,qd,[tau;0]);
    funcs.qddABA = Function('qddABA',{q,qd,tau},{qddABA});
       
    
    %Some derivs and useful relations 
    jac_q = jacobian(qddABA, q); 
    jac_qd = jacobian(qddABA,qd);
%     jac_x = jacobian(qddABA,[q;qd]);

    
    
    % Direct Jacobians of FD
    funcs.q = Function('Diff_ABA_q', {q,qd,tau}, {jac_q} ); % O(N^2) function
%     funcs.x = Function('Diff_ABA_q', {q,qd,tau}, {jac_x} );
    funcs.qd = Function('Diff_ABA_qd',{q,qd}, {jac_qd} ); % O(N^2) function
    funcs.tau = Function('Diff_ABA_tau',{q}, {jacobian(qddABA, tau)} ); % O(N^2) function
    funcs.all_first = Function('Diff_ABA_all',{q,qd,tau},{jacobian(qddABA, [q;qd;tau])}); % O(N^2) function
    
%      qVal= rand(size(q)); qdVal=rand(size(qd));
%      tauVal= rand(size(tau));
%     qddVal = full(FDab_casadi(model_no_gravity,qVal,qdVal,tauVal))
    
    if second_order              
        reshape_jac_q = reshape(jac_q,[numel(jac_q) ,1]); 
        reshape_jac_qd = reshape(jac_qd,[numel(jac_qd) ,1]);
    
        By_qq = jacobian(reshape_jac_q ,q); 
        By_qqd =jacobian(reshape_jac_q,qd);
        By_qtau =jacobian(reshape_jac_q,tau);
        By_qdqd = jacobian(reshape_jac_qd,qd);  
        
        funcs.qq = Function('Hess_ABA_qq',{q,qd,tau},{By_qq});
        funcs.qqd = Function('Hess_ABA_qqd',{q,qd},{By_qqd});
        funcs.qtau = Function('Hess_ABA_qqd',{q,tau},{By_qtau});
        funcs.qdqd = Function('Hess_ABA_q',{q,qd},{By_qdqd});
        funcs.all_second = Function('Hess_ABA_all',{q,qd,tau},{By_qq,By_qqd,By_qtau,By_qdqd});
        
        
        funcs.all = Function('Diff_Hess_ABA_all',{q,qd,tau},{jacobian(qddABA, [q;qd;tau]), By_qq,By_qqd,By_qtau,By_qdqd }); % O(N^2) function
        
        1==1;
    end 
end