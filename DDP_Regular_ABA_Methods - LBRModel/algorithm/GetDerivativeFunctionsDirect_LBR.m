function funcs = GetDerivativeFunctionsDirect_LBR(model, second_order, direction)
import casadi.*

    if nargin <=2
        direction  =0;
    end
    if nargin == 1
        second_order = true;
    end
    
    function J = jaco(f, x)
        if direction == 1
            J = jtimes(f,x,eye(length(x)));
        elseif direction == -1
            J = jtimes(f,x,eye(length(f)),true)';
        else
            J = jacobian(f,x);
        end
    end
    
    % Casadi symbolics for usual stuff
    q = MX.sym('q',[model.NB 1]);
    qd = MX.sym('qd',[model.NB 1]);
    eta = MX.sym('eta',[model.NB 1]);
    tau = MX.sym('tau',[model.NB 1]);
    
    % Call articulated body algorithm to get expression for qdd
    qddABA = FDab_casadi(model,q,qd,tau);
    funcs.qddABA = Function('qddABA',{q,qd,tau},{qddABA});

    % Direct Jacobians
    funcs.q = Function('Diff_ABA_q', {q,qd,tau}, {jaco(qddABA, q)} ); % O(N^2) function
    funcs.qd = Function('Diff_ABA_qd',{q,qd}, {jaco(qddABA, qd)} ); % O(N^2) function
    funcs.tau = Function('Diff_ABA_tau',{q}, {jaco(qddABA, tau)} ); % O(N^2) function
    funcs.all_first = Function('Diff_ABA_all',{q,qd,tau},{jaco(qddABA, [q;qd;tau])}); % O(N^2) function

    
   
    if second_order 
        
        % This call uses reverse mode AD to compute the gradient of eta'*qdd
        H_q  = jtimes( qddABA ,q, eta,true );
        H_qd = jtimes( qddABA ,qd, eta,true );
        H_tau = jtimes( qddABA ,tau, eta,true );
        
        % Rather than calculating the first/second order partials separately, this line groups
        % them together so that hopefully casadi can be smarter about it
        %H_grouped = [eta'*qddABA ; H_q; H_qd ; H_tau];
        H_grouped_hess = hessian(eta'*qddABA, [q;qd;tau]);

        % Get the second order partials individually
        H_q_q =  jaco(H_q, q);    % Each column of H_q_q can be eavluated in O(N) time, thus this function takes O(N^2) time
        H_qd_qd = jaco(H_qd, qd); 
        H_q_qd = jaco(H_q, qd);  
        H_q_tau = jaco(H_q, tau);

        % Functions for the second order partials individually
        funcs.q_q = Function('H_ABA_q_q', {q,qd,tau,eta}, {H_q_q} );
        funcs.qd_qd = Function('H_ABA_qd_qd', {q,qd,eta}, {H_qd_qd} );
        funcs.q_qd = Function('H_ABA_q_qd', {q,qd,eta}, {H_q_qd} );
        funcs.q_tau = Function('H_ABA_q_tau', {q,tau,eta}, {H_q_tau} );

        % Functions for the second-order partials all together (also
        % outputs some first-order partials that we don't really need)
        funcs.all_second = Function('H_ABA_all', {q,qd,tau,eta}, {H_q_q,H_qd_qd,H_q_qd',H_q_tau} );
        funcs.all_second_v2 = Function('H_ABA_all_v2',{q,qd,tau,eta},{H_grouped_hess});
    end
return 
end
