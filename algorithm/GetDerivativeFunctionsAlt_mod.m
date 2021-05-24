function funcs = GetDerivativeFunctionsAlt_mod(model, second_order , direction)
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

    
    % Casadi symbolics
    q = MX.sym('q',[model.NB 1]);
    qd = MX.sym('qd',[model.NB 1]);
    qdd = MX.sym('qdd',[model.NB 1]);    
    tau = MX.sym('tau',[model.NB 1]);
    
    % Call articulated body algorithm to get expression for qdd
    tauH=MX.sym('tau',[model.NB-1 1]);
    qddABA = FDab_casadi(model,q,qd,[tauH;0]);
    funcs.qddABA = Function('qddABA',{q,qd,tauH},{qddABA});
        
    
    % Get expression for the inverse dynamics model
    tauID = ID_casadi(model,q,qd,qdd);
    
    % We'll use a model with no gravity to compute M*nu for any vector nu
    model_no_gravity = model;
    model_no_gravity.gravity = [0;0;0]';
    
    % Function for multiplying by M (using RNEA) [O(N)]
    M_mult = Function('M_mult', {q, qd}, {ID_casadi(model_no_gravity, q, 0*qd, qd)});
    
    % Function for multiplying by M inverse (using ABA) [O(N)]
    InvM_Mult = Function('InvM_Mult',{q,tau}, {FDab_casadi(model_no_gravity, q, qd*0, tau)}); % O(N) function
    
    % Partials of RNEA
    Diff_ID_q = Function('Diff_ID_q',{q,qd,qdd}, {jaco(tauID, q)}); % O(N^2) function
    Diff_ID_qd = Function('Diff_ID_qd',{q,qd}, {jaco(tauID, qd)}); % O(N^2) function
    
    % All together
    Diff_ID_all = Function('Diff_ID_all',{q,qd,qdd}, {[jaco(tauID, [q;qd]) -MX.eye(model.NB)] }); % O(N^2) function

    
    funcs.M_mult = M_mult;
    funcs.InvM_Mult = InvM_Mult;
    
    % Apply the identify to relate partials of InvDyn to partials of Fwd Dyn
    funcs.q  = Function('Diff_FD_q',{q,qd,qdd},{-InvM_Mult(q, Diff_ID_q(q,qd,qdd))}); % O(N^2) function
    funcs.qd  = Function('Diff_FD_q',{q,qd},{-InvM_Mult(q, Diff_ID_qd(q,qd))}); % O(N^2) function
    funcs.tau = Function('Diff_FD_tau',{q},{InvM_Mult(q, MX.eye(model.NB) )}); % O(N^2) function
    
    % Group all first order partials together to see if Casadi can do it faster in this case
    funcs.all_first   = Function('Diff_FD_all',{q,qd,qdd},{-InvM_Mult(q, Diff_ID_all(q,qd,qdd))});

    Nb = model.NB;
    if second_order
        % More casdi symbolics
        mu = MX.sym('mu',Nb);
        lambda = MX.sym('lambda',Nb);
        
        nu = MX.sym('nu',[Nb,Nb]);
        nu_q = MX.sym('nu_q',[Nb,Nb]);
        nu_qd = MX.sym('nu_qd',[Nb,Nb]);
        nu_tau = MX.sym('nu_tau',[Nb,Nb]);
        
        out = modID_casadi(model,q,qd,qdd,mu);
        out_big = modID_casadi(model_no_gravity,q,0*qd,mu,[nu_q nu_qd nu_tau]);
        
                
        full_hessian_mod = hessian(out,[q;qd]);
        jac_big_mod      = jacobian(out_big, q);
        
        funcs.mod_hessian = Function('Mod_all_all',{q,qd,qdd,mu},{full_hessian_mod});
        funcs.mod_jac = Function('Mod_jac',{q,qd,qdd,mu,nu_q,nu_qd,nu_tau},{jac_big_mod});    

        funcs.mod_second = Function('Mod_second',{q,qd,qdd,mu,nu_q,nu_qd,nu_tau},{full_hessian_mod, jac_big_mod});
        
        A = - 1/2*full_hessian_mod(1:Nb,1:Nb) - jac_big_mod(1:Nb,:);
        H_qq_mod = (A+A');
        H_qd_qd_mod = (-full_hessian_mod(Nb+1:end, Nb+1:end));
        H_q_qd_mod = (-full_hessian_mod(1:Nb, Nb+1:end) - jac_big_mod(Nb+1:2*Nb,:)');
        H_q_tau = (-jac_big_mod(2*Nb+1:end,:)');
       
        funcs.mod_second_all = Function('Mod_second_all',{q,qd,qdd,mu,nu_q,nu_qd,nu_tau},{H_qq_mod,H_qd_qd_mod,H_q_qd_mod,H_q_tau, full_hessian_mod, jac_big_mod,out_big});
        
        FOP = [nu_q,nu_qd,nu_tau];
        ID_second_hess = [H_qq_mod H_q_qd_mod ; H_q_qd_mod' H_qd_qd_mod];
        second_tmp = Function('second_tmp',{q,qd,qdd,mu,FOP},{ID_second_hess,H_q_tau});        
          
        %
        % Reverse mode AD to get gradient of mu' * tauID
        ID_q = jtimes(tauID,q,mu,true);
        ID_qd = jtimes(tauID,qd,mu,true);

        % Relates second-order partials of InvDyn to FwdDyn
        H_q_q = -1/2 * jaco(ID_q, q) -  jtimes(M_mult(q,mu) , q, nu_q,true);
        H_q_q = H_q_q+H_q_q';

        H_qd_qd = -jaco(ID_qd, qd);
        H_q_qd  = -jaco(ID_q,qd) - jtimes(M_mult(q,mu) , q, nu_qd,true);
        H_q_tau = -jtimes(M_mult(q,mu),q, nu_tau, true);
        
        % Note: I don't currently have a "common" way to compute all these
        % partials together at the same time. That could help speed it up
        % if we could find the commonality.

        funcs.helper1 = Function('Helper1',{q,mu,nu_q},{jtimes(M_mult(q,mu) , q, nu_q,true)} );
        funcs.helper2 = Function('Helper2',{q,mu,nu_q},{jacobian(nu_q'*M_mult(q,mu),q) } );
        
        funcs.q_q = Function('H_ID_q_q', {q,qd,qdd,mu,nu_q}, {H_q_q} );
        funcs.qd_qd = Function('H_ID_qd_qd', {q,qd,mu}, {H_qd_qd} );
        funcs.q_qd = Function('H_ID_q_qd', {q,qd,mu,nu_qd}, {H_q_qd} );
        funcs.q_tau = Function('H_ID_q_tau', {q,qdd, mu,nu_tau}, {H_q_tau} );

        funcs.all_second = Function('H_ID_all', {q,qd,qdd,mu,nu_q,nu_qd,nu_tau}, {H_q_q,H_qd_qd,H_q_qd,H_q_tau} );
        %}
       
        % All at once
        
        tau    = MX.sym('tau',[model.NB 1]);
        lambda = MX.sym('lambda',[model.NB 1]);
        qdd    = FDab_casadi(model,q,qd,tau);
        mu     = InvM_Mult(q,lambda);
        ID_all_first = -InvM_Mult(q, Diff_ID_all(q,qd,qdd));        
        [ID_second_hess,H_q_tau]= second_tmp(q,qd,qdd,mu,ID_all_first);
        funcs.mod_all = Function('Mod_all',{q,qd,tau,lambda},{ID_all_first,ID_second_hess  ,H_q_tau});
    end
end