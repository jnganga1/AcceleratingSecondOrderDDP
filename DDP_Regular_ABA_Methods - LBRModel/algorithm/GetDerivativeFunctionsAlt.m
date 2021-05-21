function funcs = GetDerivativeFunctionsAlt(model, second_order)
    import casadi.*
    
    % Casadi symbolics
    q = MX.sym('q',[model.NB 1]);
    qd = MX.sym('qd',[model.NB 1]);
    qdd = MX.sym('qdd',[model.NB 1]);    
    tau = MX.sym('tau',[model.NB 1]);
    
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
    Diff_ID_q = Function('Diff_ID_q',{q,qd,qdd}, {jacobian(tauID, q)}); % O(N^2) function
    Diff_ID_qd = Function('Diff_ID_qd',{q,qd}, {jacobian(tauID, qd)}); % O(N^2) function

    % All together
    Diff_ID_all = Function('Diff_ID_all',{q,qd,qdd}, {[jacobian(tauID, [q;qd]) -MX.eye(model.NB)] }); % O(N^2) function

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
        
        nu = MX.sym('nu',[Nb,Nb]);
        nu_q = MX.sym('nu_q',[Nb,Nb]);
        nu_qd = MX.sym('nu_qd',[Nb,Nb]);
        nu_tau = MX.sym('nu_tau',[Nb,Nb]);

        % Reverse mode AD to get gradient of mu' * tauID
        ID_q = jtimes(tauID,q,mu,true);
        ID_qd = jtimes(tauID,qd,mu,true);

        % See notes :) Relates second-order partials of InvDyn to FwdDyn
        H_q_q = -1/2 * jacobian(ID_q, q) -  jtimes(M_mult(q,mu) , q, nu_q,true);
        H_q_q = H_q_q+H_q_q';

        H_qd_qd = -jacobian(ID_qd, qd);
        H_q_qd  = -jacobian(ID_q,qd) - jtimes(M_mult(q,mu) , q, nu_qd,true);
        H_q_tau = -jtimes(M_mult(q,mu),q, nu_tau, true);
        
        % Note: I don't currently have a "common" way to compute all these
        % partials together at the same time. That could help speed it up
        % if we could find the commonality.

        funcs.q_q = Function('H_ID_q_q', {q,qd,qdd,mu,nu_q}, {H_q_q} );
        funcs.qd_qd = Function('H_ID_qd_qd', {q,qd,mu}, {H_qd_qd} );
        funcs.q_qd = Function('H_ID_q_qd', {q,qd,mu,nu_qd}, {H_q_qd} );
        funcs.q_tau = Function('H_ID_q_tau', {q,qdd, mu,nu_tau}, {H_q_tau} );

        funcs.all_second = Function('H_ID_all', {q,qd,qdd,mu,nu_q,nu_qd,nu_tau}, {H_q_q,H_qd_qd,H_q_qd,H_q_tau} );
    end
end