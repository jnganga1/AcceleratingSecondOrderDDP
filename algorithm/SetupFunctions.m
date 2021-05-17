function params = SetupFunctions(x, u, l, lf, f, params)
    addpath([pwd '/support']);
    
    x_size = length(x);
    u_size = length(u);
    
    
    params.x_size = x_size;
    params.u_size = u_size;
    
    lambda = sym('lambda',[x_size 1]','real');
    
    H = l + lambda'*f;
    Hx = jacobian(H,x);
    Hu = jacobian(H,u);

    Vxx = sym('Vxx',[x_size x_size],'real');

    lfx  = jacobian(lf,x);
    lfxx = hessian(lf,x);

    if params.iLQR
        Hxx = hessian(l,x);
        Huu = hessian(l,u);
        Hux = jacobian(jacobian(l,u),x);
    else
        Hxx = jacobian(Hx, x);
        Huu = jacobian(Hu, u);
        Hux = jacobian(Hu, x);
    end
    
    fx = jacobian(f,x);
    fu = jacobian(f,u);
    
    Qxx = Hxx + fx'*Vxx*fx;
    Quu = Huu + fu'*Vxx*fu;
    Qux = Hux + fu'*Vxx*fx;

    matlabFunction(l,  'vars',{x, u},'file','support/running_cost');
    matlabFunction(lf, 'vars',{x}   ,'file','support/final_cost');
    matlabFunction(f,  'vars',{x u} ,'file','support/dynamics');

    matlabFunction(Hx', Hu', Qxx, Quu, Qux,'file','support/QInfo','vars',{x u lambda Vxx},'optimize',1==1);
    matlabFunction(lfx, lfxx ,'vars',{x},'file','support/lfInfo','optimize',1==1);
end
