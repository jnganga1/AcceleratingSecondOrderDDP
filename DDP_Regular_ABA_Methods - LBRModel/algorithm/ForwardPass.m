function [V,x,u] = ForwardPass(xbar, ubar, du, K, eps,params)
    V = 0;
    x = 0*xbar;
    u = 0*ubar;
    x(:,1) = xbar(:,1);
    xi = xbar(:,1);
    for i = 1: (params.N-1)
        dxi = xi - xbar(:,i);
        
        % Update with stepsize and feedback
        ui = ubar(:,i) + eps*du(:,i) + K(:,:,i)*dxi;
        u(:,i) = ui;
        if  any(isnan(ui))
            1==1;
        end
        
        % Add up cost
        V = V + full(params.L.L(xi, ui));
        
        % Propagate dynamics
        dt = params.dt;
        xi = CasadiDyn(params.F.qddABA,xi,ui, dt);
        x(:,i+1) = xi;
    end
    V = V + full(params.L.LF(x(:,params.N)));
   
end