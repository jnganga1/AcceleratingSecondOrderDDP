function [dV, Vx, Vxx, du, K, success] = BackwardPass_RNEA(xbar, ubar, params,regularization)
    %global RNEAbackTime RNEAbackIters
    success = 1;
    
    % Initialization
    Vx = zeros(params.x_size, params.N);
    Vxx = zeros(params.x_size, params.x_size, params.N);
    
    du = zeros(params.u_size, params.N);
    
    K = zeros(params.u_size, params.x_size, params.N);
    
%     Quu_store =  zeros(params.u_size, params.u_size, params.N);
    
    % For CodeGen
    %K = zeros(1,2, params.N);
    
    xf = xbar(:,params.N);
    [m, n] =params.L.LFderivs(xf);
    Vx(:   ,params.N)= full(m);
    Vxx(:,:,params.N)=full(n);

    dV = 0;
    for i = (params.N-1):-1:1
        xi = xbar(:,i);
        ui = ubar(:,i);
        Vxi  = Vx (:  ,i+1);
        Vxxi = Vxx(:,:,i+1);
               
        
        % If we are doing full second order DDP, add in regularization
        if params.iLQR == 0
%             [Qx, Qu,Ham_fxx,Qxx, Quu, Qux,Fxx]=NewQinfO_DDP(xi,ui, Vxi, Vxxi,params.dt);
            tic
            [Qx, Qu,Qxx, Quu, Qux]= CasadiQinfo_RNEA(xi,ui, Vxi, Vxxi,params);
            %RNEAbackTime = RNEAbackTime + toc;
            %RNEAbackIters = RNEAbackIters + 1;
           
%             [Quu,Qxx]=Regularizer(Quu,Qxx);
            
            
            
            Qxx = Qxx + eye(params.x_size)*regularization;
            
            % FOR CODEGEN
            %Quu = Quu + eye(1)*regularization;
            Quu_Before = Quu;
            Quu = Quu + eye(params.u_size)*regularization;
            
            % Make sure Quu is PD, if not, exit and increase regularization
            [~, p] = chol(Quu-eye(params.u_size)*1e-9);
            if p ~= 0
                success = 0;
                break
            end
%             Quu_store(:,:,i) = Quu;
        else
            [Qx, Qu, Qxx, Quu, Qux] = CasadiQinfoILQR_RNEA(xi, ui, Vxi, Vxxi,params);
            warning('iLQR loop.')
        end
        
        % Standard equations
        du(:,i)  = -Quu\Qu;
        K(:,:,i) = -Quu\Qux;
        Vx(:,i)    = Qx  - (Qux')*(Quu\ Qu );
        Vxx(:,:,i) = Qxx - (Qux')*(Quu\ Qux);
        dV = dV - Qu'*(Quu\Qu);
    end
end