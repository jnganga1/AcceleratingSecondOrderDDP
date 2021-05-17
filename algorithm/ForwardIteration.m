function [V, x, u,rho,ThreeNumbers] = ForwardIteration(xbar, ubar, Vprev, du, K, dV, params)
    V = 0;
    eps = 1;
    rho=0; 
    ThreeNumbers.FwdFailTrkCnt=0;
    ThreeNumbers.FwdSuc = 0;
    ThreeNumbers.FwdNonSuc =0;
    ThreeNumbers.Flag = 0;
    
    x = 0*xbar;
    u = 0*ubar;
    
    while eps > 1e-20
        
        % Try to take the step
        FwdStrTime = tic;
        [V,x,u] = ForwardPass(xbar, ubar, du, K, eps, params);
        FwdEndTime = toc(FwdStrTime);
        if params.Debug    
%             fprintf('\t eps=%.3e \t DV=%.3e \t min=%.3e\n',eps, V-Vprev, params.gamma* eps*(1-eps/2)*dV );
        end
        
        % Check if V is small enough to accept the step
        if V < Vprev + params.gamma* eps*(1-eps/2)*dV
            % if so, exit
            change = V - Vprev; 
            rho =change / (eps*(1-eps/2)*dV);  
            ThreeNumbers.FwdSuc = FwdEndTime;            
            ThreeNumbers.Flag = 1; %success
            break
        end
        
        % Else, backtrack
        eps = params.beta * eps;
        ThreeNumbers.FwdFailTrkCnt = ThreeNumbers.FwdFailTrkCnt + 1;
        ThreeNumbers.FwdNonSuc = ThreeNumbers.FwdNonSuc + FwdEndTime;
    end
end
