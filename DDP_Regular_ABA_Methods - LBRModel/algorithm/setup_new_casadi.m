function [allFuncs,params] = setup_new_casadi(x, u, l,Hmat,bVec,params)
   import casadi.*

   x_size = length(x);
   u_size = length(u);
        
   params.x_size = x_size;
   params.u_size = u_size;
   
   allFuncs.Hmat = Function('Hmat',{x,u},{Hmat})
   allFuncs.Hmat_x = Function('Hmat',{x,u},{jacobian(Hmat,x)})
   here =  jacobian(Hmat,x)
   here2 = jacobian(here,x)
   allFuncs.Hmat_xx = Function('Hmat',{x,u},{here2})
    
    %lf Partials (done outside of this script)
%     lfx  = jacobian(lf,x);
%     lfxx = hessian(lf,x);
    %l partials
    lx = simplify(jacobian(l,x)); 
    lu = simplify(jacobian(l,u));
%     lxx = hessian(l,x);
    lxx = simplify(jacobian(jacobian(l,x),x));
    lux = simplify(jacobian(jacobian(l,u),x));
    luu = simplify(jacobian(jacobian(l,u),u));
    
    %Hmat partials 
    [Hmat_x,Hmat_xx] = Mpartials(Hmat,x);  
    

    
    %bvec partials
    bVec_x = simplify((jacobian(bVec,x))); dimX = size(bVec_x);
    bVec_u = simplify((jacobian(bVec,u))); dimU = size(bVec_u);
    
    bVec_xx = sym(zeros(length(bVec),x_size,x_size));
    bVec_uu = sym(zeros(length(bVec),u_size,u_size));
    bVec_ux = sym(zeros(length(bVec),u_size,x_size));
    for i =1:dimX(1)
        for j=1:dimX(2)
           bVec_xx(i,j,:)=simplify(jacobian(bVec_x(i,j),x));            
        end
    end
    for i =1:dimU(1)
        for j=1:dimU(2)
           bVec_uu(i,j,:)=simplify(jacobian(bVec_u(i,j),u));            
        end
    end
    
     for i =1:dimU(1)
        for j=1:dimU(2)
           bVec_ux(i,j,:)=simplify(jacobian(bVec_u(i,j),x));            
        end
     end  
    
     
    funcs.l = Function('Lfunc',{x, u},{l})
%     funcs.bVec_part
    
    matlabFunction(l,  'vars',{x, u},'file',[fPath,'/LCost'],'optimize',1==1);
%     matlabFunction(lf, 'vars',{x}   ,'file','Partials/Lfcost');
%     matlabFunction(lfx, lfxx ,'vars',{x},'file','Partials/LfPartials','optimize',1==1);
    matlabFunction(lx,lxx,lu,luu,lux,'vars',{x,u},'file',[fPath,'/LPartials'],'optimize',1==1);
    
    matlabFunction(Hmat, 'vars',{x},'file',[fPath,'/HmatCalc'],'optimize',1==1);
    matlabFunction(Hmat_x,Hmat_xx,'vars',{x},'file',[fPath,'/HmatCalcDeriv'],'optimize',1==1);
    
    matlabFunction(bVec, 'vars',{x,u},'file',[fPath,'/bVecCalc'],'optimize',1==1);
    matlabFunction(bVec_x,bVec_u,bVec_xx,bVec_uu,bVec_ux,...
        'vars',{x,u},'file',[fPath,'/bVecPartialsCalc'],'optimize',1==1);
    b = [fPath,'/params.mat'];
    save(b,'params');
    1==1; %EndOfLine
end