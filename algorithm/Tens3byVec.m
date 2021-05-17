function [Out] = Tens3byVec(Tensor,Vec,strArg)
    dim = size ( Tensor );
%     Out= sym(zeros ( dim (1) ,dim (2) ));
%     Out = zeros(dim(1),dim(2));
    for i=1:dim(3)      
        switch strArg 
            case 'post'
                Out(:,i) = Tensor(: ,: ,i)*Vec;
            case 'pre'
                Out(:,i) = Vec * Tensor(: ,: ,i);
        end     
%         Out = Out + Tensor(:,:,i)*Vec(i);
    end
end 