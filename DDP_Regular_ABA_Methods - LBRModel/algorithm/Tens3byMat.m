function [Out] = Tens3byMat(Tensor,Mat,strArg)
    %strArg
    dim = size ( Tensor );
%     Out=sym(zeros(dim));
    for i=1:dim(3)
           switch strArg 
               case 'post'
                Out(: ,: ,i) = Tensor(: ,: ,i)*Mat;
               case 'pre'
               Out(: ,: ,i) = Mat* Tensor(: ,: ,i);
           end 
    end
end 