function [A] = StiffnessMatrix(P,Pb,T,Tb,coff_func,deg,...
    x_der_index_trial,y_der_index_trial,z_der_index_trial,x_der_index_test,y_der_index_test,z_der_index_test)
% generate the stiffness matrix
% input: P,Pb,T,Tb
% output: stiffness matrix A

nofbasis = size(Pb,2);% number of basis function
A = sparse(nofbasis,nofbasis);% initialize the stiffness matrix

nofelem = size(T,2);% number of element
for i=1:nofelem %loop with element
    vertices = P(:,T(:,i));
    [Gauss_weights,Gauss_nodes] = Generate_Gauss(vertices);
    for j = 1:size(Tb,1)%loop with trial basis function
        for k=1:size(Tb,1)%loop with test basis function
           A(Tb(k,i),Tb(j,i)) = A(Tb(k,i),Tb(j,i)) + sum(Gauss_weights.*feval(coff_func,Gauss_nodes(:,1),Gauss_nodes(:,2),Gauss_nodes(:,3)).*...
            BasisFunc(vertices,Gauss_nodes,deg,j,x_der_index_trial,y_der_index_trial,z_der_index_trial).*...
            BasisFunc(vertices,Gauss_nodes,deg,k,x_der_index_test,y_der_index_test,z_der_index_test));
        end
    end
end
end

