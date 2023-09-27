function [b] = LoadVector(P,Pb,T,Tb,deg,x_der_index,y_der_index)
% generate the load vector
% input: P,Pb,T,Tb
% output: b

nofbasis = size(Pb,2);% number of basis function
b = zeros(nofbasis,1);%initial

nofelem = size(T,2);%number of element
for i=1:nofelem
    vertices = P(:,T(:,i));
    [Gauss_weights,Gauss_nodes] = Generate_Gauss(vertices);
   for k=1:size(Tb,1)
      b(Tb(k,i))=b(Tb(k,i))+sum(Gauss_weights.*RHF(Gauss_nodes(:,1),Gauss_nodes(:,2)).*...
                 BasisFunc(vertices,Gauss_nodes,deg,k,x_der_index,y_der_index)) ;
   end
end
end

