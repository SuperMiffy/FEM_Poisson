function [error_L2] = ErrorEstimateL2(P,Pb,T,Tb,x_FEM1D,deg)
% L2 error estimate
% input: Pb, Tb, x_FEM1D,deg

nofelem = size(T,2);%number of element

error_L2 = 0;
for i = 1:nofelem
   vertices = P(T(:,i)); 
   [Gauss_weights,Gauss_nodes] = Generate_Gauss(vertices);
   x_local = 0;
   for j = 1:size(Tb,1)
      x_local =  x_local + x_FEM1D(Tb(j,i))*BasisFunc(vertices,Gauss_nodes,deg,j,0);
   end
   error_L2 = error_L2 + sum(Gauss_weights.*(x_true(Gauss_nodes)-x_local).^2);
end
error_L2 = sqrt(error_L2);
end

