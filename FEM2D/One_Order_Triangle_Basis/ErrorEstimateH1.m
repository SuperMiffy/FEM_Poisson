function [error_H1] = ErrorEstimateH1(P,Pb,T,Tb,x_FEM1D,deg)
% H1 error estimate
% input: Pb, Tb, x_FEM1D, deg

nofelem = size(T,2);%number of element

error_H1 = 0;
for i = 1:nofelem
   vertices = P(:,T(:,i)); 
   [Gauss_weights,Gauss_nodes] = Generate_Gauss(vertices);
   x_local_x_der = 0;
   x_local_y_der = 0;
   for j = 1:size(Tb,1)
      x_local_x_der =  x_local_x_der + x_FEM1D(Tb(j,i))*BasisFunc(vertices,Gauss_nodes,deg,j,1,0);
      x_local_y_der =  x_local_y_der + x_FEM1D(Tb(j,i))*BasisFunc(vertices,Gauss_nodes,deg,j,0,1);
   end
   error_H1 = error_H1 + sum(Gauss_weights.*(x_true_x_der(Gauss_nodes(:,1),Gauss_nodes(:,2))-x_local_x_der).^2)+...
       sum(Gauss_weights.*(x_true_y_der(Gauss_nodes(:,1),Gauss_nodes(:,2))-x_local_y_der).^2);
end
error_H1 = sqrt(error_H1);
end


