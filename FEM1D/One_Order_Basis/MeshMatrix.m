function [P,Pb,T,Tb] = MeshMatrix(left,right,n)
% generate the Mesh information
% input: left,right: the start and end point of the interval
%        n: the number of element
% output: P, Pb, T and Tb

%% generate the Mesh information for one-order basis function
P = left:(right-left)/n:right;
T = [1:n;2:n+1];

Pb = P;
Tb = T;

end

