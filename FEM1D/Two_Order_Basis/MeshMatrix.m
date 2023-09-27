function [P,Pb,T,Tb] = MeshMatrix(left,right,n,deg)
% generate the Mesh information
% input: left,right: the start and end point of the interval
%        n: the number of element
%        deg
% output: P, Pb, T and Tb

if deg == 1
%% generate the Mesh information for one-order basis function
    P = left:(right-left)/n:right;
    T = [1:n;2:n+1];
    
    Pb = P;
    Tb = T;
    
elseif deg == 2
%% generate the Mesh information for two-order basis function
    P = left:(right-left)/n:right;
    T = [1:n;2:n+1];
    
    Pb = left:(right-left)/n/2:right;
    Tb = [1:2:2*n-1;2:2:2*n;3:2:2*n+1];
end
end

