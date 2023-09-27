function [boundary_matrix] = BoundaryMatrix(Pb)
%generate the boundary matrix
% input: Pb
% output: boundary_matrix

boundary_matrix = zeros(3,2);% initialize

boundary_matrix(1,1) = -1;
boundary_matrix(1,2) = 0;
boundary_matrix(2,1) = 1;
boundary_matrix(2,2) = size(Pb,2);
boundary_matrix(3,1) = 0;
boundary_matrix(3,2) = -pi;
end

