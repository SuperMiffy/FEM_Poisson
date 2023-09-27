function [A,b] = TreatBoundaryCondition(A,b,boundary_matrix)
% treat the boundary condition
% input: A,b,boundary_matrix
% output: A,b

for i = 1:size(boundary_matrix,2)
   if boundary_matrix(1,i) == -1
       A(boundary_matrix(2,i),:) = 0;
       A(boundary_matrix(2,i),boundary_matrix(2,i)) = 1;
       
       b(boundary_matrix(2,i)) = boundary_matrix(3,i);
   elseif boundary_matrix(1,i) == 0
       b(boundary_matrix(2,i)) = b(boundary_matrix(2,i)) + boundary_matrix(3,i);
   elseif boundary_matrix(1,i) == 1
       % treat specially
   end
end

end

