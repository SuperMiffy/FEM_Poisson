function [boundary_matrix] = BoundaryMatrix(N1,N2)
% generate the boundary matrix
% note: boundary type: -1-Dirichlet; 0-Neumann; 1-Robin
% input: N1,N2
% output: boundary_matrix

boundary_matrix = zeros(3,2*(N1+N2));% initialize

boundary_matrix(1,:) = -1;% if Neumann or Robin, attention

%bottom node
for i = 1:N1
    boundary_matrix(2,i) = 1+(i-1)*(N2+1);
end
%right node
for i = N1+1:N1+N2
    boundary_matrix(2,i) = N1*(N2+1)+i-N1;
end

%top node
for i = N1+N2+1:N1+N2+N1
    boundary_matrix(2,i) = N2+1+(N1+N2+N1+1-i)*(N2+1);
end

%left node
for i = N1+N2+N1+1:N1+N2+N1+N2
    boundary_matrix(2,i) = N2+1-(i-(N1+N2+N1+1));

end

%Dirichlet boundary condition value
% for i = 1:2*(N1+N2)
%    boundary_matrix(3,i) = x_true(); 
% end


end

