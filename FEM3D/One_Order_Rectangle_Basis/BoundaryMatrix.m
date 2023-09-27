function [boundary_matrix] = BoundaryMatrix(N1,N2,N3)
% generate the boundary matrix
% note: boundary type: -1-Dirichlet; 0-Neumann; 1-Robin
% input: N1,N2
% output: boundary_matrix

boundary_matrix = zeros(3,2*(N1+N2)*(N3+1)+2*(N1-1)*(N2-1));% initialize

boundary_matrix(1,:) = -1;% if Neumann or Robin, attention

% bounary edge nodes on plane z = pde.bottom
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

% loop by z = h_z*k
for k = 2:N3+1 % k=1 have generated
   boundary_matrix(2,2*(N1+N2)*(k-1)+1:2*(N1+N2)*k) = boundary_matrix(2,1:2*(N1+N2))+(N1+1)*(N2+1)*(k-1);
end

% generate the interal node on z=pde.bottom && z=pde.top
boundary_internal_node = zeros(1,2*(N1-1)*(N2-1));
for i = 2:N1
    for j = 2:N2
        boundary_internal_node(1,(i-2)*(N2-1)+j-1) = (i-1)*(N2+1)+j;
    end
end

boundary_internal_node((N1-1)*(N2-1)+1:2*(N1-1)*(N2-1)) = boundary_internal_node(1:(N1-1)*(N2-1)) + N3*(N1+1)*(N2+1);

% fill into the boundary matrix
boundary_matrix(2,2*(N1+N2)*(N3+1)+1:2*(N1+N2)*(N3+1)+2*(N1-1)*(N2-1)) = boundary_internal_node;

%Dirichlet boundary condition value
% for i = 1:2*(N1+N2)
%    boundary_matrix(3,i) = x_true(); 
% end


end

