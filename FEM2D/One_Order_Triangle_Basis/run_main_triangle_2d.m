function [error_L2,error_H1] = run_main_triangle_2d(N1,N2)
% main function of the 2D-finite element method
% the mesh is uniform triangle mesh
% N1:number of partition on X-axis; N2:number of partition on Y-axis

pde.left = 0;
pde.right = 1;
pde.bottom = 0;
pde.top = 1;

deg = 1;% order of basis function
% n = 10;% number of elements

%% Generate the Mesh information
[P,Pb,T,Tb] = MeshMatrix(pde,N1,N2,deg);

%% Generate the stiffness matrix
A = StiffnessMatrix(P,Pb,T,Tb,'p',deg,1,0,1,0)+...
    StiffnessMatrix(P,Pb,T,Tb,'p',deg,0,1,0,1)+...
    StiffnessMatrix(P,Pb,T,Tb,'q',deg,1,0,0,0)+...
    StiffnessMatrix(P,Pb,T,Tb,'q',deg,0,1,0,0)+...
    StiffnessMatrix(P,Pb,T,Tb,'r',deg,0,0,0,0);

%% Generate the load vector
b = LoadVector(P,Pb,T,Tb,deg,0,0);

%% Treat the boundary condition
boundary_matrix = BoundaryMatrix(N1,N2);%Dirichlet boundary condition
[A,b] = TreatBoundaryCondition(A,b,boundary_matrix);

%% solve the linear equation systems
x_FEM1D = A\b;

%% Error estimate
error_L2 = ErrorEstimateL2(P,Pb,T,Tb,x_FEM1D,deg);
error_H1 = ErrorEstimateH1(P,Pb,T,Tb,x_FEM1D,deg);

%% Error order
end