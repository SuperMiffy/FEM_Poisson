function [error_L2,error_H1] = run_main(n)
% main function of the finite element method

left = 0;
right = 1;

deg = 2;% order of basis function
% n = 10;% number of elements

%% Generate the Mesh information
[P,Pb,T,Tb] = MeshMatrix(left,right,n,deg);

%% Generate the stiffness matrix
A = StiffnessMatrix(P,Pb,T,Tb,'p',deg,1,1)+...
    StiffnessMatrix(P,Pb,T,Tb,'q',deg,1,0)+...
    StiffnessMatrix(P,Pb,T,Tb,'r',deg,0,0);

%% Generate the load vector
b = LoadVector(P,Pb,T,Tb,deg,0);

%% Treat the boundary condition
boundary_matrix = BoundaryMatrix(Pb);
[A,b] = TreatBoundaryCondition(A,b,boundary_matrix);

%% solve the linear equation systems
x_FEM1D = A\b;

%% Error estimate
error_L2 = ErrorEstimateL2(P,Pb,T,Tb,x_FEM1D,deg);
error_H1 = ErrorEstimateH1(P,Pb,T,Tb,x_FEM1D,deg);

%% Error order
end
