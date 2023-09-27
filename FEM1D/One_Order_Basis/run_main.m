function [error_L2,error_H1] = run_main(n)
% main function of the finite element method

left = 0;
right = 1;

% n = 10;% number of elements

%% Generate the Mesh information
[P,Pb,T,Tb] = MeshMatrix(left,right,n);

%% Generate the stiffness matrix
A = StiffnessMatrix(P,Pb,T,Tb,'p',1,1)+...
    StiffnessMatrix(P,Pb,T,Tb,'q',1,0)+...
    StiffnessMatrix(P,Pb,T,Tb,'r',0,0);

%% Generate the load vector
b = LoadVector(P,Pb,T,Tb,0);

%% Treat the boundary condition
boundary_matrix = BoundaryMatrix(Pb);
[A,b] = TreatBoundaryCondition(A,b,boundary_matrix);

%% solve the linear equation systems
x_FEM1D = A\b;

%% Error estimate
error_L2 = ErrorEstimateL2(P,Pb,T,Tb,x_FEM1D);
error_H1 = ErrorEstimateH1(P,Pb,T,Tb,x_FEM1D);

%% Error order
end