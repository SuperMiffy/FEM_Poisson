

vertices = [-1,1,1,-1;-1,-1,1,1];
% vertices = [0,1,1,0;0,0,1,1];

[Gauss_weights,Gauss_nodes] = Generate_Gauss(vertices);

I = sum(Gauss_weights.*feval('p',Gauss_nodes(:,1),Gauss_nodes(:,2)));