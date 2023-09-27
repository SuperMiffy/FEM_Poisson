function [Gauss_weights,Gauss_nodes] = Generate_Gauss(vertices)
% Gauss integral with five points
% input: vertices
% output: Gauss_weights and Gauss_nodes
Gauss_weights = [0.236926885056189, 0.478628670499366, 0.568888888888889,...
                 0.478628670499366, 0.236926885056189];
Gauss_nodes = [-0.906179845938664, -0.538469310105683, 0, ...
                0.538469310105683, 0.906179845938664];
            
Gauss_nodes = (vertices(2)-vertices(1))/2.*Gauss_nodes + (vertices(1)+vertices(2))/2;
Gauss_weights = Gauss_weights.*(vertices(2)-vertices(1))./2;
end

