function [P,Pb,T,Tb] = MeshMatrix(pde,N1,N2,N3,deg)
% generate the Mesh information
% input: pde: pde model problem
%        N1,N2£¬N3: the number of partition on x,y,z-axis
%        deg
% output: P, Pb, T and Tb

if deg == 1
%% generate the Mesh information for one-order basis function
    P = zeros(3,(N1+1)*(N2+1)*(N3+1));
    for k = 1:N3+1
        for j = 1:N1+1
            for i = 1:N2+1
                P(1,i+(j-1)*(N2+1)+(k-1)*(N1+1)*(N2+1)) = pde.front + (j-1)*(pde.back-pde.front)/N1;
                P(2,i+(j-1)*(N2+1)+(k-1)*(N1+1)*(N2+1)) = pde.left + (i-1)*(pde.right-pde.left)/N2;
                P(3,i+(j-1)*(N2+1)+(k-1)*(N1+1)*(N2+1)) = (k-1)*(pde.top-pde.bottom)/N3;
            end
        end
    end
    
    T = zeros(8,N1*N2*N3);
    for k=1:N3
        for j = 1:N1
            for i = 1:N2
                T(1,i+(j-1)*N2+(k-1)*N1*N2) = i+(j-1)*(N2+1)+(k-1)*(N1+1)*(N2+1);
                T(2,i+(j-1)*N2+(k-1)*N1*N2) = i+(j)*(N2+1)+(k-1)*(N1+1)*(N2+1);
                T(3,i+(j-1)*N2+(k-1)*N1*N2) = i+(j)*(N2+1)+1+(k-1)*(N1+1)*(N2+1);
                T(4,i+(j-1)*N2+(k-1)*N1*N2) = i+1+(j-1)*(N2+1)+(k-1)*(N1+1)*(N2+1);
                T(5,i+(j-1)*N2+(k-1)*N1*N2) = i+(j-1)*(N2+1)+(k)*(N1+1)*(N2+1);
                T(6,i+(j-1)*N2+(k-1)*N1*N2) = i+(j)*(N2+1)+(k)*(N1+1)*(N2+1);
                T(7,i+(j-1)*N2+(k-1)*N1*N2) = i+(j)*(N2+1)+1+(k)*(N1+1)*(N2+1);
                T(8,i+(j-1)*N2+(k-1)*N1*N2) = i+1+(j-1)*(N2+1)+(k)*(N1+1)*(N2+1);
            end
        end
    end
    Pb = P;
    Tb = T;
    
elseif deg == 2
%% generate the Mesh information for two-order basis function
    P = zeros(2,(N1+1)*(N2+1));
    for j = 1:N1+1
        for i = 1:N2+1
            P(1,i+(j-1)*(N2+1)) = pde.left + (j-1)*(pde.right-pde.left)/N1;
            P(2,i+(j-1)*(N2+1)) = pde.bottom + (i-1)*(pde.top-pde.bottom)/N2;
        end
    end
    
    T = zeros(3,2*N1*N2);
    for j = 1:N1
        for i = 1:N2
            T(1,2*i-1+(j-1)*(2*N2)) = i+(j-1)*(N2+1);
            T(2,2*i-1+(j-1)*(2*N2)) = i+(j)*(N2+1);
            T(3,2*i-1+(j-1)*(2*N2)) = i+1+(j-1)*(N2+1);
            T(1,2*i+(j-1)*(2*N2)) = i+1+(j-1)*(N2+1);
            T(2,2*i+(j-1)*(2*N2)) = i+(j)*(N2+1);
            T(3,2*i+(j-1)*(2*N2)) = i+1+(j)*(N2+1);
        end
    end
    Pb = P;
    Tb = T;
    
end
end

