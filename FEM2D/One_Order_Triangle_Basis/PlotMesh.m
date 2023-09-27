function PlotMesh(N1,N2,P,T)
% plot the finite element method mesh
% input: N1,N2: the number of partition on x,y axis
%        P-dim*nof_nodes ！！ node information matrix
%        T-element_type*nof_elemments ！！ element information matrix
% output: mesh with P and T

element_type = size(T,1);
if element_type == 3 %uniform triangle mesh partition
    x = reshape(P(1,:),N2+1,N1+1);
    y = reshape(P(2,:),N2+1,N1+1);
    
    trimesh(T',x,y);
    hold on
    
    % note tag
    for i = 1:size(P,2) % node
       temp = num2str(i);
       text(P(1,i)+0.01*1/N1,P(2,i)-0.01*1/N2,temp,'Color','blue');
    end
    
    for i = 1:size(T,2)
       temp = num2str(i);
       text(sum(P(1,T(:,i)))/3,sum(P(2,T(:,i)))/3,temp,'Color','m');
    end
elseif element_type == 4 %uniform rectangle mesh partition
    x = reshape(P(1,:),N2+1,N1+1);
    y = reshape(P(2,:),N2+1,N1+1);
    
    plot(x,y,x',y');
    
    hold on
    
    % note tag
    for i = 1:size(P,2) % node
       temp = num2str(i);
       text(P(1,i)+0.01*1/N1,P(2,i)-0.01*1/N2,temp,'Color','blue');
    end
    
    for i = 1:size(T,2)
       temp = num2str(i);
       text(sum(P(1,T(:,i)))/4,sum(P(2,T(:,i)))/4,temp,'Color','m');
    end
end

end
