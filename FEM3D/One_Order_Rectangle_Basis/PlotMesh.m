function PlotMesh(N1,N2,N3,P,T)
% plot the finite element method mesh
% input: N1,N2,N3: the number of partition on x,y,z axis
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
    
    % element tag
    for i = 1:size(T,2)
       temp = num2str(i);
       text(sum(P(1,T(:,i)))/4,sum(P(2,T(:,i)))/4,temp,'Color','m');
    end
    
elseif element_type == 8 %uniform cubic mesh in 3D
    pdeback = min(P(1,:));
    pdefront = max(P(1,:));
    pdeleft = min(P(2,:));
    pderight = max(P(2,:));
    pdebottom = min(P(3,:));
    pdetop = max(P(3,:));
    
    for i = pdeback:(pdefront-pdeback)/N1:pdefront
        for j = pdeleft:(pderight-pdeleft)/N2:pderight
            plot3([i i], [j j], [pdebottom, pdetop]);
            hold on
        end
    end
    
    for i = pdeback:(pdefront-pdeback)/N1:pdefront
        for k = pdebottom:(pdetop-pdebottom)/N3:pdetop
            plot3( [i i], [pdeleft pderight],[k k]);
            hold on
        end
    end
    
    for j = pdeleft:(pderight-pdeleft)/N2:pderight
        for k = pdebottom:(pdetop-pdebottom)/N3:pdetop
            plot3([pdeback pdefront], [j j], [k k]);
            hold on
        end
    end
    
    % node tag
    for i = 1:size(P,2) % node
        temp = num2str(i);
        text(P(1,i)+0.01*1/N1,P(2,i)-0.01*1/N2,P(3,i)+0.01*1/N3,temp,'Color','b');
    end
    xlabel('x');
    ylabel('y');
    zlabel('z');
end

end
