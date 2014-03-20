function VonMises = calc_von_mises(StressNode, pois, plane)

VonMises = zeros(size(StressNode, 1), 4);
if (plane==1)
    for i=1:4
        VonMises(:,i) = StressNode(:,2+(i-1)*3).^2 - ...
                StressNode(:,2+(i-1)*3).*StressNode(:,3+(i-1)*3) + ...
                StressNode(:,3+(i-1)*3).^2 + 3*StressNode(:,4+(i-1)*3).^2;
        VonMises(:,i) = VonMises(:,i).^0.5;
    end
else
    for i=1:4
        VonMises(:,i) = (1-pois+pois^2)*(StressNode(:,2+(i-1)*3).^2 + ...
                        StressNode(:,3+(i-1)*3).^2) - (1+pois-pois^2)* ...
                        StressNode(:,2+(i-1)*3).*StressNode(:,3+(i-1)*3)+...
                        3*StressNode(:,4+(i-1)*3).^2;
        VonMises(:,i) = VonMises(:,i).^0.5;
    end
end    
