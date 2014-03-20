function Tresca = calc_tresca(StressNode, pois, plane)

nelem = size(StressNode, 1);

Tresca = zeros(nelem, 4);

if plane == 1
    for j=1:nelem
        for i=1:4
            s = [StressNode(j,2+(i-1)*3) StressNode(j,4+(i-1)*3) 0
                 StressNode(j,4+(i-1)*3) StressNode(j,3+(i-1)*3) 0
                0 0 0];
            principal = eig(s);
            Tresca(j,i) = max(principal)-min(principal);
        end
    end
else
    for j=1:nelem
        for i=1:4
            s = [StressNode(j,2+(i-1)*3) StressNode(j,4+(i-1)*3) 0
                 StressNode(j,4+(i-1)*3) StressNode(j,3+(i-1)*3) 0
                0 0 pois*(StressNode(j,2+(i-1)*3)+StressNode(j,3+(i-1)*3))];
            principal = eig(s);
            Tresca(j,i) = max(principal)-min(principal);
        end
    end
end
