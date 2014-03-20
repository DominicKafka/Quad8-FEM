function [StressOut, StressNode] = nodal_stresses(elnodes, stress)

StressOut = [elnodes(:,1), stress(:, [1:6, 10:12, 7:9])];

factors = [(1+sqrt(3)/2), -0.5, -0.5, (1-sqrt(3)/2)]';

cols = [[2, 5, 11, 8]
        [5, 2, 8, 11]
        [8, 5, 11, 2]
        [11, 2, 8, 5]];

StressNode = elnodes(:,1);

for i = 0:2
    for row = cols'
        StressNode(:, row(1) + i) = StressOut(:, row + i)*factors;
    end
end