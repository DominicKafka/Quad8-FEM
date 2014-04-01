function [coord,displ,elnode,node] = node_coord_displ(h,l,m,n,load_opt)
%FIXME: displ is only calculated if load_opt == 2, otherwise it is never defined.

node = 0;
deltx = l/(2*n);
for i = 1:2*n+1;
    multiplier = mod(i, 2) + 1;
    delty = h/(multiplier*m);
    for j = 1:(multiplier*m + 1);
        node = node + 1;
        x = (i - 1)*deltx;
        y = (j - 1)*delty - h/2;
        coord(node, :) = [node x y];
        if any(load_opt == [2, 3, 4])
            u = (R - y)*sin(x/R) - x;
            v = R - (R - y)*cos(x/R) - y;
            if load_opt == 2
                displ(node, :) = [node u v];
                coord(node, :) = [node x y];
            else
                coord(node, :) = [node (x+u) (y+v)];
            end
        end
    end
end

el = 0;
for i=1:n
    for j=1:m
        elnode1 = (3*m+2)*(i-1) + 2*(j-1) + 1;
        elnode2 = (3*m+2)*(i-1) + 2*(j-1) + 3;
        elnode3 = (3*m+2)*i + 2*(j-1) + 1;
        elnode4 = (3*m+2)*i + 2*(j-1) + 3;
        elnode5 = 2*m+1+(3*m+2)*(i-1) + j;
        elnode6 = (3*m+2)*i + 2*(j-1) + 2;
        elnode7 = 2*m+1+(3*m+2)*(i-1) + j+1;
        elnode8 = (3*m+2)*(i-1) + 2*(j-1) + 2;
        el1 = el+1;
        elnode(el1,:) = [el1 elnode1 elnode3 elnode4 elnode2 ...
            elnode5 elnode6 elnode7 elnode8];
        el=el+1;
    end
end
