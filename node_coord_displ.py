def node_coord_displ(h,l,m,n,load_opt):

    import math
    coord = []
    displ = []
    elnode = []
    node = 0
    deltx = l / (2. * n)
    for i in range(2 * n + 1):
        multiplier = ((i + 1) % 2) + 1
        delty = h / (float(multiplier) * m)
        for j in range(multiplier * m + 1):
            node = node + 1
            x = i * deltx
            y = j * delty - h / 2.
            coord.append([node, x, y])
            if ((load_opt == 2) or (load_opt == 3) or (load_opt == 4)):
                u = (R - y) * math.sin(x / R) - x
                v = R - (R - y) * math.cos(x / R) - y
                if load_opt == 2:
                    displ.append([node, u, v])
                    coord.append([node, x, y])
                else:
                    coord.append([node(x + u), (y + v)])

    el = 0
    for i in range(n):
        for j in range(m):
            elnode1 = (3 * m + 2) * (i) + 2 * (j) + 1
            elnode2 = (3 * m + 2) * (i) + 2 * (j) + 3
            elnode3 = (3 * m + 2) * (i + 1) + 2 * (j) + 1
            elnode4 = (3 * m + 2) * (i + 1) + 2 * (j) + 3
            elnode5 = 2 * m + 1 + (3 * m + 2) * (i) + (j + 1)
            elnode6 = (3 * m + 2) * (i + 1) + 2 * (j) + 2
            elnode7 = 2 * m + 1 + (3 * m + 2) * (i) + (j + 1) + 1
            elnode8 = (3 * m + 2) * (i) + 2 * (j) + 2
            el1 = el + 1
            elnode.append([el1, elnode1, elnode3, elnode4, elnode2, elnode5, elnode6, elnode7, elnode8])
            el = el + 1

    return coord, displ, elnode, node, el
