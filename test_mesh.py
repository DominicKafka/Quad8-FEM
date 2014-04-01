import unittest

from mesh_utils import node_coord_displ

class MeshTestCase(unittest.TestCase):
    def testNodeCoordDispl(self):
        h = 1.
        l = 10.
        m = 10
        n = 100
        option = 1.
        E = 200000.
        nu = 0.3
        t = 1.
        load_opt = 0
        Mag = 0
        R = 0
        RadD = 0
        V0 = 100.

        coord, displ, elnode, node = node_coord_displ(h, l, m, n, load_opt, R)

        self.assertEqual(node, 3221, "Incorrect number of nodes generated")
        self.assertEqual(elnode.shape[0], node, "Elnode is the wrong shape")

if __name__ == '__main__':
    unittest.main()