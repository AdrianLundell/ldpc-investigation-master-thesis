import unittest
import qc_ldpc.qc_graph_tools as qc
import numpy as np 

class PathSearchTest(unittest.TestCase):

    def test_no_edges(self):
        """
        o

        ■
        """
        G = qc.QC_tanner_graph(1,1,1)
        
        self.assertDictEqual(qc.shortest_distances(G, 0), {})
        self.assertDictEqual(qc.shortest_cycles(G, 0), {})


        """   
        o o o o o

        ■ ■ ■ ■ ■
        """
        G = qc.QC_tanner_graph(1,1,5)
        
        self.assertDictEqual(qc.shortest_distances(G, 0), {})
        self.assertDictEqual(qc.shortest_cycles(G, 0), {})
        
    def test_no_cycles(self):
        """
        o
        |
        ■
        """
        G = qc.QC_tanner_graph(1,1,1)
        G.add_edges([(0,0)])

        self.assertDictEqual(qc.shortest_distances(G, 0), {1 : 1})
        self.assertDictEqual(qc.shortest_cycles(G, 0), {})


        """
        o o o o o
        | | | | |
        ■ ■ ■ ■ ■
        """
        G = qc.QC_tanner_graph(1,1,5)
        G.add_cyclical_edge_set(0,0)

        self.assertDictEqual(qc.shortest_distances(G, 0), {5 : 1})
        self.assertDictEqual(qc.shortest_cycles(G, 0), {})

        """
        o o
        |/|
        ■ ■ 
        """
        G = qc.QC_tanner_graph(2,2,1)
        G.add_edges([(0,0), (0,1), (1,1)])

        self.assertDictEqual(qc.shortest_distances(G, 0), {1:2, 2:1, 3:1})
        self.assertDictEqual(qc.shortest_cycles(G, 0), {})

    def test_small_cycles(self):
        """
        o o o o
        |x| |x|
        ■ ■ ■ ■ 
        """
        G = qc.QC_tanner_graph(2,2,2)
        G.add_edges([(0,0), (0,1), (1,0), (1,1)])
        G.add_edges([(2,2), (2,3), (3,2), (3,3)])

        #First cycle
        self.assertDictEqual(qc.shortest_distances(G, 0), {1 : 2, 
                                                       4 : 1, 
                                                       5 : 1})
        self.assertDictEqual(qc.shortest_cycles(G, 0), {4 : 4, 
                                                        5 : 4}) 

        #Second cycle    
        self.assertDictEqual(qc.shortest_distances(G, 7), {6 : 2, 
                                                       2 : 1, 
                                                       3 : 1})
        self.assertDictEqual(qc.shortest_cycles(G, 7), {2 : 4, 
                                                        3 : 4}) 

    def test_bottle_neck(self):
        """
        o o o
         /|x|
        ■ ■ ■
        """
        G = qc.QC_tanner_graph(3,3,1)
        G.add_edges([(0,1), (1,1), (1,2), (2,1), (2,2)])

        self.assertDictEqual(qc.shortest_distances(G, 0), {1 : 2, 
                                                       2 : 2, 
                                                       4 : 1,
                                                       5 : 3})
        self.assertDictEqual(qc.shortest_cycles(G, 0), {}) 

        """
        o o o
         /|x|
        ■ ■ ■ (+ edge (0,2))
        """
        G = qc.QC_tanner_graph(3,3,1)
        G.add_edges([(0,1), (1,1), (1,2), (2,1), (2,2), (0,2)])

        self.assertDictEqual(qc.shortest_distances(G, 0), {1 : 2, 
                                                       2 : 2, 
                                                       4 : 1,
                                                       5 : 1})
        self.assertDictEqual(qc.shortest_cycles(G, 0), {4 : 4, 5 : 4})

    def test_fully_connected(self):
        """
        All cns and vns connnected
        """
        G = qc.QC_tanner_graph(4,4,1)
        G.add_edges([(0,0),(0,1),(0,2), (0,3)])
        G.add_edges([(1,0),(1,1),(1,2), (1,3)])
        G.add_edges([(2,0),(2,1),(2,2), (2,3)])
        G.add_edges([(3,0),(3,1),(3,2), (3,3)])

        self.assertDictEqual(qc.shortest_distances(G, 0), {1 : 2, 2 : 2, 3 : 2,
                                                           4 : 1, 5 : 1, 6  : 1, 7 : 1})
        self.assertDictEqual(qc.shortest_cycles(G, 0), {4 : 4, 5 : 4, 6 : 4, 7 : 4})

    def test_long_cycle(self):
        """
        o o o o 
        |/|/|/|
        ■ ■ ■ ■ (+ edge (3,0))
        """
        G = qc.QC_tanner_graph(4,4,1)
        G.add_edges([(0,0),(0,1),
                     (1,1),(1,2),
                     (2,2),(2,3),
                     (3,3),(3,0)])
                     
        self.assertDictEqual(qc.shortest_distances(G, 0), {1: 2, 
                                                           2: 4, 
                                                           3: 2, 
                                                           4:1, 
                                                           5:1, 
                                                           6:3, 
                                                           7: 3})

        self.assertDictEqual(qc.shortest_cycles(G, 0), {4 : 8, 5 : 8})