from cmath import inf
import numpy as np 

class Message:
    
    bits : list

    def __init__(self, bits) -> None:
        """Initilises a message from an iterable of ones and zeros"""

        self.bits = list(bits)

    def validate(self, recovered_bits):
        """Returns the rate of correct bits for a recovered message"""

        assert len(self.bits) == len(recovered_bits), "Messages must be of equal length"
        return sum([bit1 == bit2 for bit1, bit2 in zip(self.bits, recovered_bits)]) / len(self.bits)

    def syndrome(self):
        pass

class Detector:
    
    variance = 0.5

    def read_bit(self, bit):
        """
        Returns the simulated log likelyhood ratio 
            llr(b) = log(P(b = 0)/P(b = 1)
        of bit b
        """
        pass 

    def read_llg(self, message):   
        """Returns the complete llg list for a message"""     
        return [self.read_bit(bit) for bit in message]

    def read_hard(self, message):
        """Returns the hard message"""
        llg_list = self.read_llg(message)
        message = [int (llg > 0) for llg in llg_list]
        return message 


class TannerGraph:

    message : Message
    bit_nodes : list 
    check_nodes : list

    correct_bit_ratio : float 
    n_iterations : float 

    max_iterations = 5 

    def __init__(self, parity_check_matrix) -> None:
        """Initialise a tanners graph from a parity_check_matrix H"""
        m, n = parity_check_matrix.size
        self.bit_nodes = [self.BitNode() for row in range(m)]
        self.check_nodes = [self.BitNode() for row in range(n)]

        for bit_node_index in range(m):
            self.bit_nodes[bit_node_index].init()

        for check_node_index in range(n):
            pass
            
    def init_message(self, message, detector):
        pass
    
    def upward_pass(self):
        for bit_node in self.bit_nodes:
            bit_node.upward_pass

    def downward_pass(self):
        for check_node in self.check_nodes:
            check_node.downward_pass

    def read_hard(self):
        pass

    def validate(self):
        pass

    class BitNode:

        init_value : float
        recieved_values_dict : dict         

        def upward_pass(self):
            """Perform one upward pass of the standard LDPC demo"""
            base_value = self.init_value + sum(self.recieved_values_dict.values())
            
            for check_node, recieved_value in self.recieved_values_dict.items():
                check_node.recieved_values_dict[self] = base_value - recieved_value

    class CheckNode:
        
        recieved_values_dict : dict 

        def downward_pass(self):
            """Perform one downward pass of the standard LDPC demo"""

            least_reliable_bit = inf
            next_least_reliable_bit = inf

            for recieved_value in self.recieved_values_dict.values():
                if abs(recieved_value) <= least_reliable_bit:
                    next_least_reliable_bit = least_reliable_bit
                    least_reliable_bit = recieved_value


            for bit_node, recieved_value in self.recieved_values_dict.items():

                send_value = least_reliable_bit                
                if recieved_value == least_reliable_bit:
                    send_value = next_least_reliable_bit

                bit_node.recieved_values_dict[self] = send_value

#Initialise parts 
detector = Detector()
message = Message([1,0])
graph = TannerGraph([[1,1],
                     [0,1]])

#Initialise the algorithm, reading the message into the graph
graph.init_message(message, detector)

#Print true and recieved message
print(message.bits)
print(graph.read_hard(message))

#Run algorithm unitl conditions are fullfilled 
while graph.correct_bit_ratio < 1 and graph.n_iterations < graph.max_iterations:

    graph.upward_pass()
    graph.downward_pass()
    graph.validate()

    graph.n_iterations += 1

#Print corrected message
print(graph.read_hard)

#%%
import numpy as np 

a = np.array([[1,1],
             [0,1]])

for r in a:
    print(r)