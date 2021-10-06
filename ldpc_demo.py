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

    class bit_node:
        pass 

    class check_node:
        pass 

