import numpy as np

from protein import Protein
import math

class TCR(Protein):
    LENGTH = 13  # Default length for TCR
    NAME = "TCR"
    K = 10.0  # Spring constant for compression energy calculation
    DIFFUSION_CONSTANT = 10000  # nm^2/s

    def __init__(self):
        """
        Initialize the TCR with a name and length.

        :param name: Name of the TCR.
        :param length: Length of the TCR.
        """
        super().__init__(self.NAME, self.LENGTH, self.K, self.DIFFUSION_CONSTANT)
        self._is_bound = False


    def set_bound(self):
        """
        Set the TCR to a bound state.
        This method can be extended to include specific behavior for bound state.
        """
        self._bound = True

    def get_energy(self, h_next):
        """
        Calculate the energy for the TCR.

        :param h_prev: Previous height.
        :param h_next: Next height.
        :return: Energy of the TCR.
        """
        if self._is_bound:
            return -math.inf
        return self.get_compression_energy(h_next)





