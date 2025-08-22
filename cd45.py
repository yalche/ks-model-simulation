from protein import Protein
import math
import numpy as np

class CD45(Protein):
    LENGTH = 50  # Default length for TCR
    NAME = "CD45"
    K = 10  # Spring constant for compression energy calculation
    BINDING_ENERGY = -10.0
    BINDING_RANGE = [LENGTH - 10, LENGTH + 10] # Range for binding
    DIFFUSION_CONSTANT = 11000  # nm^2/s

    def __init__(self):
        """
        Initialize the TCR with a name and length.

        :param name: Name of the TCR.
        :param length: Length of the TCR.
        """
        super().__init__(self.NAME, self.LENGTH, self.K, self.DIFFUSION_CONSTANT)


    def get_bind_energy(self, h_next):
        """
        Calculate the binding energy for the TCR.

        :param h_prev: Previous height.
        :param h_next: Next height.
        :return: Binding energy of the TCR.
        """
        if h_next > self.BINDING_RANGE[0] or h_next < self.BINDING_RANGE[1]:
            return self.BINDING_ENERGY
        return 0.0

    def get_energy(self, h_next):
        """
        Calculate the energy for the TCR.

        :param h_prev: Previous height.
        :param h_next: Next height.
        :return: Energy of the TCR.
        """
        return self.get_compression_energy(h_next) + self.get_bind_energy(h_next)





