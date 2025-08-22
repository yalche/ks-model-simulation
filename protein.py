import math

import numpy as np


class Protein:
    def __init__(self, name, length, k, diffusion):
        """
        Initialize the Protein with a name and sequence.
        :param name: Name of the protein.
        :param sequence: Amino acid sequence of the protein.
        """
        self._name = name
        self._length = length
        self._k = k
        self._diffusion = diffusion
        self._location = (0, 0)

    def set_bound(self):
        pass


    def get_diffusion(self):
        return self._diffusion

    def get_location(self):
        return self._location

    def set_location(self, location):
        self._location = location

    def get_length(self):
        return self._length

    def get_name(self):
        return self._name

    def get_compression_energy(self, h_next):
        """
        Calculate the compression energy for the TCR.

        :param h_prev: Previous height.
        :param h_next: Next height.
        :return: Compression energy.
        """
        if h_next < self._length:
            return 0
        return 0.5 * (h_next - self._length) ** 2 * self._k

    @staticmethod
    def get_random_step(molecules_num, diffusion_constant):
        thetas = np.random.uniform(0, 2 * np.pi, size=molecules_num)
        steps = np.abs(np.random.normal(0, np.sqrt(0.01*diffusion_constant*4), size=molecules_num))
        dx = np.cos(thetas) * steps
        dy = np.sin(thetas) * steps
        return dx, dy

    def get_energy(self,h_next):
        """
        Calculate the energy for the protein.

        :param h_prev: Previous height.
        :param h_next: Next height.
        :return: Energy of the protein.
        """
        raise NotImplementedError("This method should be implemented in subclasses.")