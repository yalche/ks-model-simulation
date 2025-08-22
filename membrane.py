import numpy as np
from distributor import Distributor
from tcr import TCR
from cd45 import CD45

class TCellMembrane:
    K = 25  # Membrane bending rigidity in kBT

    def __init__(self, size):
        self._default_inter_membrane = 70
        self._size = size
        self._grid = [[0 for _ in range(size)] for _ in range(size)]
        self._grid_inter_membranes = np.ones((self._size, self._size)) * self._default_inter_membrane
        self._cell_size = 10
        self._proteins = []


    def init_grid(self, cd45_num, tcr_num, distributor, radius=None, center=None):
        tcr_coords = Distributor.circle(self._grid, n_molecules=tcr_num, radius=self._size / 16)
        for x, y in tcr_coords:
            self._grid[int(x)][int(y)] = TCR()
            self._proteins.append(self._grid[int(x)][int(y)])
            #self._grid_height[int(x)][int(y)] = TCR.LENGTH
        cd45_coords = []
        if distributor == "uniform":
            cd45_coords = Distributor.uniform(self._grid, n_molecules=cd45_num)
        elif distributor == "circle":
            cd45_coords = Distributor.circle(self._grid, cd45_num, center, radius)
        elif distributor == "semicircle":
            cd45_coords = Distributor.semicircle(self._grid, cd45_num, center, radius)
        for x, y in cd45_coords:
            self._grid[int(x)][int(y)] = CD45()
            self._proteins.append(self._grid[int(x)][int(y)])
            self._grid_inter_membranes[int(x)][int(y)] = CD45.LENGTH

    def get_location_after_step(self, x, y, dx, dy):
        new_x = (x + dx // self._cell_size).astype(int) % self._size
        new_y = (y + dy // self._cell_size).astype(int) % self._size
        return new_x, new_y

    def get_molecule_at(self, x, y):
        if 0 <= x < self._size and 0 <= y < self._size:
            return self._grid[int(x)][int(y)]
        return None

    def set_molecule_at(self, x, y, new_x, new_y):
        if self._grid[int(new_x)][int(new_y)] == 0:
            self._grid[int(new_x)][int(new_y)] = self._grid[int(x)][int(y)]
            self._grid[int(x)][int(y)] = 0

    def get_proteins(self):
        return self._proteins

    def get_inter_membrane_at(self, x, y):
        return self._grid_inter_membranes[x][y]

    def set_height_at(self, x, y, height):
        self._grid_inter_membranes[int(x)][int(y)] = height

    def get_bending_energy(self, x, y, new_height):
        h = new_height
        h_up = self._grid_inter_membranes[int((x - 1) % self._size)][int(y)]
        h_down = self._grid_inter_membranes[int((x + 1) % self._size)][int(y)]
        h_left = self._grid_inter_membranes[int(x)][int((y - 1) % self._size)]
        h_right = self._grid_inter_membranes[int(x)][int((y + 1) % self._size)]
        bending_energy = (h_up + h_down + h_left + h_right - 4 * h) ** 2 * self.K / (2 * self._cell_size ** 2)
        return bending_energy

    def get_inter_membrane_step(self):
        # dz - random change in height - Gaussian with mean 0 and std 1
        return np.random.normal(0, 1, (self._size, self._size))


class APCMembrane:
    def __init__(self, size):
        self._default_inter_membrane = 70
        self._size = size
        self._grid = np.zeros((size, size))

    def init_grid(self, pmhc_num):
        pmhc_coords = Distributor.uniform(self._grid, n_molecules=pmhc_num)
        for x, y in pmhc_coords:
            self._grid[int(x)][int(y)] = 1

    def get_molecule_at(self, x, y):
        if 0 <= x < self._size and 0 <= y < self._size:
            return self._grid[int(x)][int(y)]
        return None


