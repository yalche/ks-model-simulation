import math

import numpy as np
from membrane import TCellMembrane, APCMembrane
from protein import Protein


class Simulation:
    """
    A class to represent a simulation.
    """

    def __init__(self, ):
        """
        Initialize the simulation with default parameters.
        """
        self._time_step = 0.01
        self._total_time = 100.0
        self._grid_size = 200
        self._cd45_num = 500
        self._tcr_num = 125
        self._pmhc_num = 300
        #init grid
        self._membrane = TCellMembrane(self._grid_size)
        self._membrane.init_grid(self._cd45_num, self._tcr_num, distributor="uniform")
        self._energy_grid = np.zeros((self._grid_size, self._grid_size))
        self._apc_membrane = APCMembrane(self._grid_size)
        self._apc_membrane.init_grid(self._pmhc_num)

    def _is_step_legal(self, x, y, new_x, new_y, i):
        for j in range(len(x)):
            if i != j and new_x[i] == new_x[j] and new_y[i] == new_y[j]:
                return False
            if i != j and x[j] == new_x[i] and y[j] == new_y[i]:
                return False
        return True

    def _get_step_new_energy(self, protein, x, y, new_x, new_y):
        return [self._energy_grid[x][y] - protein.get_energy(self._membrane.get_inter_membrane_at(x,y)),
                self._energy_grid[new_x][new_y] + protein.get_energy(self._membrane.get_inter_membrane_at(new_x,new_y))]

    def _tcr_binding(self, protein, x, y, new_x, new_y):
        self._energy_grid[x][y] -= protein.get_energy(self._membrane.get_inter_membrane_at(x,y))
        self._energy_grid[new_x][new_y] = -math.inf
        self._membrane.set_height_at(new_x, new_y, protein.get_length())
        directions = [(-1, 0), (0, -1), (0, 1), (1, 0)]
        for d in directions:
            adj_x, adj_y = new_x + d[0], new_y + d[1]
            if self._membrane.get_molecule_at(adj_x, adj_y) == 0:
                self._energy_grid[adj_x][adj_y] = self._membrane.get_bending_energy(adj_x, adj_y, self._membrane.get_inter_membrane_at(adj_x, adj_y))
            else:
                self._energy_grid[adj_x][adj_y] = (self._membrane.get_molecule_at(adj_x, adj_y).get_energy(self._membrane.get_inter_membrane_at(adj_x, adj_y))
                                                   + self._membrane.get_bending_energy(adj_x, adj_y, self._membrane.get_inter_membrane_at(adj_x, adj_y)))


    def _get_membrane_new_energy(self, i, j, dz):
        new_height = self._membrane.get_inter_membrane_at(i, j) + dz
        if self._membrane.get_molecule_at(i, j) != 0:
            new_energy = self._membrane.get_molecule_at(i, j).get_energy(new_height) + self._membrane.get_bending_energy(i, j, new_height)
        else:
            new_energy = self._membrane.get_bending_energy(i, j, new_height)
        return new_energy

    def _is_step_accepted(self, delta_e):
        if delta_e < 0:
            return True
        else:
            prob = np.exp(-delta_e)
            rand = np.random.uniform(0, 1)
            if rand < prob:
                return True
        return False

    def run(self):
        """
        Run the simulation.
        """
        num_steps = int(self._total_time / self._time_step)
        proteins = self._membrane.get_proteins()
        diffusion_constants = np.array([protein.get_diffusion() for protein in proteins])

        for i in range(self._grid_size):
            for j in range(self._grid_size):
                if self._membrane.get_molecule_at(i, j) != 0:
                    self._energy_grid[i][j] = self._membrane.get_molecule_at(i, j).get_energy(
                        self._membrane.get_inter_membrane_at(i, j))
                else:
                    self._energy_grid[i][j] = self._membrane.get_inter_membrane_at(i, j)

        for protein in proteins:
            x_i, y_i = protein.get_location()
            if protein.get_name() == "TCR" and self._apc_membrane.get_molecule_at(x_i, y_i) != 0:
                print("Initial binding at ({}, {})".format(x_i, y_i))
                self._tcr_binding(protein, x_i, y_i, x_i, y_i)
                protein.set_bound()

        for step in range(num_steps):
            dz = self._membrane.get_inter_membrane_step()
            for i in range(self._grid_size):
                for j in range(self._grid_size):
                    new_energy = self._get_membrane_new_energy(i, j, dz[i][j])
                    delta_energy = new_energy - self._energy_grid[i][j]
                    if self._is_step_accepted(delta_energy):
                        new_height = self._membrane.get_inter_membrane_at(i, j) + dz[i][j]
                        self._energy_grid[i][j] = new_energy
                        self._membrane.set_height_at(i, j, new_height)

            dx, dy = Protein.get_random_step(len(proteins), diffusion_constants)
            # new list of new positions
            x = [protein.get_location()[0] for protein in proteins]
            y = [protein.get_location()[1] for protein in proteins]

            new_x, new_y = self._membrane.get_location_after_step(x, y, dx, dy)
            for i, protein in enumerate(proteins):
                x_i, y_i = protein.get_location()
                new_x_i, new_y_i = new_x[i], new_y[i]
                if self._is_step_legal(x, y, new_x, new_y, i):
                    if protein.get_name() == "TCR" and self._apc_membrane.get_molecule_at(new_x_i, new_y_i) != 0:
                        self._tcr_binding(protein, x_i, y_i, new_x_i, new_y_i)
                        protein.set_location((new_x_i, new_y_i))
                        self._membrane.set_molecule_at(x_i, y_i, new_x_i, new_y_i)
                        protein.set_bound()
                    else:
                        new_energy = self._get_step_new_energy(protein, x_i, y_i, new_x_i, new_y_i)
                        delta_energy = new_energy[0] - self._energy_grid[x_i][y_i] + new_energy[1] - self._energy_grid[new_x_i][new_y_i]
                        if self._is_step_accepted(delta_energy):
                            protein.set_location((new_x_i, new_y_i))
                            self._membrane.set_molecule_at(x_i, y_i, new_x_i, new_y_i)
                            self._energy_grid[x_i][y_i] = new_energy[0]
                            self._energy_grid[new_x_i][new_y_i] = new_energy[1]

            print("Step {}/{} completed".format(step + 1, num_steps))
            print([protein.get_location() for protein in proteins])

def main():
    sim = Simulation()
    sim.run()

if __name__ == "__main__":
    main()




