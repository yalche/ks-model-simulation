import numpy as np

CD45_NUM = 1000
TCR_NUM = 1000
grid_NUM = 400

class Simulation:
    """
    A class to represent a simulation.
    """

    def __init__(self, cd45_num=CD45_NUM, tcr_num=TCR_NUM, grid_num=grid_NUM):
        """
        Initialize the simulation with default parameters.
        """
        self.cd45_num = cd45_num
        self.tcr_num = tcr_num
        self.grid = np.zeros((grid_num, grid_num))
        self.step_size = 1
        self.grid_size = grid_num
        self.cd45, self.tcr = self.initialize_positions()

    def initialize_positions(self):
        """
        Initialize the positions of CD45 and TCR molecules randomly within the grid.
        """
        x, y = np.meshgrid(np.arange(self.grid_size), np.arange(self.grid_size))
        all_positions = np.vstack([x.ravel(), y.ravel()]).T  # shape: (grid_size^2, 2)

        np.random.shuffle(all_positions)
        cd45 = all_positions[:self.cd45_num].T
        tcr = all_positions[self.cd45_num:self.cd45_num + self.tcr_num].T

        return cd45, tcr

    def step(self):
        """
        Simulate a step in the protein dynamics.
        """
        occupied_locations = set()

        # Move CD45 molecules
        cd45_new_locations = self.move_molecules(self.cd45, occupied_locations, molecule_type=1)
        # Move TCR molecules
        tcr_new_locations = self.move_molecules(self.tcr, occupied_locations, molecule_type=2)

        return self.grid

    def move_molecules(self, molecules_locations: np.array, occupied_locations, molecule_type):
        """
        Move molecules to new locations while avoiding occupied positions.
        """
        new_locations = np.copy(molecules_locations)
        for i in range(len(molecules_locations)):
            new_x, new_y = None, None
            while (new_x, new_y) in occupied_locations or new_x is None or new_y is None:
                # Randomly choose a direction to move
                dx, dy = np.random.choice([-1, 0, 1], size=2)
                new_x = (molecules_locations[0,i] + dx) % self.grid_size
                new_y = (molecules_locations[1,i] + dy) % self.grid_size
            occupied_locations.add((new_x, new_y))
            new_locations[:, i] = (new_x, new_y)
        return new_locations




