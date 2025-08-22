import random

class Distributor:
    @staticmethod
    def uniform(grid, n_molecules):
        empty_coords = [(x, y) for y, row in enumerate(grid) for x, val in enumerate(row) if val == 0]

        if n_molecules > len(empty_coords):
            raise ValueError("Too many molecules for available empty cells")

        return random.sample(empty_coords, n_molecules)

    @staticmethod
    def circle(grid, n_molecules, center=None, radius=None):
        height = len(grid)
        width = len(grid[0])

        if center is None:
            center = (width // 2, height // 2)
        if radius is None:
            radius = min(width, height) // 4

        empty_coords = [
            (x, y) for y in range(height) for x in range(width)
            if grid[y][x] == 0 and (x - center[0])**2 + (y - center[1])**2 <= radius**2
        ]

        if n_molecules > len(empty_coords):
            raise ValueError("Too many molecules for circle area")

        return random.sample(empty_coords, n_molecules)

    @staticmethod
    def semicircle(grid, n_molecules, center=None, radius=None):
        height = len(grid)
        width = len(grid[0])

        if center is None:
            center = (width // 2, height // 2)
        if radius is None:
            radius = min(width, height) // 4

        empty_coords = [
            (x, y) for y in range(height) for x in range(width)
            if grid[y][x] == 0
            and (x - center[0])**2 + (y - center[1])**2 <= radius**2
            and y >= center[1]
        ]

        if n_molecules > len(empty_coords):
            raise ValueError("Too many molecules for semicircle area")

        return random.sample(empty_coords, n_molecules)
