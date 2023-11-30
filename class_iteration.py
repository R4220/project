# class_iteration.py

from class_group import group

class iteration:
    n_at = 0
    ax = [0, 0, 0]
    ay = [0, 0, 0]
    az = [0, 0, 0]
    t = 0

    block = []

    pass

    def first_two_lines(self):
        return f'{self.n_at}\nLattice = \"{self.ax[0]}, {self.ax[1]}, {self.ax[2]}, {self.ay[0]}, {self.ay[1]}, {self.ay[2]}, {self.az[0]}, {self.az[1]}, {self.az[2]}\"'
    
    def store(self, line):
        self.block = self.block.append(line)
        return
    
    def single_frame(self):
        return