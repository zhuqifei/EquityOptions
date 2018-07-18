from abc import ABC, abstractmethod


class Solver(ABC):

    def __init__(self):
        """
        Constructor
        """

    @abstractmethod
    def solve(self):
        return NotImplemented
