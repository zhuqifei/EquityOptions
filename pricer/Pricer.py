from abc import ABC, abstractmethod


class Pricer(ABC):

    def __init__(self):
        """
        Constructor
        """

    @abstractmethod
    def price(self):
        return NotImplemented
