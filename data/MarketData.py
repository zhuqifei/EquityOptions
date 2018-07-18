from abc import ABC, abstractmethod


class MarketData(ABC):

    def __init__(self):
        """
        Constructor
        """

    @abstractmethod
    def getData(self, ticker):
        return NotImplemented
