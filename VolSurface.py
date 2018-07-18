import pandas as pd

class VolSurface:

    def __init__(self):
        """
        Constructor
        """

    def __init__(self, rowLabels, columnLabels, values):
        self.values = values
        self.rowLabels = rowLabels
        self.columnLabels = columnLabels
        self.surface = pd.DataFrame(data=self.values, index=self.rowLabels, columns=self.columnLabels)

    # static method to calibrate a surface given market data object
    @staticmethod
    def calibrate():
        return NotImplemented
