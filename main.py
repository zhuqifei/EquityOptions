"""
Created on Jul 17, 2018

@author: Qifei
"""

import os
import numpy as np
import pandas as pd

from VolSurface import VolSurface
from data.MarketDataEquityOptions import MarketdataEquityOptions
from pricer.BSPricer import BlackScholesSolver

if __name__ == '__main__':

    # 1. Load market data prices
    dataMS = MarketdataEquityOptions().getData('MS')

    # 2. Calibrate a surface
    surface = VolSurface.calibrate()


    exit