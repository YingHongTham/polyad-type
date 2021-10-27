import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt


exp = pd.read_csv('Expression-matrix-Jan-2020.csv',index_col=0)
exp = exp.fillna(0) ##empty entries mean 0
