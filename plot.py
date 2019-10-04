import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
data = pd.read_csv('Convergence.csv',header=0)
style.use('ggplot')
y = data['Pbest']
x = data['It']
plt.plot(x,y)
plt.title('Convergence')
plt.xlabel('Iterations')
plt.ylabel('Pbest')
plt.show()
