from math import log10
import matplotlib.pyplot as plt
import numpy as np

signals_3 = np.array([2380749.470003745,36487.46171314418,65316.07025585285,563246.5695573135,179950.1537963712,80127.0630163737])
mags_3 = np.array([7.48,12.42,11.89,10.21,11,11.02])

signals_4 = np.array([1494987.5618783291,52529.82517136498,81846.70188109255,10228.746426281861,842598.8019151615,10691.346366423219])
mags_4 = np.array([8.21,11.85,11.25,13.09,9.61,13.06])

signals_6 = np.array([62708.46312241099,12408.604100044024,109841.52334134144,13389.80601068848,30548.733055026125,11466.245693156969])
mags_6 =

def logT(signals):
    for i in range(signals.size):
        signals[i] = log10(signals[i])
    return signals


x_3 = logT(signals_3)
y_3 = mags_3

m_3,b_3 = np.polyfit(x_3,y_3,1)

print(m_3,b_3)

plt.plot(x_3,y_3,"ro",x_3,m_3 * x_3 + b_3)

plt.show()
