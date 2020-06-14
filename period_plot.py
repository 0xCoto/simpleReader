import numpy as np
import matplotlib.pyplot as plt
import csv

plt.rcParams["figure.figsize"] = (17,10)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

xx=[]
y=[]
with open('results.txt', 'r') as csvfile:
    plots= csv.reader(csvfile, delimiter=',')
    for row in plots:
        xx.append(float(row[2])*1000)
        y.append(float(row[3]))

plt.plot(xx,y, c='royalblue', label='Measured data (TLM-18)', marker='o')
plt.ticklabel_format(useOffset=False)

plt.axvline(x=xx[y.index(max(y))], color='brown', linestyle='--', linewidth=2)
plt.text(xx[y.index(max(y))]+0.02, max(y)-1, str('%f' % xx[y.index(max(y))])+' ms', color = 'brown', fontsize=16)

#plt.text(271, 2, '(Observation data have been incoherently\n de-dispersed with a DM of $26.7787 \ pc/cm^{3})$', color = 'gray', fontsize=14)

plt.title('$P_{0}-\mathrm{Search}$ | $'+str(xx[0])+'-'+str(xx[-1])+'\mathrm{ \ ms \ (Step: \ }\Delta P = '+str(float('%f' % (xx[1]-xx[0])))+' \mathrm{ \ ms})$', fontsize=20, y=1.01)

plt.xlabel('$\mathrm{Period \ } (P_{0})$', fontsize=16)
plt.ylabel('$\mathrm{Signal-to-Noise \ Ratio \ (S/N)}$', fontsize=16)

plt.xlim(min(xx),max(xx))

#plt.legend(bbox_to_anchor=(1,1), loc="upper right")
#plt.xticks(rotation=90)

plt.grid()

plt.savefig('psearch.png', bbox_inches='tight')
