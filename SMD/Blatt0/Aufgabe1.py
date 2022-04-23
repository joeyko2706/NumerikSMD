import numpy as np
import matplotlib.pyplot as plt

#Aufgabenteil a)
height, weight = np.genfromtxt('height_weight.txt', unpack = True)

bins1 = 5
bins2 = 10
bins3 = 15
bins4 = 20
bins5 = 30
bins6 = 50

#Histogramme für die Größe/Height mit verschiedenen Bins

plt.subplot(3,2,1)
plt.hist(height, label = '5 Bins', bins = bins1, color = 'm')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,2)
plt.hist(height, label = '10 Bins', bins = bins2, color = 'm')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,3)
plt.hist(height, label = '15 Bins', bins = bins3, color = 'm')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,4)
plt.hist(height, label='20 Bins', bins=bins4, color = 'm')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,5)
plt.hist(height, label = '30 Bins', bins = bins5, color = 'm')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,6)
plt.hist(height, label = '50 Bins', bins = bins6, color = 'm')#, histtype = 'step')
plt.legend(loc = 'best')

#plt.title('Histogramme für die Größe')
plt.tight_layout()
plt.savefig('HeightPlot.pdf')
plt.clf()
#-------------
#Histogramme für das Gewicht mit verschiedenen Bins

plt.subplot(3,2,1)
plt.hist(weight, label = '5 Bins', bins = bins1, color = 'lawngreen')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,2)
plt.hist(weight, label = '10 Bins', bins = bins2, color = 'lawngreen')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,3)
plt.hist(weight, label = '15 Bins', bins = bins3, color = 'lawngreen')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,4)
plt.hist(weight, label='20 Bins', bins=bins4, color = 'lawngreen')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,5)
plt.hist(weight, label = '30 Bins', bins = bins5, color = 'lawngreen')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,6)
plt.hist(weight, label = '50 Bins', bins = bins6, color = 'lawngreen')#, histtype = 'step')
plt.legend(loc = 'best')

plt.tight_layout()
plt.savefig('weightPlot.pdf')
plt.clf()


#Aufgabenteil c)
x = np.log(np.random.randint(1,101,10**5))

plt.subplot(3,2,1)
plt.hist(x, label = '5 Bins', bins = bins1, color = 'deeppink')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,2)
plt.hist(x, label = '10 Bins', bins = bins2, color = 'deeppink')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,3)
plt.hist(x, label = '15 Bins', bins = bins3, color = 'deeppink')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,4)
plt.hist(x, label='20 Bins', bins=bins4, color = 'deeppink')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,5)
plt.hist(x, label = '30 Bins', bins = bins5, color = 'deeppink')#, histtype = 'step')
plt.legend(loc = 'best')

plt.subplot(3,2,6)
plt.hist(x, label = '50 Bins', bins = bins6, color = 'deeppink')#, histtype = 'step')
plt.legend(loc = 'best')

plt.tight_layout()
plt.savefig('Aufgabec.pdf')