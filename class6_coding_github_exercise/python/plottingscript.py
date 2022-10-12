# Plotting packages
import matplotlib.pyplot as plt   # Plots
from matplotlib import rc         # Figure fonts
from OrthogDemo.py import OrthogDemo.py

# Choosing n and m
n,m =2,3
# Calling OrthogDemo function
orthog = OrthogDemo(n,m)

# Testing if the Integral is 0 or 1/2 
boolean= round(orthog.integral,2) == 0.0

#Creating plot
fig =plt.figure(figsize=(16,16))
rc('font',size=12)    
rc('font',weight='bold') 
rc('xtick',labelsize=12)  
rc('ytick',labelsize=12)

# ax1
ax = fig.add_subplot(311)
plt.title('n, m, integral',weight='bold',fontsize=12)
plt.plot(orthog.time,orthog.wave1,label='Wave 1')
plt.xlabel('Time',weight='bold',fontsize=12)
plt.legend(loc='best')
plt.grid()
# ax2
ax = fig.add_subplot(312)
plt.plot(orthog.time,orthog.wave2,label='Wave 2')
plt.xlabel('Time',weight='bold',fontsize=12)
plt.legend(loc='best')
plt.grid()
# ax3
ax = fig.add_subplot(313)
if boolean == True:
	string = 'Integral is 0'
plt.plot(orthog.time,orthog.product,label='Product')
plt.xlabel('Time',weight='bold',fontsize=12)
plt.legend(loc='best')
plt.grid()

plt.show()


