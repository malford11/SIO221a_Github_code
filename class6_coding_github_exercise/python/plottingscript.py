# Plotting packages
import matplotlib.pyplot as plt   # Plots
from matplotlib import rc         # Figure fonts
from OrthogDemo import OrthogDemo

# Choosing n and m
n,m =2,3
# Calling OrthogDemo function
orthog = OrthogDemo(n,m)

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
plt.title('Integral is '+str(round(orthog.integral,1)),weight='bold',fontsize=12)
plt.plot(orthog.time,orthog.product,label='Product')
plt.xlabel('Time',weight='bold',fontsize=12)
plt.legend(loc='best')
plt.grid()

plt.show()


