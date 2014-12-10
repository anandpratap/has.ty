from pylab import *

data = np.loadtxt("clustering.txt")

x = data[:,0];
y = data[:,1];
cx = data[:,2];
cy = data[:,3];
color = data[:,4];

scatter(x, y, c=color, marker="+");
plot(cx, cy, 'kx');
show()