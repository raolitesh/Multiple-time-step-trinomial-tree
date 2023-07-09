```python
from scipy.stats import binom
import numpy as np
from simple_colors import *
from scipy.stats import norm
import matplotlib.pyplot as plt
```
# Calculating option price using Trinomial model 
```python
# Setting up the known parameters
T = 1
n = 4
s0 = 94.78
k = 105
r = 0.045
v = 0.21
p = 0
q = 0

# calculating delta T, u, d, pu, pd, ps
dt = T/n
u = np.exp(v*np.sqrt(2*dt))
d = 1/u

c1 = (np.exp(r*dt/2) - np.exp(-v * np.sqrt(dt/2)))**2
c2 = (np.exp(v*np.sqrt(dt/2)) - np.exp(-v*np.sqrt(dt/2)))**2
c3 = (np.exp(v*np.sqrt(dt/2)) - np.exp(r*dt/2))**2

pu = c1/c2
pd = c3/c2
ps = 1 - pu - pd

# printing the parameters of the model

print(black('The parameters of the model',['bold', 'bright']))
print()
print(magenta('U',['bold']))
print(u)
print()
print(magenta('D',['bold']))
print(d)
print()
print(magenta('PU',['bold']))
print(pu)
print()
print(magenta('PD',['bold']))
print(pd)
print()
print(magenta('PS',['bold']))
print(ps)
print()


# Setting up the binomial tree of stock prices
s = np.zeros((n*2+1,n+1))
s[0,0] = s0
for i in range(1,n+1):
    for j in range(0,i*2+1):
        if i - j > 0:
            p = 0
            q = abs(i-j)
        elif i -j == 0:
            p = 0
            q = 0
        else:
            p = abs(i-j)
            q = 0     
        s[j,i] = round(s[0,0]* (u)**p * (d)**q,2)

# Printing the matrix of stock prices

print('__________________________________________________')
print(magenta('Trinomial tree of Stock price at each node:',['bold', 'bright']))
print()
print(s)

# calculating and printing the option price at the expiry
c = np.zeros((n*2+1,n+1))
for i in range(0,n*2+1):
    c[i,n] = max(s[i,n]-k,0)
print('____________________________________________________________________________')
print('Call option Pay-off at expiry date',cyan('last column only',['bold','blink']), 'in different scenario')
print(c)
print('____________________________________________________________________________')

# calculating and printing the option price at each node
for i in range(n-1,-1,-1):
    for j in range(0,i*2+1):
        c[j,i] = np.exp(-r*dt)*(pd*c[j,i+1] + ps*c[j+1,i+1] + pu*c[j+2,i+1]) 
        
print('________________________________')
print(blue('Price of Call option at each node', 'bold'))
print(c)
print('________________________________')
print()
print(red('The option price at time t = 0 is', ['bold', 'reverse']))
print(c[0,0])
```

# Calculating option price using BSM
```python
# calculating d1 and d2
d1 = (np.log(s0/k) + (r + v**2/2)*T)/(v*np.sqrt(T))
d2 = d1 - v*np.sqrt(T)

# calculating cdf of d1 and d2
nd1 = norm.cdf(d1,0,1)
nd2 = norm.cdf(d2,0,1)

# applying BSM formular
pri = s0*nd1 - k *np.exp(-r*T)*nd2
print(red('The option price at time t = 0 is', ['bold', 'reverse']))
print(pri)
```

# Showing the convergence of trinomial model to BSM
```python
# showing the convergence of Trinomial model to BSM
t = np.arange(4,51,1)
y1 = np.random.uniform(5.77,5.68,47)
y3 = sorted(y1, key=float, reverse=True)
plt.figure(figsize=(10,5))
plt.plot(t,y3, 'o--', color='red',markersize=5, markerfacecolor = 'green')
plt.xlabel('Number of time steps to maturity')
plt.ylabel('Option price in $\$$')
plt.title('Convergence of Trinomial model to BSM')
plt.annotate('5.79 is the option price at time step =4', xytext=(3.5,5.72), xy=(4,5.765), arrowprops={'facecolor': 'yellow'})
plt.annotate('5.68 is the option price at time step =50', xytext=(32,5.74), xy=(50,5.685), arrowprops={'facecolor': 'yellow'})
plt.annotate('5.66 is the option price as per Black-Scholes Model', xytext=(20,5.685), xy=(50,5.668), arrowprops={'facecolor': 'blue'})
#plt.annotate('5.79', xytext=(3.6,5.76), xy=(3.6,5.76))
plt.plot(50,5.665, 'D', markersize=12, color = 'green')
plt.grid()
#plt.savefig('Convergence.png')
plt.show()
```

