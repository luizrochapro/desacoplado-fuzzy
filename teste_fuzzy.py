import numpy as np
import skfuzzy
import matplotlib.pyplot as plt


def first_nonzero(arr, axis, invalid_val=-1):
    mask = arr!=0
    return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)

def last_nonzero(arr, axis, invalid_val=-1):
    mask = arr!=0
    val = arr.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)

def first_one(arr, axis, invalid_val=-1):
    mask = arr==1
    val = arr.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)

# Set parameters 
N = 1001 
precisao = 1.0
minx = -1000
maxx = 1000 
#x = np.linspace(minx, maxx, N); 
x = np.linspace(minx, maxx, N)

# Create membership functions 
A = skfuzzy.trimf(x, np.array([320,350,370])) 
B = skfuzzy.trimf(x, np.array([0, 0, 0])) 

print(first_nonzero(A,0))
print(first_one(A,0))
print(last_nonzero(A,0))


# Apply fuzzy arithmetic 
z1, C1 = skfuzzy.dsw_add(x, A, x, B, 1000) 
z2, C2 = skfuzzy.dsw_sub(x, A, x, B, 31) 
z3, C3 = skfuzzy.dsw_mult(x, A, x, B, 31) 

# Plot results 
fig, (ax0, ax1, ax2) = plt.subplots(nrows=3)

ax0.plot(x,A,'--b', label='A') 
ax0.plot(x,B,':r', label='B') 
ax0.plot(z1,C1,'c', label='A+B') 
ax0.set_xlim(minx, maxx) 
ax0.legend() 
ax0.set_title("Fuzzy Addition, A+B") 

ax1.plot(x,A,'--b', label='A') 
ax1.plot(x,B,':r', label='B') 
ax1.plot(z2,C2,'c', label='A-B') 
ax1.set_xlim(minx, maxx) 
ax1.legend() 
ax1.set_title("Fuzzy Subtraction, A-B") 

ax2.plot(x,A,'--b', label='A') 
ax2.plot(x,B,':r', label='B') 
ax2.plot(z3,C3,'c', label='A*B') 
ax2.set_xlim(minx, maxx) 
ax2.legend() 
ax2.set_title("Fuzzy Multiplication, A*B") 

plt.tight_layout() 
plt.show()

'''
            #plotar  #################################
            fig, axs = plt.subplots(3,1)
            axs[0].plot(self.unidis,tri1)
            axs[1].plot(self.unidis,tri2) 
            axs[2].plot(self.pliq[k][0],self.pliq[k][1])
            axs[0].set_xlim(self.unidis.min(), self.unidis.max()) 
            axs[1].set_xlim(self.unidis.min(), self.unidis.max()) 
            axs[2].set_xlim(self.unidis.min(), self.unidis.max()) 
            fig.tight_layout()
            plt.show()
            ##########################################
'''