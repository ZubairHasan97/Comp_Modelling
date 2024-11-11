import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Modelling 3-gene Oscillatory_Network: w/Negative Feedback response

# initial state
print("Enter intital state of 3 genes: e.g. 2 3 4")
y0 = [float(x) for x in input().split()]
# print(y0)
t = np.linspace(0,200, num=100)

# key parameters
k_1 = float(input("production_rate_G1 (K1): "))
gamma_1 = float(input("degradation_rate_G1 (gamma_1): "))
k_2 = float(input("production_rate_G2 (K2): "))
gamma_2 = float(input("degradation_rate_G2 (gamma_2): "))
k_3 = float(input("production_rate_G3 (K3): "))
gamma_3 = float(input("degradation_rate_G3 (gamma_3): "))
n = float(input("constant n: (e.g. within 1-9) "))
c = float(input("constant C: "))
param = [k_1, gamma_1, k_2, gamma_2, k_3, gamma_3, n, c]

# simlation function
def sim3gene(init_var, timepoints, params):

  g1 = init_var[0]
  g2 = init_var[1]
  g3 = init_var[2]

  #params
  k1 = params[0]
  gam1 = params[1]
  k2 = params[2]
  gam2 = params[3]
  k3 = params[4]
  gam3 = params[5]
  n = params[6]
  c = params[7]

  #ode system for g1,g2,g3
  dg1dt = ((c**n / (c**n + g3**n)) * k1) - (gam1 * g1) # gene1 repressed by the expression of gene3 -> repression hill
  dg2dt = ((g1**n / (c**n + g1**n)) * k2) - (gam2 * g2) # gene2 gets activated by the expression of gene1 -> activation hill function
  dg3dt = ((g2**n / (c**n + g2**n)) * k3) - (gam3 * g3) # gene3 gets activated by the expression of gene2 -> activation hill function

  return dg1dt, dg2dt, dg3dt

# invoke ode solver
y = odeint(sim3gene, y0, t, args=(param,))

# printing out the defferential results of 3 genes
print(y)

# plotting the results
fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=False)

ax1.plot(t, y[:,0], label="Gene1", color="b")
ax2.plot(t, y[:,1], label="Gene2", color="r")
ax3.plot(t, y[:,2], label="Gene3", color="g")

ax1.set_xlabel("time")
ax1.set_ylabel("G1-exp-lvl")

ax2.set_xlabel("time")
ax2.set_ylabel("G2-exp-lvl")

ax3.set_xlabel("time")
ax3.set_ylabel("G3-exp-lvl")

ax1.legend()
ax2.legend()
ax3.legend()

plt.tight_layout()
plt.show()
