

import math
import sympy as symp
import numpy as nu
import matplotlib.pyplot as plt


plt.xlim(-65, 5)
#Step 1(draw poles)
zeros=[]
poles=[0, -25, complex(-50, 10), complex(-50, -10)]
S = symp.symbols('s')
Eqn= S ** 4 + 125 * S ** 3 + 5100 * S ** 2 + 65000 * S
plt.scatter([0, -25, -50, -50], [0, 0, 10, -10], color="red",
            marker= "X", s=150)
plt.title('Root locus')


# draw root locus branch (real axis locus)
plt.plot([0, -25], [0, 0], color="black", linewidth = 4)

# Step 2 (draw the asymptotes (dotted lines))
NumOfPoles=4
NumOfZeros=0
NumOfAsym = NumOfPoles - NumOfZeros
AngleOfAsym= 180 / NumOfAsym
IntersectionWithX= (poles[0] + poles[1] + poles[2] + poles[3]) / NumOfAsym
IntersectionWithX=IntersectionWithX.real
y = 70 * math.sin(math.radians(AngleOfAsym))
x = 70 * math.cos(math.radians(AngleOfAsym))

plt.plot([IntersectionWithX, IntersectionWithX + x], [0, +y], color='grey', linestyle='dotted', linewidth = 2)
plt.plot([IntersectionWithX, IntersectionWithX + x], [0, -y], color='grey', linestyle='dotted', linewidth = 2)
plt.plot([IntersectionWithX, IntersectionWithX - x], [0, +y], color='grey', linestyle='dotted', linewidth = 2)
plt.plot([IntersectionWithX, IntersectionWithX - x], [0, -y], color='grey', linestyle='dotted', linewidth = 2)




# Step 3 ( find break away points )
bpoint=[]
# here i make differentiation of the ch. equation with right to S
diffEqn=symp.diff(Eqn, S)
g = symp.Poly(diffEqn, S)
All = nu.roots(g.all_coeffs())
for x in range(len(All)):   # Accept only the real ones
    if(All[x].imag ==0):
        bpoint.append(All[x].real)

plt.scatter(bpoint, [0, 0, 0], color="blue",
            marker= "o", s=150)




# Step 4 (Find the intersections with the imaginary axis using Routh)
k = symp.symbols('K')
t= (4580 * 65000 - 125 * k) / 4580
All = symp.solve(t)
intersectwithimag=math.sqrt(All[0] / 4580)


plt.scatter([0, 0], [intersectwithimag, -intersectwithimag], color="green",
            marker= "x", s=150)





# Step 5 (Find Departure angles for complex poles)
setad= 180 - math.degrees(math.atan(10 / (-25))) - math.degrees(math.atan(10 / -50)) - 90
ey = 4 * math.sin(math.radians(180 - setad))
ex = 4 * math.cos(math.radians(180 - setad))
plt.plot([-50, -50 + 3], [10, 10], color='purple', linestyle='dashed', linewidth = 1)
plt.plot([-50, -50 + 3], [-10, -10], color='purple', linestyle='dashed', linewidth = 1)
plt.plot([-50, -50 - ex], [10, 10 + ey], color='purple', linestyle='dashed', linewidth = 1)
plt.plot([-50, -50 - ex], [-10, -10 - ey], color='purple', linestyle='dashed', linewidth = 1)
deg= str(round(setad, 5)) + '°'
plt.text(-48, 12, deg, fontsize=7)
plt.text(-48, -14, deg, fontsize=7)


#Step 6 ( plot root locus )
for k in nu.linspace(-20000, 20000000, 600):
  g = symp.Poly(Eqn + k, S)
  All = nu.roots(g.all_coeffs())
  for x in range(len(All)):
      plt.scatter([All[x].real], [All[x].imag], color="yellow", marker= "o", s=13)
                  
plt.show()
