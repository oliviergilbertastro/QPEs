

func = lambda x1, x2: 4*x1+x2

print(func(2,3))
print((lambda x1, x2: 4*x1+x2)(4,5))

class EmptyClass():
    def __init__(self, var):
        self.var = var

test = EmptyClass({"potato": 1})
var = test.var
res = var["potato"]
var2 = var

test.var["potato"] = 2
print(var)
print(var2)
print(res)

print("\x1b[33mYELLOW\x1b[0m")


import matplotlib.pyplot as plt
x = 10
y = 4

# Define the uncertainties in x and y directions
xerr = [2, 4]  # This means x error ranges from x-2 to x+4
yerr = [0.5, 0.3]  # This means y error ranges from y-0.5 to y+0.3

# Plot the point with error bars
plt.errorbar(x, y, xerr=[[xerr[0]], [xerr[1]]], yerr=[[yerr[0]], [yerr[1]]], fmt='o', capsize=5)

plt.show()