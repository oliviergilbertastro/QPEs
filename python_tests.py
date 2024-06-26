
if False:
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
    import astropy.io.fits as pyfits

    image = pyfits.open(f"data/science/tde0/rings.v3.skycell.2044.052.stk.r.unconv.fits")

    for i in image:
        print(i)
        print(i.header)
        print(i.data)
    plt.imshow(image[1].data)
    plt.show()

import matplotlib.pyplot as plt
import numpy as np
# example data
x = np.arange(0.1, 4, 0.5)
y = np.exp(-x)

# example error bar values that vary with x-position
error = 0.1 + 0.2 * x

fig, (ax0, ax1) = plt.subplots(nrows=2, sharex=True)
ax0.errorbar(x, y, yerr=error, fmt='-o')
ax0.set_title('variable, symmetric error')

# error bar values w/ different -/+ errors that
# also vary with the x-position
lower_error = 0.4 * error
upper_error = error
asymmetric_error = [lower_error, upper_error]

print(asymmetric_error)
ax1.errorbar(x, y, xerr=asymmetric_error, fmt='o')
ax1.set_title('variable, asymmetric error')
ax1.set_yscale('log')
plt.show()