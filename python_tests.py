
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