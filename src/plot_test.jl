
using PyPlot

x = range(0, 2π, length=20)
y = cos.(x)

p = scatter(x, y)
show(p)
#savefig("myPPplot.png", dpi=800)