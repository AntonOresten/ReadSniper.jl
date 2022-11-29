
using PyPlot
using CSV, DataFrames

df = CSV.read("test/output/coverage.csv", DataFrame)

x, y = df.i, df.frequency
n = length(x)

bar(x, y, width=1)
show()