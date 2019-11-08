include("loadData.jl")
include("SolveBalances.jl")

#DefineTime
tStart = 0.0;
tStop = 8.0;
tStep = 0.1;

#loadParameterValues
Param = Float64[]

print("model")
#SolverResults
(t, x) = SolveBalances(tStart, tStop, tStep, Param)

#Define species from solver
xGlc = x[:,7 + 1];
xSuc = x[:,8 + 1];
xFor = x[:,9 + 1];
xLac = x[:,10 + 1];
xAce = x[:,11 + 1];
xEth = x[:,12 + 1];
xBio = x[:,13 + 1];
#
#
using PyPlot
figure(1)
plot(t, xGlc, color = "k", label = "HCM FBA", linewidth = 2)
plot(t, xSuc, color = "k", linewidth = 2)
plot(t, xFor, color = "k", linewidth = 2)
plot(t, xLac, color = "k", linewidth = 2)
plot(t, xAce, color = "k", linewidth = 2)
plot(t, xEth, color = "k", linewidth = 2)
plot(t, xBio, color = "k", linewidth = 2)

legend(fontsize = 18)
xlabel("Time (hr)", fontsize = 20)
ylabel("Abundance (mM)", fontsize = 20)
#axis([0,8,0,2.1;])
display(gcf())
