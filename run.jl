include("loadData.jl")
include("SolveBalances.jl")

#DefineTime
tStart = 0.0;
tStop = 8.0;
tStep = 0.1;

#loadParameterValues
Param = Float64[]

#SolverResults
(t,x) = SolveBalances(tStart,tStop,tStep,Param)

#Define species from solver
xGlc = x[:,7]
xSuc = x[:,8]
xFor = x[:,9]
xLac = x[:,10]
xAce = x[:,11]
xEth = x[:,12]
xBio = x[:,13]


using PyPlot
figure(1)
plot(t,xGlc,color="k",label="HCM FBA",linewidth=2)
plot(t,xSuc, color="k",linewidth=2)
plot(t,xFor, color="k",linewidth=2)
plot(t,xLac, color="k",linewidth=2)
plot(t,xAce, color="k",linewidth=2)
plot(t,xEth, color="k",linewidth=2)
plot(t,xBio, color="k",linewidth=2)

legend(fontsize=18)
xlabel("Time (hr)",fontsize=20)
ylabel("Abundance (mM)",fontsize=20)
#axis([0,8,0,2.1;])
gcf()
