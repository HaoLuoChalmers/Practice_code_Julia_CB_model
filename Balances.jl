include("Kinetics.jl")
include("Control.jl")
using LinearAlgebra
function Balances(t,x,dxdt_vector,DF)

  # Alias the species vector -
  e1 = x[1];
  e2 = x[2];
  e3 = x[3];
  e4 = x[4];
  e5 = x[5];
  e6 = x[6];
  e7 = x[7];
  M_glc = x[8];
  M_suc = x[9];
  M_for = x[10];
  M_lac = x[11];
  M_ace = x[12];
  M_eth = x[13];
  M_bio = x[14];

  #Define rate vector/Users/lhao/Documents/Git/ji_HCM-FBA-Model-Code-master/1AnaerobicCase/Kinetics.jl
  kinetics_dict = Kinetics(t,x,DF);
  rM = kinetics_dict["rM_vector"] #Metabolite Rate Vector
  rE = kinetics_dict["rE_vector"] #Enzyme Rate Vector
  rG = kinetics_dict["rG_vector"] #Growth Rate

  #Get Parameters from DF
  alpha = DF["EnzymeSynthesis"]
  beta = DF["Degradation"]
  Z = DF["ModeMatrix"]
  S = DF["MetaboliteMatrix"]
  (num_reations,num_modes) = size(Z)


  #Define Control Vector
  (u,v) = Control(t,x,rM,rE,DF)
  # println('u')
  # println(u)
  # println('v')
  # println(v)

  #V = eye(length(v)).*v
  #V = Array(Diagonal(fill(1.0, (length(v)))))*v
  V = Matrix{Float64}(I, length(v), length(v)).*v



  #Redine rates (rate*control)
  Growth_rate = rG.*v
  mu = sum(Growth_rate)


  #Enzyme Balances
  dxdt_vector[1] = alpha + rE[1]*u[1] - (beta+mu)*e1;
  dxdt_vector[2] = alpha + rE[2]*u[2] - (beta+mu)*e2;
  dxdt_vector[3] = alpha + rE[3]*u[3] - (beta+mu)*e3;
  dxdt_vector[4] = alpha + rE[4]*u[4] - (beta+mu)*e4;
  dxdt_vector[5] = alpha + rE[5]*u[5] - (beta+mu)*e5;
  dxdt_vector[6] = alpha + rE[6]*u[6] - (beta+mu)*e6;
  dxdt_vector[7] = alpha + rE[7]*u[7] - (beta+mu)*e7;

  # println("rM",rM)
  # println("M_bio",M_bio)
  #
  # #println(alpha)
  # #println(beta)
  # println("rE",rE)
  # println("u",u)
  # println("v",v)

  #
  # println("mu",mu)
  # println("e",x[1:6])
  #
  #
  # println('\n')

  #Metabolite Balances
  dxdt_vector[8:14] = S[1:7,:]*Z*V*rM*M_bio

  return dxdt_vector
end
