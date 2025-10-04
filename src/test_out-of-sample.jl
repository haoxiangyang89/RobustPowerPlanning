# test different Γ
using Distributed;
addprocs(4);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

@everywhere caseList = [14,118]; 
@everywhere ci = 1;
@everywhere T = 4;

αmax = [0.128,0.140];
dαbot = αmax[ci];
dαtop = 0.2*dαbot;
hαtop = 0.2;
hαbot = 0.2;
expansion_factor = 0.2;

fData, bData, hData, dData = readInData(ci,caseList,1e4,1,24);
uData = makeUAData(fData, hData, dData, T, dαtop, dαbot, hαtop, hαbot);
groupData = load("../data/groupDict.jld");
groupDict = groupData["groupDict"];
xLimit = 1;

# initialize the bounds
vmaxT = Dict();
vminT = Dict();
θDmaxT = Dict();
θDminT = Dict();
for t in 1:T
    vmaxT[t] = Dict();
    vminT[t] = Dict();
    θDmaxT[t] = Dict();
    θDminT[t] = Dict();
    for i in fData.IDList
        vmaxT[t][i] = fData.Vmax[i];
        vminT[t][i] = fData.Vmin[i];
    end
    for k in fData.brList
        θDmaxT[t][k] = pi/6;
        θDminT[t][k] = -pi/6;
    end
end

gamma_data = load("Gamma_Test.jld");
gamma_results = Dict();
for gd in 1:3
    for gh in 1:3
        gamma_results[gd,gh] = gamma_data["results"][gd,gh];
    end
end
gamma_0_data = load("deterministic_solution.jld");
gamma_results[0,0] = gamma_0_data["results"];

Γ = Dict("d" => 1, "h" => 1);

oos_results = pmap(i -> test_second(fData,uData,hData,bData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,gamma_results), 1:1000);

save("OoS_Test.jld","results",oos_results);