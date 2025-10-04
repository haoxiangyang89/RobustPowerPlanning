# test the sensitivity analysis

using Distributed;
addprocs(8);
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

# renewable sensitivity
for i in hData.hList
    hData.cost[i] = 6.25;
end
results_H = Dict();
for i in 1:8
    results_H[i] = solve_process(fData, uData, hData, bData, T, vmaxT, vminT, θDmaxT, θDminT, xLimit, groupDict[ci], 1, 3, 0.05*i,10);
end
save("Sensitivity_Test_H_series.jld","results",results_H);

# battery sensitivity
for i in hData.hList
    hData.cost[i] = 100;
end
bData_list = [];
for i in fData.IDList
    bData.cap[i] = 2;
    bData.cost[i] = 1.3;
end
for u in 1:10
    bData_copy = deepcopy(bData);
    for i in fData.IDList
        bData_copy.uCap[i] = (u-1)*0.2 + 0.1;
    end
    push!(bData_list, bData_copy);
end

for i in fData.IDList
    bData.uCap[i] = 0.5;
end
for u in 1:6
    bData_copy = deepcopy(bData);
    for i in fData.IDList
        bData_copy.cap[i] = (u-1)*0.2 + 0.5;
    end
    push!(bData_list, bData_copy);
end
results_B = pmap(i -> solve_process(fData, uData, hData, bData_list[i], T, vmaxT, vminT, θDmaxT, θDminT, xLimit, groupDict[ci], 2, 2, expansion_factor,4), 1:16);
save("Sensitivity_Test_B.jld","results",results_B);