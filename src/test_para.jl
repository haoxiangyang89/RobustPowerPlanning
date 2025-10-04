using Distributed;
addprocs(32);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

@everywhere caseList = [14,118]; 
@everywhere ci = 1;
@everywhere T = 6;

αmax = [0.128,0.140];
dαbot = αmax[ci];
dαtop = 0.2*dαbot;
hαtop = 0.2;
hαbot = 0.2;
expansion_factor = 0.2;

fData, bData, hData, dData = readInData(ci,caseList,1e4,1,16);
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

ΓSet = Iterators.product(1:3,1:3);

k_list = [k for k in fData.brList if (k != (7,8,1))&&(k[1] < k[2])];
t_list = 2:6;
outDict = solve_process_para(fData, uData, hData, bData, T, vmaxT, vminT, θDmaxT, θDminT, xLimit, groupDict[ci], 3, 3, expansion_factor, k_list[1], t_list[1]);
