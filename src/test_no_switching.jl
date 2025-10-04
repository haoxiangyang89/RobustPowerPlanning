using Distributed;
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
Γd = 1;
Γh = 3;
Γ = Dict("d" => Γd, "h" => Γh);

# create the first-stage model
mp = createFirst(fData,uData,hData,bData,T,vmaxT,vminT,θDmaxT,θDminT,xLimit,no_threads);
@constraint(mp,[k in fData.brList, t in 2:T],mp[:x][k,t] == 1);

# append the all-zero nominal scenario
uDict = Dict();
uDict["u_dp"] = Dict();
uDict["u_dm"] = Dict();
uDict["u_hp"] = Dict();
uDict["u_hm"] = Dict();
for t in 2:T
    for m in eachindex(groupDict[ci][2])
        uDict["u_dp"][m,t] = 0;
        uDict["u_dm"][m,t] = 0;
    end
    for i in hData.hList
        uDict["u_hp"][i,t] = 0;
        uDict["u_hm"][i,t] = 0;
    end
end
solveInfo = Dict("Γ" => [Γd, Γh], "iter" => []);
uList = [];
keepIter = true;

while keepIter
    # append the scenario to the first stage
    tFirst_start = time();
    push!(uList, uDict);
    mp = appendScen(mp,fData,uData,hData,T,groupDict[ci],vmaxT,vminT,θDmaxT,θDminT,expansion_factor,[uDict]);

    # solve the first stage solution to obtain the initial solution
    objhat,sphat,sqhat,xhat,yhat,zhat = solve_first(mp,fData,hData);
    tFirst_end = time();
    tFirst_elapsed = tFirst_end - tFirst_start;

    # feed the first-stage solution to the second stage problem
    # subp = createSecond(fData,uData,hData,T,groupDict[ci],Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uDict);
    tSecond_start = time();
    subp_dual = createSecond_dual(fData,uData,hData,T,groupDict[ci],Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,no_threads);

    # solve the second-stage problem and obtain the worst-case scenario
    uDict, subp_obj = solve_second(subp_dual,hData);
    tSecond_end = time();
    tSecond_elapsed = tSecond_end - tSecond_start;

    # tell when to stop: if the scenario has been generated
    if uDict in uList
        keepIter = false;
    end

    push!(solveInfo["iter"], [objhat, tFirst_elapsed, subp_obj, tSecond_elapsed]);
end

# obtain the optimal value and optimal solution
objhat,sphat,sqhat,xhat,yhat,zhat = solve_first(mp,fData,hData);
results_no_switching = [objhat,sphat,sqhat,xhat,yhat,zhat];
