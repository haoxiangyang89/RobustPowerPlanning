# create the dual formulation of the inner min of the min-max-min problem
function createSecond_dual(fData,uData,hData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,no_threads = 1)
    # first-stage model without any scenarios
    θu = Dict();
    for t in 1:T
        θu[t] = Dict();
        for k in fData.brList
            θu[t][k] = max(abs(θDmaxT[t][k]),abs(θDminT[t][k]));
        end
    end
    # scaling issue exists for this primal problem.
    # "NumericFocus" => 0, "BarConvTol" => 1e-6, "MIPGap" => 1e-6, "BarQCPConvTol" => 1e-6, "OptimalityTol" => 1e-6, "IntFeasTol" => 1e-6, "FeasibilityTol" => 1e-6,
    sprob = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 1, "Threads" => no_threads,"TimeLimit" => 600));

    # obtain the pairs that are connected
    connectPair = [];
    connectDict = Dict();
    branchDict1 = Dict();
    branchDict2 = Dict();
    for i in fData.IDList
        connectDict[i] = [];
        branchDict1[i] = [];
        branchDict2[i] = [];
    end
    for k in fData.brList
        push!(branchDict1[k[1]],k);
        push!(branchDict2[k[2]],k);
        if !((k[1],k[2]) in connectPair)
            push!(connectPair,(k[1],k[2]));
            push!(connectDict[k[1]],k[2]);
        end
    end
    kpDict = Dict();
    for k in fData.brList
        if (k[1],k[2]) in keys(kpDict)
            push!(kpDict[(k[1],k[2])],k);
        else
            kpDict[(k[1],k[2])] = [k];
        end
    end

    # add the strengthened bounds on cs/ss
    csmax = Dict();
    csmin = Dict();
    ssmax = Dict();
    ssmin = Dict();
    csConst = Dict();
    ssConst = Dict();
    for t in 2:T
        for k in fData.brList
            if (θDmaxT[t][k] >= 0)&(θDminT[t][k] >= 0)
                csmax[k,t] = cos(θDminT[t][k]);
                csmin[k,t] = cos(θDmaxT[t][k]);
            elseif (θDmaxT[t][k] < 0)&(θDminT[t][k] < 0)
                csmax[k,t] = cos(θDmaxT[t][k]);
                csmin[k,t] = cos(θDminT[t][k]);
            else
                csmax[k,t] = 1;
                csmin[k,t] = min(cos(θDmaxT[t][k]),cos(θDminT[t][k]));
            end
        end
        for k in fData.brList
            ssmax[k,t] = sin(θDmaxT[t][k]);
            ssmin[k,t] = sin(θDminT[t][k]);
            csConst[k,t] = (cos(θDmaxT[t][k]) - cos(θDminT[t][k]))/(θDmaxT[t][k] - θDminT[t][k]);
            ssConst[k,t] = (sin(θDmaxT[t][k]) - sin(θDminT[t][k]))/(θDmaxT[t][k] - θDminT[t][k]);
        end
    end

    M = Dict();
    for k in fData.brList
        if fData.rateA[k] < Inf
            M[k] = fData.rateA[k]*2;
        else
            M[k] = 100;
        end
    end

    vδ = Dict();
    θϕ = Dict();
    θδ = Dict();
    for t in 2:T
        for i in fData.IDList
            vδ[i,t] = vmaxT[t][i]+vminT[t][i];
        end
        for k in fData.brList
            θϕ[k,t] = (θDmaxT[t][k] + θDminT[t][k])/2;
            θδ[k,t] = (θDmaxT[t][k] - θDminT[t][k])/2;
        end
    end

    groupList = 1:length(groupDict[2]);
    group_rev = Dict();
    for i in fData.IDList
        for m in groupList
            if i in groupDict[2][m]
                group_rev[i] = m;
            end
        end
    end

    # create dual variables
    # dual variables for voltage constraints
    @variable(sprob,λv[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob,λvu[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob,λvl[i in fData.IDList, t in 2:T] <= 0);

    # dual variables for θ difference bounds
    @variable(sprob,λθ1[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λθ2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λθRef[t in 2:T]);

    # dual variables for generations
    @variable(sprob,λspu[i in fData.genIDList, t in 1:T] >= 0);
    @variable(sprob,λspl[i in fData.genIDList, t in 1:T] <= 0);
    @variable(sprob,λsqu[i in fData.genIDList, t in 1:T] >= 0);
    @variable(sprob,λsql[i in fData.genIDList, t in 1:T] <= 0);
    @variable(sprob,λspIni[i in fData.genIDList]);
    @variable(sprob,λsqIni[i in fData.genIDList]);
    @variable(sprob,λrampu[i in fData.genIDList, t in 2:T] >= 0);
    @variable(sprob,λrampl[i in fData.genIDList, t in 2:T] <= 0);

    # dual variables for cosine relaxation
    @variable(sprob,λcs1[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λcs2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λcs3[k in fData.brList, t in 2:T] <= 0);
    # dual variables for sine relaxation
    @variable(sprob,λss1[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λss2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λss4[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λss5[k in fData.brList, t in 2:T] <= 0);
    # dual variables for wc/ws/cs/ss/vv/l equalities of reverse flows
    @variable(sprob,λwce[k in fData.brList, t in 2:T]);
    @variable(sprob,λwse[k in fData.brList, t in 2:T]);
    @variable(sprob,λcse[k in fData.brList, t in 2:T]);
    @variable(sprob,λsse[k in fData.brList, t in 2:T]);
    @variable(sprob,λvve[k in fData.brList, t in 2:T]);

    # dual variables for power flow equations
    @variable(sprob,λptrans1[fData.brList, t in 2:T] >= 0);
    @variable(sprob,λptrans2[fData.brList, t in 2:T] <= 0);
    @variable(sprob,λqtrans1[fData.brList, t in 2:T] >= 0);
    @variable(sprob,λqtrans2[fData.brList, t in 2:T] <= 0);

    # dual variables for the SOC inequalities
    @variable(sprob,μ1[k in fData.brList,j in 1:2, t in 2:T; !(fData.rateA[k]==Inf)]);
    @variable(sprob,μ2[k in fData.brList,j in 1:3, t in 2:T]);
    @variable(sprob,μ3[k in fData.brList,j in 1:2, t in 2:T]);
    @variable(sprob,μ4[i in fData.IDList,j in 1:2, t in 2:T]);
    @variable(sprob,μ5[k in fData.brList,j in 1:4, t in 2:T]);
    @variable(sprob,ν1[k in fData.brList, t in 2:T; !(fData.rateA[k]==Inf)]>= 0);
    @variable(sprob,ν2[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,ν3[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,ν4[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob,ν5[k in fData.brList, t in 2:T] >= 0);

    # dual variables for vv McCormick
    @variable(sprob,λvv1[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λvv2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λvv3[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λvv4[k in fData.brList, t in 2:T] >= 0);
    # dual variables for wc McCormick
    @variable(sprob,λwc1[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λwc2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λwc3[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λwc4[k in fData.brList, t in 2:T] >= 0);
    # dual variables for ws McCormick
    @variable(sprob,λws1[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λws2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λws3[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λws4[k in fData.brList, t in 2:T] >= 0);

    # dual variables for tangent constraints
    @variable(sprob,λtangent1[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λtangent2[k in fData.brList, t in 2:T] <= 0);
    # dual variables for LNC constraints
    @variable(sprob,λlnc1[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λlnc2[k in fData.brList, t in 2:T] <= 0);

    # dual variables for battery dynamics
    @variable(sprob,λbat_lim[i in fData.IDList, t in 1:T] >= 0);
    @variable(sprob,λbat_trans[i in fData.IDList, t in 2:T]);
    @variable(sprob,μbat_out[i in fData.IDList, j in 1:2, t in 2:T]);
    @variable(sprob,νbat_out[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob,λbat_eff[i in fData.IDList, l in 1:length(bData.ηα[i]), t in 2:T] >= 0);

    # dual variables for power balance
    @variable(sprob, -fData.cz <= λpi[i in fData.IDList, t in 2:T] <= fData.cz);
    @variable(sprob, -fData.cz <= λqi[i in fData.IDList, t in 2:T] <= fData.cz);

    # dual variables for demand variables
    @variable(sprob, λdp[i in fData.IDList, t in 2:T]);
    @variable(sprob, λdq[i in fData.IDList, t in 2:T]);
    @variable(sprob, u_dp[m in groupList, t in 2:T], Bin);
    @variable(sprob, u_dm[m in groupList, t in 2:T], Bin);
    @variable(sprob, rdpp[i in fData.IDList, t in 2:T]);
    @variable(sprob, rdpm[i in fData.IDList, t in 2:T]);
    @variable(sprob, rdqp[i in fData.IDList, t in 2:T]);
    @variable(sprob, rdqm[i in fData.IDList, t in 2:T]);

    # dual variables for renewable generation variables
    @variable(sprob, λhp[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzp1[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzp2[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzp3[i in fData.IDList, t in 2:T] <= 0);
    @variable(sprob, λhzm1[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzm2[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzm3[i in fData.IDList, t in 2:T] <= 0);
    @variable(sprob, u_hp[i in hData.hList, t in 2:T], Bin);
    @variable(sprob, u_hm[i in hData.hList, t in 2:T], Bin);
    @variable(sprob, rhpp[i in hData.hList, t in 2:T]);
    @variable(sprob, rhpm[i in hData.hList, t in 2:T]);
    @variable(sprob, rhppz1[i in hData.hList, t in 2:T]);
    @variable(sprob, rhppz2[i in hData.hList, t in 2:T]);
    @variable(sprob, rhpmz1[i in hData.hList, t in 2:T]);
    @variable(sprob, rhpmz2[i in hData.hList, t in 2:T]);

    # dual variables for non-anticipitavity constraints
    @variable(sprob, λxIni[k in fData.brList, t in 2:T]);
    @variable(sprob, λyIni[i in fData.IDList]);
    @variable(sprob, λzIni[i in hData.hList]);

    # create dual constraints
    # constraints on the extreme points
    @constraint(sprob, uncertain_budget_d, sum(u_dp[m,t] + u_dm[m,t] for m in groupList for t in 2:T) <= Γ["d"]);
    @constraint(sprob, one_extreme_pt_d[m in groupList, t in 2:T], u_dp[m,t] + u_dm[m,t] <= 1);
    @constraint(sprob, uncertain_budget_h, sum(u_hp[i,t] + u_hm[i,t] for i in hData.hList for t in 2:T) <= Γ["h"]);
    @constraint(sprob, one_extreme_pt_h[i in hData.hList, t in 2:T], u_hp[i,t] + u_hm[i,t] <= 1);

    # linearizing the bilinear terms
    @constraint(sprob, rdppCons1[i in fData.IDList, t in 2:T], rdpp[i,t] <= fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdppCons2[i in fData.IDList, t in 2:T], rdpp[i,t] >= -fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdppCons3[i in fData.IDList, t in 2:T], rdpp[i,t] <= λdp[i,t] + 5 * (1 - u_dp[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdppCons4[i in fData.IDList, t in 2:T], rdpp[i,t] >= λdp[i,t] - 5 * (1 - u_dp[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rdpmCons1[i in fData.IDList, t in 2:T], rdpm[i,t] <= fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdpmCons2[i in fData.IDList, t in 2:T], rdpm[i,t] >= -fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdpmCons3[i in fData.IDList, t in 2:T], rdpm[i,t] <= λdp[i,t] + 5 * (1 - u_dm[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdpmCons4[i in fData.IDList, t in 2:T], rdpm[i,t] >= λdp[i,t] - 5 * (1 - u_dm[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rdqpCons1[i in fData.IDList, t in 2:T], rdqp[i,t] <= fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdqpCons2[i in fData.IDList, t in 2:T], rdqp[i,t] >= -fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdqpCons3[i in fData.IDList, t in 2:T], rdqp[i,t] <= λdq[i,t] + 5 * (1 - u_dp[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdqpCons4[i in fData.IDList, t in 2:T], rdqp[i,t] >= λdq[i,t] - 5 * (1 - u_dp[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rdqmCons1[i in fData.IDList, t in 2:T], rdqm[i,t] <= fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdqmCons2[i in fData.IDList, t in 2:T], rdqm[i,t] >= -fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdqmCons3[i in fData.IDList, t in 2:T], rdqm[i,t] <= λdq[i,t] + 5 * (1 - u_dm[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdqmCons4[i in fData.IDList, t in 2:T], rdqm[i,t] >= λdq[i,t] - 5 * (1 - u_dm[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rhppCons1[i in hData.hList, t in 2:T], rhpp[i,t] <= fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppCons2[i in hData.hList, t in 2:T], rhpp[i,t] >= -fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppCons3[i in hData.hList, t in 2:T], rhpp[i,t] <= λhp[i,t] + 5 * (1 - u_hp[i,t]) * fData.cz);
    @constraint(sprob, rhppCons4[i in hData.hList, t in 2:T], rhpp[i,t] >= λhp[i,t] - 5 * (1 - u_hp[i,t]) * fData.cz);

    @constraint(sprob, rhpmCons1[i in hData.hList, t in 2:T], rhpm[i,t] <= fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmCons2[i in hData.hList, t in 2:T], rhpm[i,t] >= -fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmCons3[i in hData.hList, t in 2:T], rhpm[i,t] <= λhp[i,t] + 5 * (1 - u_hm[i,t]) * fData.cz);
    @constraint(sprob, rhpmCons4[i in hData.hList, t in 2:T], rhpm[i,t] >= λhp[i,t] - 5 * (1 - u_hm[i,t]) * fData.cz);

    @constraint(sprob, rhppzCons11[i in hData.hList, t in 2:T], rhppz1[i,t] <= fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons12[i in hData.hList, t in 2:T], rhppz1[i,t] >= -fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons13[i in hData.hList, t in 2:T], rhppz1[i,t] <= λhzp1[i,t] + 5 * (1 - u_hp[i,t]) * fData.cz);
    @constraint(sprob, rhppzCons14[i in hData.hList, t in 2:T], rhppz1[i,t] >= λhzp1[i,t] - 5 * (1 - u_hp[i,t]) * fData.cz);

    @constraint(sprob, rhpmzCons11[i in hData.hList, t in 2:T], rhpmz1[i,t] <= fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons12[i in hData.hList, t in 2:T], rhpmz1[i,t] >= -fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons13[i in hData.hList, t in 2:T], rhpmz1[i,t] <= λhzm1[i,t] + 5 * (1 - u_hm[i,t]) * fData.cz);
    @constraint(sprob, rhpmzCons14[i in hData.hList, t in 2:T], rhpmz1[i,t] >= λhzm1[i,t] - 5 * (1 - u_hm[i,t]) * fData.cz);

    @constraint(sprob, rhppzCons21[i in hData.hList, t in 2:T], rhppz2[i,t] <= fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons22[i in hData.hList, t in 2:T], rhppz2[i,t] >= -fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons23[i in hData.hList, t in 2:T], rhppz2[i,t] <= λhzp3[i,t] + 5 * (1 - u_hp[i,t]) * fData.cz);
    @constraint(sprob, rhppzCons24[i in hData.hList, t in 2:T], rhppz2[i,t] >= λhzp3[i,t] - 5 * (1 - u_hp[i,t]) * fData.cz);

    @constraint(sprob, rhpmzCons21[i in hData.hList, t in 2:T], rhpmz2[i,t] <= fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons22[i in hData.hList, t in 2:T], rhpmz2[i,t] >= -fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons23[i in hData.hList, t in 2:T], rhpmz2[i,t] <= λhzm3[i,t] + 5 * (1 - u_hm[i,t]) * fData.cz);
    @constraint(sprob, rhpmzCons24[i in hData.hList, t in 2:T], rhpmz2[i,t] >= λhzm3[i,t] - 5 * (1 - u_hm[i,t]) * fData.cz);

    # set up the SOC constraints
    for t in 2:T
        for k in fData.brList
            if fData.rateA[k]<Inf
                μ1List = [μ1[k,1,t],μ1[k,2,t]];
                @constraint(sprob,[ν1[k,t]; μ1List] in SecondOrderCone());
            end
            μ2List = [μ2[k,1,t],μ2[k,2,t],μ2[k,3,t]];
            @constraint(sprob,[ν2[k,t]; μ2List] in SecondOrderCone());
            μ3List = [μ3[k,1,t],μ3[k,2,t]];
            @constraint(sprob,[ν3[k,t]; μ3List] in SecondOrderCone());
            μ5List = [μ5[k,1,t],μ5[k,2,t],μ5[k,3,t],μ5[k,4,t]];
            @constraint(sprob,[ν5[k,t]; μ5List] in SecondOrderCone());
        end
        for i in fData.IDList
            μ4List = [μ4[i,1,t],μ4[i,2,t]];
            @constraint(sprob,[ν4[i,t]; μ4List] in SecondOrderCone());
            μbatList = [μbat_out[i,1,t],μbat_out[i,2,t]];
            @constraint(sprob,[νbat_out[i,t]; μbatList] in SecondOrderCone());
        end
    end

    # set up the constraints for OPF variables
    @constraint(sprob, pConstr1[k in fData.brList, t in 2:T; fData.rateA[k]<Inf], λptrans1[k,t] + λptrans2[k,t] - μ1[k,1,t] + λpi[k[1],t] == 0);
    @constraint(sprob, pConstr2[k in fData.brList, t in 2:T; fData.rateA[k]==Inf], λptrans1[k,t] + λptrans2[k,t] + λpi[k[1],t] == 0);
    @constraint(sprob, qConstr1[k in fData.brList, t in 2:T; fData.rateA[k]<Inf], λqtrans1[k,t] + λqtrans2[k,t] - μ1[k,2,t] + λqi[k[1],t] == 0);
    @constraint(sprob, qConstr2[k in fData.brList, t in 2:T; fData.rateA[k]==Inf], λqtrans1[k,t] + λqtrans2[k,t] + λqi[k[1],t] == 0);
    @constraint(sprob, vConstr[i in fData.IDList, t in 2:T], -μ4[i,1,t] - (vmaxT[t][i] + vminT[t][i]) * λv[i,t] + λvu[i,t] + λvl[i,t] +
                sum(-vminT[t][k[2]]*λvv1[k,t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[2]]*λvv2[k,t]/(fData.τ1[k]*fData.τ2[k]) -
                    vmaxT[t][k[2]]*λvv3[k,t]/(fData.τ1[k]*fData.τ2[k]) - vminT[t][k[2]]*λvv4[k,t]/(fData.τ1[k]*fData.τ2[k]) for k in branchDict1[i]) +
                sum(-vminT[t][k[1]]*λvv1[k,t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*λvv2[k,t]/(fData.τ1[k]*fData.τ2[k]) -
                    vminT[t][k[1]]*λvv3[k,t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*λvv4[k,t]/(fData.τ1[k]*fData.τ2[k]) for k in branchDict2[i]) == 0);
    @constraint(sprob, vhatConstr[i in fData.IDList, t in 2:T], -μ4[i,2,t] - ν4[i,t] + λv[i,t] + λpi[i,t]*fData.gs[i] - λqi[i,t]*fData.bs[i]
                - sum(fData.g[k]*(λptrans1[k,t] + λptrans2[k,t])/(fData.τ1[k]^2) - 
                    (fData.b[k] + fData.bc[k]/2)*(λqtrans1[k,t] + λqtrans2[k,t])/(fData.τ1[k]^2) for k in branchDict1[i])
                - sum(μ2[k,2,t]/(fData.τ1[k]^2*sqrt(2)) + ν2[k,t]/(fData.τ1[k]^2*sqrt(2)) + 
                    μ5[k,3,t]/(fData.τ1[k]^2*sqrt(2)) + ν5[k,t]/(fData.τ1[k]^2*sqrt(2)) +
                    cos(θδ[k,t])*vδ[k[2],t]/fData.τ2[k]*(vmaxT[t][k[2]]/fData.τ2[k]*λlnc1[k,t] + 
                    vminT[t][k[2]]/fData.τ2[k]*λlnc2[k,t])/(fData.τ1[k]^2) for k in branchDict1[i])
                - sum(μ2[k,3,t]/(fData.τ2[k]^2*sqrt(2)) + ν2[k,t]/(fData.τ2[k]^2*sqrt(2)) + 
                    μ5[k,4,t]/(fData.τ2[k]^2*sqrt(2)) + ν5[k,t]/(fData.τ2[k]^2*sqrt(2)) +
                    cos(θδ[k,t])*vδ[k[1],t]/fData.τ1[k]*(vmaxT[t][k[1]]/fData.τ1[k]*λlnc1[k,t] + 
                    vminT[t][k[1]]/fData.τ1[k]*λlnc2[k,t])/(fData.τ2[k]^2) for k in branchDict2[i]) == 0);
    @constraint(sprob, vvConstr[k in fData.brList, t in 2:T],
                λvv1[k,t]+λvv2[k,t]+λvv3[k,t]+λvv4[k,t] + λvve[k,t] - λvve[(k[2],k[1],k[3]),t] - μ2[k,1,t]
                - csmin[k,t]*λwc1[k,t] - csmax[k,t]*λwc2[k,t] - csmax[k,t]*λwc3[k,t] - csmin[k,t]*λwc4[k,t]
                - ssmin[k,t]*λws1[k,t] - ssmax[k,t]*λws2[k,t] - ssmax[k,t]*λws3[k,t] - ssmin[k,t]*λws4[k,t] == 0);
    @constraint(sprob, wcConstr[k in fData.brList, t in 2:T],
                fData.g[k]*(λptrans1[k,t] + λptrans2[k,t]) - fData.b[k]*(λqtrans1[k,t] + λqtrans2[k,t]) - μ5[k,1,t] +
                λwc1[k,t] + λwc2[k,t] + λwc3[k,t] + λwc4[k,t] + λwce[k,t] - λwce[(k[2],k[1],k[3]),t] -
                λtangent1[k,t]*tan(θDmaxT[t][k]) - λtangent2[k,t]*tan(θDminT[t][k]) +
                vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*cos(θϕ[k,t])*(λlnc1[k,t] + λlnc2[k,t]) == 0);
    @constraint(sprob, wsConstr[k in fData.brList, t in 2:T],
                fData.b[k]*(λptrans1[k,t] + λptrans2[k,t]) + fData.g[k]*(λqtrans1[k,t] + λqtrans2[k,t]) + 
                λtangent1[k,t] + λtangent2[k,t] - μ5[k,2,t] +
                λws1[k,t] + λws2[k,t] + λws3[k,t] + λws4[k,t] + λwse[k,t] + λwse[(k[2],k[1],k[3]),t] +
                vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*sin(θϕ[k,t])*(λlnc1[k,t] + λlnc2[k,t]) == 0);
    @constraint(sprob, csConstr[k in fData.brList, t in 2:T],
                λcs1[k,t] + λcs2[k,t] + λcs3[k,t] - μ3[k,2,t] + ν3[k,t] + λcse[k,t] - λcse[(k[2],k[1],k[3]),t]
                - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc1[k,t] 
                - vmaxT[t][k[1]]*vmaxT[t][k[2]]*λwc2[k,t]/(fData.τ1[k]*fData.τ2[k])
                - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc3[k,t] 
                - vmaxT[t][k[1]]*vmaxT[t][k[2]]*λwc4[k,t]/(fData.τ1[k]*fData.τ2[k]) == 0);
    @constraint(sprob, ssConstr[k in fData.brList, t in 2:T],
                λss1[k,t] + λss2[k,t] + λss4[k,t] + λss5[k,t] + λsse[k,t] + λsse[(k[2],k[1],k[3]),t]
                - vminT[t][k[1]]*vminT[t][k[2]]*λws1[k,t]/(fData.τ1[k]*fData.τ2[k]) 
                - vmaxT[t][k[1]]*vmaxT[t][k[2]]*λws2[k,t]/(fData.τ1[k]*fData.τ2[k])
                - vminT[t][k[1]]*vminT[t][k[2]]*λws3[k,t]/(fData.τ1[k]*fData.τ2[k]) 
                - vmaxT[t][k[1]]*vmaxT[t][k[2]]*λws4[k,t]/(fData.τ1[k]*fData.τ2[k]) == 0);
    refBusInd = Dict();
    for i in fData.IDList
        if fData.bType[i] == 3
            refBusInd[i] = 1;
        else
            refBusInd[i] = 0;
        end
    end
    @constraint(sprob, θConstr[i in fData.IDList, t in 2:T],
                λθRef[t]*refBusInd[i] + sum(λθ1[k,t] + λθ2[k,t] - csConst[k,t]*λcs3[k,t] 
                - cos(θu[t][k]/2)*λss4[k,t] - cos(θu[t][k]/2)*λss5[k,t]
                - sqrt(1-cos(θu[t][k]))/θu[t][k]*μ3[k,1,t] for k in branchDict1[i])
                + sum(-λθ1[k,t] - λθ2[k,t] + csConst[k,t]*λcs3[k,t] 
                + cos(θu[t][k]/2)*λss4[k,t] + cos(θu[t][k]/2)*λss5[k,t]
                + sqrt(1-cos(θu[t][k]))/θu[t][k]*μ3[k,1,t] for k in branchDict2[i]) == 0);

    # set up constraints for battery variables
    @constraint(sprob, ep_varConstr[i in fData.IDList, t in 2:T], -μbat_out[i,1,t] + sum(λbat_eff[i,l,t] for l in 1:length(bData.ηα[i])) - λpi[i,t] == 0);
    @constraint(sprob, eq_varConstr[i in fData.IDList, t in 2:T], -μbat_out[i,2,t] - λqi[i,t] == 0);
    @constraint(sprob, f_varConstr[i in fData.IDList, t in 2:T], fData.Δt * λbat_trans[i,t] - sum(bData.ηα[i][l] * λbat_eff[i,l,t] for l in 1:length(bData.ηα[i])) == 0);
    @constraint(sprob, I_varConstr[i in fData.IDList, t in 2:(T-1)], λbat_lim[i,t] + λbat_trans[i,t] - λbat_trans[i,t+1] >= 0);
    @constraint(sprob, I_varConstrT[i in fData.IDList], λbat_lim[i,T] + λbat_trans[i,T] >= 0);
    @constraint(sprob, I_varConstr1[i in fData.IDList], λbat_lim[i,1] - λbat_trans[i,2] >= 0);

    # set up constraints for demand/renewable variables
    @constraint(sprob, dp_varConstr[i in fData.IDList, t in 2:T], λpi[i,t] + λdp[i,t] == 0);
    @constraint(sprob, dq_varConstr[i in fData.IDList, t in 2:T], λqi[i,t] + λdq[i,t] == 0);
    @constraint(sprob, h_varConstr[i in fData.IDList, t in 2:T], -λpi[i,t] + λhp[i,t] >= 0);
    @constraint(sprob, uzhpConstr[i in hData.hList, t in 2:T], λhzp1[i,t] + λhzp2[i,t] + λhzp3[i,t] - uData[i].RESP0[t] * expansion_factor * uData[i].RESPmax * λhp[i,t] >= 0);
    @constraint(sprob, uzhmConstr[i in hData.hList, t in 2:T], λhzm1[i,t] + λhzm2[i,t] + λhzm3[i,t] + uData[i].RESP0[t] * expansion_factor * uData[i].RESPmin * λhp[i,t] >= 0);

    # set up constraints for x/y/z
    @constraint(sprob, xConstr1[k in fData.brList, t in 2:T; fData.rateA[k] < Inf], M[k] * (λptrans1[k,t] - λptrans2[k,t] + λqtrans1[k,t] - λqtrans2[k,t])
                 - fData.rateA[k] * ν1[k,t] + (λθ1[k,t] - λθ2[k,t])*2*pi + λxIni[k,t] == 0);
    @constraint(sprob, xConstr2[k in fData.brList, t in 2:T; fData.rateA[k] == Inf], M[k] * (λptrans1[k,t] - λptrans2[k,t] + λqtrans1[k,t] - λqtrans2[k,t])
                 + (λθ1[k,t] - λθ2[k,t])*2*pi + λxIni[k,t] == 0);
    @constraint(sprob, yConstr[i in fData.IDList], -sum(bData.cap[i]*λbat_lim[i,t] + νbat_out[i,t]*bData.uCap[i] for t in 2:T) - bData.cap[i]*λbat_lim[i,1] + λyIni[i] == 0);
    @constraint(sprob, zConstr[i in hData.hList], -sum(λhzp2[i,t] + λhzp3[i,t] + λhzm2[i,t] + λhzm3[i,t] + λhp[i,t] * uData[i].RESP0[t] * expansion_factor for t in 2:T)
                 + λzIni[i] == 0);

    # set up constraints for sp/sq
    @constraint(sprob, spConstr[i in fData.genIDList, t in 2:(T-1)], λspu[i,t] + λspl[i,t] + λrampu[i,t] - λrampu[i,t+1] + λrampl[i,t] - λrampl[i,t+1] - λpi[fData.Loc[i],t] ==
                -fData.cp[i].params[2]);
    @constraint(sprob, sqConstr[i in fData.genIDList, t in 2:T], λsqu[i,t] + λsql[i,t] - λqi[fData.Loc[i],t] == 0);
    @constraint(sprob, spConstr1[i in fData.genIDList], λspu[i,1] + λspl[i,1] + λspIni[i] - λrampu[i,2] - λrampl[i,2] == 0);
    @constraint(sprob, sqConstr1[i in fData.genIDList], λsqu[i,1] + λsql[i,1] + λsqIni[i] == 0);
    @constraint(sprob, spConstrT[i in fData.genIDList], λspu[i,T] + λspl[i,T] + λrampu[i,T] + λrampl[i,T] - λpi[fData.Loc[i],T] == -fData.cp[i].params[2]);

    # set up the objective function
    @objective(sprob, Max, -sum(sum(-μ4[i,2,t]/4 + ν4[i,t]/4 + λvu[i,t]*vmaxT[t][i] + λvl[i,t]*vminT[t][i] - vmaxT[t][i]*vminT[t][i]*λv[i,t] 
        + sum(bData.ηβ[i][l] * λbat_eff[i,l,t] for l in 1:length(bData.ηα[i])) + λdp[i,t]*uData[i].DP0[t] + λdq[i,t]*uData[i].DQ0[t]
        + (uData[i].DPmax[t] - uData[i].DP0[t])*rdpp[i,t] + (uData[i].DPmin[t] - uData[i].DP0[t])*rdpm[i,t] 
        + (uData[i].DQmax[t] - uData[i].DQ0[t])*rdqp[i,t] + (uData[i].DQmin[t] - uData[i].DQ0[t])*rdqm[i,t] 
            for t in 2:T) + yhat[i]*λyIni[i] for i in fData.IDList) - 
        sum(sum(fData.Pmax[i] * λspu[i,t] + fData.Pmin[i] * λspl[i,t] + fData.Qmax[i] * λsqu[i,t] + fData.Qmin[i] * λsql[i,t] 
        + λrampu[i,t] * fData.RU[i] + λrampl[i,t] * fData.RD[i] for t in 2:T)
        + fData.Pmax[i] * λspu[i,1] + fData.Pmin[i] * λspl[i,1] + fData.Qmax[i] * λsqu[i,1] + fData.Qmin[i] * λsql[i,1] 
        + sphat[i] * λspIni[i] + sqhat[i] * λsqIni[i] for i in fData.genIDList) -
        sum(-3/4*μ3[k,2,t] + 5/4*ν3[k,t] + csmax[k,t]*λcs1[k,t] + csmin[k,t]*λcs2[k,t] + (cos(θDminT[t][k]) - csConst[k,t]*(fData.σ[k]+θDminT[t][k]))*λcs3[k,t]
        + (sin(θu[t][k]/2) - θu[t][k]/2*cos(θu[t][k]/2))*(λss4[k,t] - λss5[k,t]) + ssmax[k,t]*λss1[k,t] + ssmin[k,t]*λss2[k,t]
        + (θDmaxT[t][k] + fData.σ[k] + 2*pi)*λθ1[k,t] + (θDminT[t][k] + fData.σ[k] - 2*pi)*λθ2[k,t] - fData.σ[k]*(λss4[k,t] + λss5[k,t])
        + M[k]*(λptrans1[k,t] - λptrans2[k,t] + λqtrans1[k,t] - λqtrans2[k,t]) + xhat[k,t]*λxIni[k,t]
        + vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]] - vmaxT[t][k[1]]*vmaxT[t][k[2]])/(fData.τ1[k]*fData.τ2[k])*λlnc1[k,t]
        - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]] - vmaxT[t][k[1]]*vmaxT[t][k[2]])/(fData.τ1[k]*fData.τ2[k])*λlnc2[k,t]
        - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv1[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv2[k,t]
        - vminT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv3[k,t] - vmaxT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv4[k,t]
        - csmin[k,t]*vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc1[k,t] - csmax[k,t]*vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc2[k,t]
        - csmax[k,t]*vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc3[k,t] - csmin[k,t]*vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc4[k,t]
        - ssmin[k,t]*vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λws1[k,t] - ssmax[k,t]*vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λws2[k,t]
        - ssmax[k,t]*vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λws3[k,t] - ssmin[k,t]*vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λws4[k,t]
           for t in 2:T for k in fData.brList) - 
        sum(sum(λhp[i,t]*uData[i].RESP0[t] - λhzp3[i,t] - λhzm3[i,t] + rhpp[i,t]*uData[i].RESPmax*uData[i].RESP0[t] - rhpm[i,t]*uData[i].RESPmin*uData[i].RESP0[t]
        + rhppz1[i,t] + rhppz2[i,t] + rhpmz1[i,t] + rhpmz2[i,t] for t in 2:T) + zhat[i]*λzIni[i] for i in hData.hList)
            );

    return sprob;
end

# given an uncertainty realization, solve the subproblem to get the scenario cost
function testSecond_dual(fData,uData,hData,T,groupDict,uhat,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,no_threads = 1)
    # first-stage model without any scenarios
    θu = Dict();
    for t in 1:T
        θu[t] = Dict();
        for k in fData.brList
            θu[t][k] = max(abs(θDmaxT[t][k]),abs(θDminT[t][k]));
        end
    end
    # scaling issue exists for this primal problem.
    # "NumericFocus" => 0, "BarConvTol" => 1e-6, "MIPGap" => 1e-6, "BarQCPConvTol" => 1e-6, "OptimalityTol" => 1e-6, "IntFeasTol" => 1e-6, "FeasibilityTol" => 1e-6,
    # sprob = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 1, "Threads" => no_threads,"TimeLimit" => 300));
    sprob = Model(COPT.Optimizer);

    # obtain the pairs that are connected
    connectPair = [];
    connectDict = Dict();
    branchDict1 = Dict();
    branchDict2 = Dict();
    for i in fData.IDList
        connectDict[i] = [];
        branchDict1[i] = [];
        branchDict2[i] = [];
    end
    for k in fData.brList
        push!(branchDict1[k[1]],k);
        push!(branchDict2[k[2]],k);
        if !((k[1],k[2]) in connectPair)
            push!(connectPair,(k[1],k[2]));
            push!(connectDict[k[1]],k[2]);
        end
    end
    kpDict = Dict();
    for k in fData.brList
        if (k[1],k[2]) in keys(kpDict)
            push!(kpDict[(k[1],k[2])],k);
        else
            kpDict[(k[1],k[2])] = [k];
        end
    end

    # add the strengthened bounds on cs/ss
    csmax = Dict();
    csmin = Dict();
    ssmax = Dict();
    ssmin = Dict();
    csConst = Dict();
    ssConst = Dict();
    for t in 2:T
        for k in fData.brList
            if (θDmaxT[t][k] >= 0)&(θDminT[t][k] >= 0)
                csmax[k,t] = cos(θDminT[t][k]);
                csmin[k,t] = cos(θDmaxT[t][k]);
            elseif (θDmaxT[t][k] < 0)&(θDminT[t][k] < 0)
                csmax[k,t] = cos(θDmaxT[t][k]);
                csmin[k,t] = cos(θDminT[t][k]);
            else
                csmax[k,t] = 1;
                csmin[k,t] = min(cos(θDmaxT[t][k]),cos(θDminT[t][k]));
            end
        end
        for k in fData.brList
            ssmax[k,t] = sin(θDmaxT[t][k]);
            ssmin[k,t] = sin(θDminT[t][k]);
            csConst[k,t] = (cos(θDmaxT[t][k]) - cos(θDminT[t][k]))/(θDmaxT[t][k] - θDminT[t][k]);
            ssConst[k,t] = (sin(θDmaxT[t][k]) - sin(θDminT[t][k]))/(θDmaxT[t][k] - θDminT[t][k]);
        end
    end

    M = Dict();
    for k in fData.brList
        if fData.rateA[k] < Inf
            M[k] = fData.rateA[k]*2;
        else
            M[k] = 100;
        end
    end

    vδ = Dict();
    θϕ = Dict();
    θδ = Dict();
    for t in 2:T
        for i in fData.IDList
            vδ[i,t] = vmaxT[t][i]+vminT[t][i];
        end
        for k in fData.brList
            θϕ[k,t] = (θDmaxT[t][k] + θDminT[t][k])/2;
            θδ[k,t] = (θDmaxT[t][k] - θDminT[t][k])/2;
        end
    end

    groupList = 1:length(groupDict[2]);
    group_rev = Dict();
    for i in fData.IDList
        for m in groupList
            if i in groupDict[2][m]
                group_rev[i] = m;
            end
        end
    end

    # create dual variables
    # dual variables for voltage constraints
    @variable(sprob,λv[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob,λvu[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob,λvl[i in fData.IDList, t in 2:T] <= 0);

    # dual variables for θ difference bounds
    @variable(sprob,λθ1[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λθ2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λθRef[t in 2:T]);

    # dual variables for generations
    @variable(sprob,λspu[i in fData.genIDList, t in 1:T] >= 0);
    @variable(sprob,λspl[i in fData.genIDList, t in 1:T] <= 0);
    @variable(sprob,λsqu[i in fData.genIDList, t in 1:T] >= 0);
    @variable(sprob,λsql[i in fData.genIDList, t in 1:T] <= 0);
    @variable(sprob,λspIni[i in fData.genIDList]);
    @variable(sprob,λsqIni[i in fData.genIDList]);
    @variable(sprob,λrampu[i in fData.genIDList, t in 2:T] >= 0);
    @variable(sprob,λrampl[i in fData.genIDList, t in 2:T] <= 0);

    # dual variables for cosine relaxation
    @variable(sprob,λcs1[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λcs2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λcs3[k in fData.brList, t in 2:T] <= 0);
    # dual variables for sine relaxation
    @variable(sprob,λss1[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λss2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λss4[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λss5[k in fData.brList, t in 2:T] <= 0);
    # dual variables for wc/ws/cs/ss/vv/l equalities of reverse flows
    @variable(sprob,λwce[k in fData.brList, t in 2:T]);
    @variable(sprob,λwse[k in fData.brList, t in 2:T]);
    @variable(sprob,λcse[k in fData.brList, t in 2:T]);
    @variable(sprob,λsse[k in fData.brList, t in 2:T]);
    @variable(sprob,λvve[k in fData.brList, t in 2:T]);

    # dual variables for power flow equations
    @variable(sprob,λptrans1[fData.brList, t in 2:T] >= 0);
    @variable(sprob,λptrans2[fData.brList, t in 2:T] <= 0);
    @variable(sprob,λqtrans1[fData.brList, t in 2:T] >= 0);
    @variable(sprob,λqtrans2[fData.brList, t in 2:T] <= 0);

    # dual variables for the SOC inequalities
    @variable(sprob,μ1[k in fData.brList,j in 1:2, t in 2:T; !(fData.rateA[k]==Inf)]);
    @variable(sprob,μ2[k in fData.brList,j in 1:3, t in 2:T]);
    @variable(sprob,μ3[k in fData.brList,j in 1:2, t in 2:T]);
    @variable(sprob,μ4[i in fData.IDList,j in 1:2, t in 2:T]);
    @variable(sprob,μ5[k in fData.brList,j in 1:4, t in 2:T]);
    @variable(sprob,ν1[k in fData.brList, t in 2:T; !(fData.rateA[k]==Inf)]>= 0);
    @variable(sprob,ν2[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,ν3[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,ν4[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob,ν5[k in fData.brList, t in 2:T] >= 0);

    # dual variables for vv McCormick
    @variable(sprob,λvv1[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λvv2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λvv3[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λvv4[k in fData.brList, t in 2:T] >= 0);
    # dual variables for wc McCormick
    @variable(sprob,λwc1[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λwc2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λwc3[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λwc4[k in fData.brList, t in 2:T] >= 0);
    # dual variables for ws McCormick
    @variable(sprob,λws1[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λws2[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λws3[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λws4[k in fData.brList, t in 2:T] >= 0);

    # dual variables for tangent constraints
    @variable(sprob,λtangent1[k in fData.brList, t in 2:T] >= 0);
    @variable(sprob,λtangent2[k in fData.brList, t in 2:T] <= 0);
    # dual variables for LNC constraints
    @variable(sprob,λlnc1[k in fData.brList, t in 2:T] <= 0);
    @variable(sprob,λlnc2[k in fData.brList, t in 2:T] <= 0);

    # dual variables for battery dynamics
    @variable(sprob,λbat_lim[i in fData.IDList, t in 1:T] >= 0);
    @variable(sprob,λbat_trans[i in fData.IDList, t in 2:T]);
    @variable(sprob,μbat_out[i in fData.IDList, j in 1:2, t in 2:T]);
    @variable(sprob,νbat_out[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob,λbat_eff[i in fData.IDList, l in 1:length(bData.ηα[i]), t in 2:T] >= 0);

    # dual variables for power balance
    @variable(sprob, -fData.cz <= λpi[i in fData.IDList, t in 2:T] <= fData.cz);
    @variable(sprob, -fData.cz <= λqi[i in fData.IDList, t in 2:T] <= fData.cz);

    # dual variables for demand variables
    @variable(sprob, λdp[i in fData.IDList, t in 2:T]);
    @variable(sprob, λdq[i in fData.IDList, t in 2:T]);
    u_dp = Dict();
    u_dm = Dict();
    for m in groupList
        for t in 2:T
            u_dp[m,t] = uhat["u_dp"][m,t];
            u_dm[m,t] = uhat["u_dm"][m,t];
        end
    end
    @variable(sprob, rdpp[i in fData.IDList, t in 2:T]);
    @variable(sprob, rdpm[i in fData.IDList, t in 2:T]);
    @variable(sprob, rdqp[i in fData.IDList, t in 2:T]);
    @variable(sprob, rdqm[i in fData.IDList, t in 2:T]);

    # dual variables for renewable generation variables
    @variable(sprob, λhp[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzp1[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzp2[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzp3[i in fData.IDList, t in 2:T] <= 0);
    @variable(sprob, λhzm1[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzm2[i in fData.IDList, t in 2:T] >= 0);
    @variable(sprob, λhzm3[i in fData.IDList, t in 2:T] <= 0);
    u_hp = Dict();
    u_hm = Dict();
    for i in hData.hList
        for t in 2:T
            u_hp[i,t] = uhat["u_hp"][i,t];
            u_hm[i,t] = uhat["u_hm"][i,t];
        end
    end
    @variable(sprob, rhpp[i in hData.hList, t in 2:T]);
    @variable(sprob, rhpm[i in hData.hList, t in 2:T]);
    @variable(sprob, rhppz1[i in hData.hList, t in 2:T]);
    @variable(sprob, rhppz2[i in hData.hList, t in 2:T]);
    @variable(sprob, rhpmz1[i in hData.hList, t in 2:T]);
    @variable(sprob, rhpmz2[i in hData.hList, t in 2:T]);

    # dual variables for non-anticipitavity constraints
    @variable(sprob, λxIni[k in fData.brList, t in 2:T]);
    @variable(sprob, λyIni[i in fData.IDList]);
    @variable(sprob, λzIni[i in hData.hList]);

    # create dual constraints
    # linearizing the bilinear terms
    @constraint(sprob, rdppCons1[i in fData.IDList, t in 2:T], rdpp[i,t] <= fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdppCons2[i in fData.IDList, t in 2:T], rdpp[i,t] >= -fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdppCons3[i in fData.IDList, t in 2:T], rdpp[i,t] <= λdp[i,t] + 5 * (1 - u_dp[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdppCons4[i in fData.IDList, t in 2:T], rdpp[i,t] >= λdp[i,t] - 5 * (1 - u_dp[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rdpmCons1[i in fData.IDList, t in 2:T], rdpm[i,t] <= fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdpmCons2[i in fData.IDList, t in 2:T], rdpm[i,t] >= -fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdpmCons3[i in fData.IDList, t in 2:T], rdpm[i,t] <= λdp[i,t] + 5 * (1 - u_dm[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdpmCons4[i in fData.IDList, t in 2:T], rdpm[i,t] >= λdp[i,t] - 5 * (1 - u_dm[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rdqpCons1[i in fData.IDList, t in 2:T], rdqp[i,t] <= fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdqpCons2[i in fData.IDList, t in 2:T], rdqp[i,t] >= -fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdqpCons3[i in fData.IDList, t in 2:T], rdqp[i,t] <= λdq[i,t] + 5 * (1 - u_dp[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdqpCons4[i in fData.IDList, t in 2:T], rdqp[i,t] >= λdq[i,t] - 5 * (1 - u_dp[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rdqmCons1[i in fData.IDList, t in 2:T], rdqm[i,t] <= fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdqmCons2[i in fData.IDList, t in 2:T], rdqm[i,t] >= -fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdqmCons3[i in fData.IDList, t in 2:T], rdqm[i,t] <= λdq[i,t] + 5 * (1 - u_dm[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdqmCons4[i in fData.IDList, t in 2:T], rdqm[i,t] >= λdq[i,t] - 5 * (1 - u_dm[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rhppCons1[i in hData.hList, t in 2:T], rhpp[i,t] <= fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppCons2[i in hData.hList, t in 2:T], rhpp[i,t] >= -fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppCons3[i in hData.hList, t in 2:T], rhpp[i,t] <= λhp[i,t] + 5 * (1 - u_hp[i,t]) * fData.cz);
    @constraint(sprob, rhppCons4[i in hData.hList, t in 2:T], rhpp[i,t] >= λhp[i,t] - 5 * (1 - u_hp[i,t]) * fData.cz);

    @constraint(sprob, rhpmCons1[i in hData.hList, t in 2:T], rhpm[i,t] <= fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmCons2[i in hData.hList, t in 2:T], rhpm[i,t] >= -fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmCons3[i in hData.hList, t in 2:T], rhpm[i,t] <= λhp[i,t] + 5 * (1 - u_hm[i,t]) * fData.cz);
    @constraint(sprob, rhpmCons4[i in hData.hList, t in 2:T], rhpm[i,t] >= λhp[i,t] - 5 * (1 - u_hm[i,t]) * fData.cz);

    @constraint(sprob, rhppzCons11[i in hData.hList, t in 2:T], rhppz1[i,t] <= fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons12[i in hData.hList, t in 2:T], rhppz1[i,t] >= -fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons13[i in hData.hList, t in 2:T], rhppz1[i,t] <= λhzp1[i,t] + 5 * (1 - u_hp[i,t]) * fData.cz);
    @constraint(sprob, rhppzCons14[i in hData.hList, t in 2:T], rhppz1[i,t] >= λhzp1[i,t] - 5 * (1 - u_hp[i,t]) * fData.cz);

    @constraint(sprob, rhpmzCons11[i in hData.hList, t in 2:T], rhpmz1[i,t] <= fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons12[i in hData.hList, t in 2:T], rhpmz1[i,t] >= -fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons13[i in hData.hList, t in 2:T], rhpmz1[i,t] <= λhzm1[i,t] + 5 * (1 - u_hm[i,t]) * fData.cz);
    @constraint(sprob, rhpmzCons14[i in hData.hList, t in 2:T], rhpmz1[i,t] >= λhzm1[i,t] - 5 * (1 - u_hm[i,t]) * fData.cz);

    @constraint(sprob, rhppzCons21[i in hData.hList, t in 2:T], rhppz2[i,t] <= fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons22[i in hData.hList, t in 2:T], rhppz2[i,t] >= -fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons23[i in hData.hList, t in 2:T], rhppz2[i,t] <= λhzp3[i,t] + 5 * (1 - u_hp[i,t]) * fData.cz);
    @constraint(sprob, rhppzCons24[i in hData.hList, t in 2:T], rhppz2[i,t] >= λhzp3[i,t] - 5 * (1 - u_hp[i,t]) * fData.cz);

    @constraint(sprob, rhpmzCons21[i in hData.hList, t in 2:T], rhpmz2[i,t] <= fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons22[i in hData.hList, t in 2:T], rhpmz2[i,t] >= -fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons23[i in hData.hList, t in 2:T], rhpmz2[i,t] <= λhzm3[i,t] + 5 * (1 - u_hm[i,t]) * fData.cz);
    @constraint(sprob, rhpmzCons24[i in hData.hList, t in 2:T], rhpmz2[i,t] >= λhzm3[i,t] - 5 * (1 - u_hm[i,t]) * fData.cz);

    # set up the SOC constraints
    for t in 2:T
        for k in fData.brList
            if fData.rateA[k]<Inf
                μ1List = [μ1[k,1,t],μ1[k,2,t]];
                @constraint(sprob,[ν1[k,t]; μ1List] in SecondOrderCone());
            end
            μ2List = [μ2[k,1,t],μ2[k,2,t],μ2[k,3,t]];
            @constraint(sprob,[ν2[k,t]; μ2List] in SecondOrderCone());
            μ3List = [μ3[k,1,t],μ3[k,2,t]];
            @constraint(sprob,[ν3[k,t]; μ3List] in SecondOrderCone());
            μ5List = [μ5[k,1,t],μ5[k,2,t],μ5[k,3,t],μ5[k,4,t]];
            @constraint(sprob,[ν5[k,t]; μ5List] in SecondOrderCone());
        end
        for i in fData.IDList
            μ4List = [μ4[i,1,t],μ4[i,2,t]];
            @constraint(sprob,[ν4[i,t]; μ4List] in SecondOrderCone());
            μbatList = [μbat_out[i,1,t],μbat_out[i,2,t]];
            @constraint(sprob,[νbat_out[i,t]; μbatList] in SecondOrderCone());
        end
    end

    # set up the constraints for OPF variables
    @constraint(sprob, pConstr1[k in fData.brList, t in 2:T; fData.rateA[k]<Inf], λptrans1[k,t] + λptrans2[k,t] - μ1[k,1,t] + λpi[k[1],t] == 0);
    @constraint(sprob, pConstr2[k in fData.brList, t in 2:T; fData.rateA[k]==Inf], λptrans1[k,t] + λptrans2[k,t] + λpi[k[1],t] == 0);
    @constraint(sprob, qConstr1[k in fData.brList, t in 2:T; fData.rateA[k]<Inf], λqtrans1[k,t] + λqtrans2[k,t] - μ1[k,2,t] + λqi[k[1],t] == 0);
    @constraint(sprob, qConstr2[k in fData.brList, t in 2:T; fData.rateA[k]==Inf], λqtrans1[k,t] + λqtrans2[k,t] + λqi[k[1],t] == 0);
    @constraint(sprob, vConstr[i in fData.IDList, t in 2:T], -μ4[i,1,t] - (vmaxT[t][i] + vminT[t][i]) * λv[i,t] + λvu[i,t] + λvl[i,t] +
                sum(-vminT[t][k[2]]*λvv1[k,t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[2]]*λvv2[k,t]/(fData.τ1[k]*fData.τ2[k]) -
                    vmaxT[t][k[2]]*λvv3[k,t]/(fData.τ1[k]*fData.τ2[k]) - vminT[t][k[2]]*λvv4[k,t]/(fData.τ1[k]*fData.τ2[k]) for k in branchDict1[i]) +
                sum(-vminT[t][k[1]]*λvv1[k,t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*λvv2[k,t]/(fData.τ1[k]*fData.τ2[k]) -
                    vminT[t][k[1]]*λvv3[k,t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*λvv4[k,t]/(fData.τ1[k]*fData.τ2[k]) for k in branchDict2[i]) == 0);
    @constraint(sprob, vhatConstr[i in fData.IDList, t in 2:T], -μ4[i,2,t] - ν4[i,t] + λv[i,t] + λpi[i,t]*fData.gs[i] - λqi[i,t]*fData.bs[i]
                - sum(fData.g[k]*(λptrans1[k,t] + λptrans2[k,t])/(fData.τ1[k]^2) - 
                    (fData.b[k] + fData.bc[k]/2)*(λqtrans1[k,t] + λqtrans2[k,t])/(fData.τ1[k]^2) for k in branchDict1[i])
                - sum(μ2[k,2,t]/(fData.τ1[k]^2*sqrt(2)) + ν2[k,t]/(fData.τ1[k]^2*sqrt(2)) + 
                    μ5[k,3,t]/(fData.τ1[k]^2*sqrt(2)) + ν5[k,t]/(fData.τ1[k]^2*sqrt(2)) +
                    cos(θδ[k,t])*vδ[k[2],t]/fData.τ2[k]*(vmaxT[t][k[2]]/fData.τ2[k]*λlnc1[k,t] + 
                    vminT[t][k[2]]/fData.τ2[k]*λlnc2[k,t])/(fData.τ1[k]^2) for k in branchDict1[i])
                - sum(μ2[k,3,t]/(fData.τ2[k]^2*sqrt(2)) + ν2[k,t]/(fData.τ2[k]^2*sqrt(2)) + 
                    μ5[k,4,t]/(fData.τ2[k]^2*sqrt(2)) + ν5[k,t]/(fData.τ2[k]^2*sqrt(2)) +
                    cos(θδ[k,t])*vδ[k[1],t]/fData.τ1[k]*(vmaxT[t][k[1]]/fData.τ1[k]*λlnc1[k,t] + 
                    vminT[t][k[1]]/fData.τ1[k]*λlnc2[k,t])/(fData.τ2[k]^2) for k in branchDict2[i]) == 0);
    @constraint(sprob, vvConstr[k in fData.brList, t in 2:T],
                λvv1[k,t]+λvv2[k,t]+λvv3[k,t]+λvv4[k,t] + λvve[k,t] - λvve[(k[2],k[1],k[3]),t] - μ2[k,1,t]
                - csmin[k,t]*λwc1[k,t] - csmax[k,t]*λwc2[k,t] - csmax[k,t]*λwc3[k,t] - csmin[k,t]*λwc4[k,t]
                - ssmin[k,t]*λws1[k,t] - ssmax[k,t]*λws2[k,t] - ssmax[k,t]*λws3[k,t] - ssmin[k,t]*λws4[k,t] == 0);
    @constraint(sprob, wcConstr[k in fData.brList, t in 2:T],
                fData.g[k]*(λptrans1[k,t] + λptrans2[k,t]) - fData.b[k]*(λqtrans1[k,t] + λqtrans2[k,t]) - μ5[k,1,t] +
                λwc1[k,t] + λwc2[k,t] + λwc3[k,t] + λwc4[k,t] + λwce[k,t] - λwce[(k[2],k[1],k[3]),t] -
                λtangent1[k,t]*tan(θDmaxT[t][k]) - λtangent2[k,t]*tan(θDminT[t][k]) +
                vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*cos(θϕ[k,t])*(λlnc1[k,t] + λlnc2[k,t]) == 0);
    @constraint(sprob, wsConstr[k in fData.brList, t in 2:T],
                fData.b[k]*(λptrans1[k,t] + λptrans2[k,t]) + fData.g[k]*(λqtrans1[k,t] + λqtrans2[k,t]) + 
                λtangent1[k,t] + λtangent2[k,t] - μ5[k,2,t] +
                λws1[k,t] + λws2[k,t] + λws3[k,t] + λws4[k,t] + λwse[k,t] + λwse[(k[2],k[1],k[3]),t] +
                vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*sin(θϕ[k,t])*(λlnc1[k,t] + λlnc2[k,t]) == 0);
    @constraint(sprob, csConstr[k in fData.brList, t in 2:T],
                λcs1[k,t] + λcs2[k,t] + λcs3[k,t] - μ3[k,2,t] + ν3[k,t] + λcse[k,t] - λcse[(k[2],k[1],k[3]),t]
                - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc1[k,t] 
                - vmaxT[t][k[1]]*vmaxT[t][k[2]]*λwc2[k,t]/(fData.τ1[k]*fData.τ2[k])
                - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc3[k,t] 
                - vmaxT[t][k[1]]*vmaxT[t][k[2]]*λwc4[k,t]/(fData.τ1[k]*fData.τ2[k]) == 0);
    @constraint(sprob, ssConstr[k in fData.brList, t in 2:T],
                λss1[k,t] + λss2[k,t] + λss4[k,t] + λss5[k,t] + λsse[k,t] + λsse[(k[2],k[1],k[3]),t]
                - vminT[t][k[1]]*vminT[t][k[2]]*λws1[k,t]/(fData.τ1[k]*fData.τ2[k]) 
                - vmaxT[t][k[1]]*vmaxT[t][k[2]]*λws2[k,t]/(fData.τ1[k]*fData.τ2[k])
                - vminT[t][k[1]]*vminT[t][k[2]]*λws3[k,t]/(fData.τ1[k]*fData.τ2[k]) 
                - vmaxT[t][k[1]]*vmaxT[t][k[2]]*λws4[k,t]/(fData.τ1[k]*fData.τ2[k]) == 0);
    refBusInd = Dict();
    for i in fData.IDList
        if fData.bType[i] == 3
            refBusInd[i] = 1;
        else
            refBusInd[i] = 0;
        end
    end
    @constraint(sprob, θConstr[i in fData.IDList, t in 2:T],
                λθRef[t]*refBusInd[i] + sum(λθ1[k,t] + λθ2[k,t] - csConst[k,t]*λcs3[k,t] 
                - cos(θu[t][k]/2)*λss4[k,t] - cos(θu[t][k]/2)*λss5[k,t]
                - sqrt(1-cos(θu[t][k]))/θu[t][k]*μ3[k,1,t] for k in branchDict1[i])
                + sum(-λθ1[k,t] - λθ2[k,t] + csConst[k,t]*λcs3[k,t] 
                + cos(θu[t][k]/2)*λss4[k,t] + cos(θu[t][k]/2)*λss5[k,t]
                + sqrt(1-cos(θu[t][k]))/θu[t][k]*μ3[k,1,t] for k in branchDict2[i]) == 0);

    # set up constraints for battery variables
    @constraint(sprob, ep_varConstr[i in fData.IDList, t in 2:T], -μbat_out[i,1,t] + sum(λbat_eff[i,l,t] for l in 1:length(bData.ηα[i])) - λpi[i,t] == 0);
    @constraint(sprob, eq_varConstr[i in fData.IDList, t in 2:T], -μbat_out[i,2,t] - λqi[i,t] == 0);
    @constraint(sprob, f_varConstr[i in fData.IDList, t in 2:T], fData.Δt * λbat_trans[i,t] - sum(bData.ηα[i][l] * λbat_eff[i,l,t] for l in 1:length(bData.ηα[i])) == 0);
    @constraint(sprob, I_varConstr[i in fData.IDList, t in 2:(T-1)], λbat_lim[i,t] + λbat_trans[i,t] - λbat_trans[i,t+1] >= 0);
    @constraint(sprob, I_varConstrT[i in fData.IDList], λbat_lim[i,T] + λbat_trans[i,T] >= 0);
    @constraint(sprob, I_varConstr1[i in fData.IDList], λbat_lim[i,1] - λbat_trans[i,2] >= 0);

    # set up constraints for demand/renewable variables
    @constraint(sprob, dp_varConstr[i in fData.IDList, t in 2:T], λpi[i,t] + λdp[i,t] == 0);
    @constraint(sprob, dq_varConstr[i in fData.IDList, t in 2:T], λqi[i,t] + λdq[i,t] == 0);
    @constraint(sprob, h_varConstr[i in fData.IDList, t in 2:T], -λpi[i,t] + λhp[i,t] >= 0);
    @constraint(sprob, uzhpConstr[i in hData.hList, t in 2:T], λhzp1[i,t] + λhzp2[i,t] + λhzp3[i,t] - uData[i].RESP0[t] * expansion_factor * uData[i].RESPmax * λhp[i,t] >= 0);
    @constraint(sprob, uzhmConstr[i in hData.hList, t in 2:T], λhzm1[i,t] + λhzm2[i,t] + λhzm3[i,t] + uData[i].RESP0[t] * expansion_factor * uData[i].RESPmin * λhp[i,t] >= 0);

    # set up constraints for x/y/z
    @constraint(sprob, xConstr1[k in fData.brList, t in 2:T; fData.rateA[k] < Inf], M[k] * (λptrans1[k,t] - λptrans2[k,t] + λqtrans1[k,t] - λqtrans2[k,t])
                 - fData.rateA[k] * ν1[k,t] + (λθ1[k,t] - λθ2[k,t])*2*pi + λxIni[k,t] == 0);
    @constraint(sprob, xConstr2[k in fData.brList, t in 2:T; fData.rateA[k] == Inf], M[k] * (λptrans1[k,t] - λptrans2[k,t] + λqtrans1[k,t] - λqtrans2[k,t])
                 + (λθ1[k,t] - λθ2[k,t])*2*pi + λxIni[k,t] == 0);
    @constraint(sprob, yConstr[i in fData.IDList], -sum(bData.cap[i]*λbat_lim[i,t] + νbat_out[i,t]*bData.uCap[i] for t in 2:T) - bData.cap[i]*λbat_lim[i,1] + λyIni[i] == 0);
    @constraint(sprob, zConstr[i in hData.hList], -sum(λhzp2[i,t] + λhzp3[i,t] + λhzm2[i,t] + λhzm3[i,t] + λhp[i,t] * uData[i].RESP0[t] * expansion_factor for t in 2:T)
                 + λzIni[i] == 0);

    # set up constraints for sp/sq
    @constraint(sprob, spConstr[i in fData.genIDList, t in 2:(T-1)], λspu[i,t] + λspl[i,t] + λrampu[i,t] - λrampu[i,t+1] + λrampl[i,t] - λrampl[i,t+1] - λpi[fData.Loc[i],t] ==
                -fData.cp[i].params[2]);
    @constraint(sprob, sqConstr[i in fData.genIDList, t in 2:T], λsqu[i,t] + λsql[i,t] - λqi[fData.Loc[i],t] == 0);
    @constraint(sprob, spConstr1[i in fData.genIDList], λspu[i,1] + λspl[i,1] + λspIni[i] - λrampu[i,2] - λrampl[i,2] == 0);
    @constraint(sprob, sqConstr1[i in fData.genIDList], λsqu[i,1] + λsql[i,1] + λsqIni[i] == 0);
    @constraint(sprob, spConstrT[i in fData.genIDList], λspu[i,T] + λspl[i,T] + λrampu[i,T] + λrampl[i,T] - λpi[fData.Loc[i],T] == -fData.cp[i].params[2]);

    # set up the objective function
    @objective(sprob, Max, -sum(sum(-μ4[i,2,t]/4 + ν4[i,t]/4 + λvu[i,t]*vmaxT[t][i] + λvl[i,t]*vminT[t][i] - vmaxT[t][i]*vminT[t][i]*λv[i,t] 
        + sum(bData.ηβ[i][l] * λbat_eff[i,l,t] for l in 1:length(bData.ηα[i])) + λdp[i,t]*uData[i].DP0[t] + λdq[i,t]*uData[i].DQ0[t]
        + (uData[i].DPmax[t] - uData[i].DP0[t])*rdpp[i,t] + (uData[i].DPmin[t] - uData[i].DP0[t])*rdpm[i,t] 
        + (uData[i].DQmax[t] - uData[i].DQ0[t])*rdqp[i,t] + (uData[i].DQmin[t] - uData[i].DQ0[t])*rdqm[i,t] 
            for t in 2:T) + yhat[i]*λyIni[i] for i in fData.IDList) - 
        sum(sum(fData.Pmax[i] * λspu[i,t] + fData.Pmin[i] * λspl[i,t] + fData.Qmax[i] * λsqu[i,t] + fData.Qmin[i] * λsql[i,t] 
        + λrampu[i,t] * fData.RU[i] + λrampl[i,t] * fData.RD[i] for t in 2:T)
        + fData.Pmax[i] * λspu[i,1] + fData.Pmin[i] * λspl[i,1] + fData.Qmax[i] * λsqu[i,1] + fData.Qmin[i] * λsql[i,1] 
        + sphat[i] * λspIni[i] + sqhat[i] * λsqIni[i] for i in fData.genIDList) -
        sum(-3/4*μ3[k,2,t] + 5/4*ν3[k,t] + csmax[k,t]*λcs1[k,t] + csmin[k,t]*λcs2[k,t] + (cos(θDminT[t][k]) - csConst[k,t]*(fData.σ[k]+θDminT[t][k]))*λcs3[k,t]
        + (sin(θu[t][k]/2) - θu[t][k]/2*cos(θu[t][k]/2))*(λss4[k,t] - λss5[k,t]) + ssmax[k,t]*λss1[k,t] + ssmin[k,t]*λss2[k,t]
        + (θDmaxT[t][k] + fData.σ[k] + 2*pi)*λθ1[k,t] + (θDminT[t][k] + fData.σ[k] - 2*pi)*λθ2[k,t] - fData.σ[k]*(λss4[k,t] + λss5[k,t])
        + M[k]*(λptrans1[k,t] - λptrans2[k,t] + λqtrans1[k,t] - λqtrans2[k,t]) + xhat[k,t]*λxIni[k,t]
        + vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]] - vmaxT[t][k[1]]*vmaxT[t][k[2]])/(fData.τ1[k]*fData.τ2[k])*λlnc1[k,t]
        - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]] - vmaxT[t][k[1]]*vmaxT[t][k[2]])/(fData.τ1[k]*fData.τ2[k])*λlnc2[k,t]
        - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv1[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv2[k,t]
        - vminT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv3[k,t] - vmaxT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λvv4[k,t]
        - csmin[k,t]*vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc1[k,t] - csmax[k,t]*vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc2[k,t]
        - csmax[k,t]*vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc3[k,t] - csmin[k,t]*vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λwc4[k,t]
        - ssmin[k,t]*vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λws1[k,t] - ssmax[k,t]*vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λws2[k,t]
        - ssmax[k,t]*vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λws3[k,t] - ssmin[k,t]*vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*λws4[k,t]
           for t in 2:T for k in fData.brList) - 
        sum(sum(λhp[i,t]*uData[i].RESP0[t] - λhzp3[i,t] - λhzm3[i,t] + rhpp[i,t]*uData[i].RESPmax*uData[i].RESP0[t] - rhpm[i,t]*uData[i].RESPmin*uData[i].RESP0[t]
        + rhppz1[i,t] + rhppz2[i,t] + rhpmz1[i,t] + rhpmz2[i,t] for t in 2:T) + zhat[i]*λzIni[i] for i in hData.hList)
            );

    return sprob;
end

# obtain the dual information for the test_subdual
function obtain_dual(sprob,fData,hData,T,groupDict)
    groupList = 1:length(groupDict[2]);
    group_rev = Dict();
    for i in fData.IDList
        for m in groupList
            if i in groupDict[2][m]
                group_rev[i] = m;
            end
        end
    end

    u_hp_coeff_mat = -(dual.(sprob[:rhppCons1]) .* fData.cz + dual.(sprob[:rhppCons2]) .* (-fData.cz) + 
                     dual.(sprob[:rhppCons3]) .* (-5 * fData.cz) + dual.(sprob[:rhppCons4]) .* (5 * fData.cz) +
                     dual.(sprob[:rhppzCons11]) .* fData.cz + dual.(sprob[:rhppzCons12]) .* (-fData.cz) + 
                     dual.(sprob[:rhppzCons13]) .* (-5 * fData.cz) + dual.(sprob[:rhppzCons14]) .* (5 * fData.cz) +
                     dual.(sprob[:rhppzCons21]) .* fData.cz + dual.(sprob[:rhppzCons22]) .* (-fData.cz) + 
                     dual.(sprob[:rhppzCons23]) .* (-5 * fData.cz) + dual.(sprob[:rhppzCons24]) .* (5 * fData.cz));
    u_hm_coeff_mat = -(dual.(sprob[:rhpmCons1]) .* fData.cz + dual.(sprob[:rhpmCons2]) .* (-fData.cz) + 
                    dual.(sprob[:rhpmCons3]) .* (-5 * fData.cz) + dual.(sprob[:rhpmCons4]) .* (5 * fData.cz) +
                    dual.(sprob[:rhpmzCons11]) .* fData.cz + dual.(sprob[:rhpmzCons12]) .* (-fData.cz) + 
                    dual.(sprob[:rhpmzCons13]) .* (-5 * fData.cz) + dual.(sprob[:rhpmzCons14]) .* (5 * fData.cz) +
                    dual.(sprob[:rhpmzCons21]) .* fData.cz + dual.(sprob[:rhpmzCons22]) .* (-fData.cz) + 
                    dual.(sprob[:rhpmzCons23]) .* (-5 * fData.cz) + dual.(sprob[:rhpmzCons24]) .* (5 * fData.cz));
    u_dp_coeff_mat = -(dual.(sprob[:rdppCons1]) .* fData.cz + dual.(sprob[:rdppCons2]) .* (-fData.cz) + 
                    dual.(sprob[:rdppCons3]) .* (-5 * fData.cz) + dual.(sprob[:rdppCons4]) .* (5 * fData.cz) +
                    dual.(sprob[:rdqpCons1]) .* fData.cz + dual.(sprob[:rdqpCons2]) .* (-fData.cz) + 
                    dual.(sprob[:rdqpCons3]) .* (-5 * fData.cz) + dual.(sprob[:rdqpCons4]) .* (5 * fData.cz));
    u_dm_coeff_mat = -(dual.(sprob[:rdpmCons1]) .* fData.cz + dual.(sprob[:rdpmCons2]) .* (-fData.cz) + 
                    dual.(sprob[:rdpmCons3]) .* (-5 * fData.cz) + dual.(sprob[:rdpmCons4]) .* (5 * fData.cz) +
                    dual.(sprob[:rdqmCons1]) .* fData.cz + dual.(sprob[:rdqmCons2]) .* (-fData.cz) + 
                    dual.(sprob[:rdqmCons3]) .* (-5 * fData.cz) + dual.(sprob[:rdqmCons4]) .* (5 * fData.cz));
    
    u_hp_coeff = Dict();
    u_hm_coeff = Dict();
    u_dp_coeff = Dict();
    u_dm_coeff = Dict();

    for t in 2:T
        for i in hData.hList
            u_hp_coeff[i,t] = u_hp_coeff_mat[i,t];
            u_hm_coeff[i,t] = u_hm_coeff_mat[i,t];
        end
        for ig in groupList
            u_dp_coeff[ig,t] = sum(u_dp_coeff_mat[i,t] for i in groupDict[2][ig]);
            u_dm_coeff[ig,t] = sum(u_dm_coeff_mat[i,t] for i in groupDict[2][ig]);
        end
    end
    
    return u_hp_coeff, u_hm_coeff, u_dp_coeff, u_dm_coeff;
end

function dual_sub_uhat(fData,uData,hData,T,groupDict,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uhat)
    groupList = 1:length(groupDict[2]);
    group_rev = Dict();
    for i in fData.IDList
        for m in groupList
            if i in groupDict[2][m]
                group_rev[i] = m;
            end
        end
    end

    # dual's dual should be similar to the primal
    θu = Dict();
    for t in 1:T
        θu[t] = Dict();
        for k in fData.brList
            θu[t][k] = max(abs(θDmaxT[t][k]),abs(θDminT[t][k]));
        end
    end
    # scaling issue exists for this primal problem.
    subd = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0, "Threads" => 1));
    # subd = Model(COPT.Optimizer);

    # obtain the pairs that are connected
    connectPair = [];
    connectDict = Dict();
    branchDict1 = Dict();
    branchDict2 = Dict();
    for i in fData.IDList
        connectDict[i] = [];
        branchDict1[i] = [];
        branchDict2[i] = [];
    end
    for k in fData.brList
        push!(branchDict1[k[1]],k);
        push!(branchDict2[k[2]],k);
        if !((k[1],k[2]) in connectPair)
            push!(connectPair,(k[1],k[2]));
            push!(connectDict[k[1]],k[2]);
        end
    end
    kpDict = Dict();
    for k in fData.brList
        if (k[1],k[2]) in keys(kpDict)
            push!(kpDict[(k[1],k[2])],k);
        else
            kpDict[(k[1],k[2])] = [k];
        end
    end

    # violation slack variables
    @variable(subd,lpplus[i in fData.IDList, t in 2:T] >= 0);
    @variable(subd,lpminus[i in fData.IDList, t in 2:T] >= 0);
    @variable(subd,lqplus[i in fData.IDList, t in 2:T] >= 0);
    @variable(subd,lqminus[i in fData.IDList, t in 2:T] >= 0);

    # OPF variables
    @variable(subd,p[k in fData.brList, t in 2:T]);
    @variable(subd,q[k in fData.brList, t in 2:T]);
    @variable(subd,vminT[t][i] <= v[i in fData.IDList, t in 2:T] <= vmaxT[t][i]);
    @variable(subd,vhat[i in fData.IDList, t in 2:T]);
    @variable(subd,θ[i in fData.IDList, t in 2:T]);

    # switching, battery setup, renewable expansion binary variables
    @variable(subd,x[k in fData.brList,t in 2:T]);
    @variable(subd,y[i in fData.IDList]);
    @variable(subd,z[i in hData.hList]);
    
    # set up the reference bus angle = 0
    for i in fData.IDList
        for t in 2:T
            if fData.bType[i] == 3
                @constraint(subd,θ[i,t] == 0);
            end
        end
    end

    # set up the variables: active/reactive power injection
    @variable(subd, fData.Pmin[i] <= sp[i in fData.genIDList, t in 1:T] <= fData.Pmax[i]);
    @variable(subd, fData.Qmin[i] <= sq[i in fData.genIDList, t in 1:T] <= fData.Qmax[i]);
    @constraint(subd, spIni[i in fData.genIDList], sp[i,1] == sphat[i]);
    @constraint(subd, sqIni[i in fData.genIDList], sq[i,1] == sqhat[i]);
    # ramping constraints
    @constraint(subd,ramp_up[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t-1] <= fData.RU[i]);
    @constraint(subd,ramp_dn[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t-1] >= fData.RD[i]);

    sphatsum = Dict();
    sqhatsum = Dict();
    for t in 2:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(subd,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end

        for i in fData.IDList
            sqhatsum[i,t] = @expression(subd,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(subd,tAux2[k in fData.brList, t in 2:T]);
    @variable(subd,tAux3[k in fData.brList, t in 2:T]);
    @variable(subd,tAux4[k in fData.brList, t in 2:T] >= 0);
    @variable(subd,tAux5[i in fData.IDList, t in 2:T]);
    @variable(subd,tAux6[i in fData.IDList, t in 2:T] >= 0);

    @variable(subd,vv[k in fData.brList, t in 2:T]);
    @variable(subd,cos(θu[t][k]) <= cs[k in fData.brList, t in 2:T] <= 1);
    @variable(subd,-sin(θu[t][k]) <= ss[k in fData.brList, t in 2:T] <= sin(θu[t][k]));
    # add the strengthened bounds on cs
    csmax = Dict();
    csmin = Dict();
    ssmax = Dict();
    ssmin = Dict();
    for t in 2:T
        for k in fData.brList
            if (θDmaxT[t][k] >= 0)&(θDminT[t][k] >= 0)
                csmax[k,t] = cos(θDminT[t][k]);
                csmin[k,t] = cos(θDmaxT[t][k]);
                @constraint(subd, cs[k,t] <= csmax[k,t]);
                @constraint(subd, cs[k,t] >= csmin[k,t]);
            elseif (θDmaxT[t][k] < 0)&(θDminT[t][k] < 0)
                csmax[k,t] = cos(θDmaxT[t][k]);
                csmin[k,t] = cos(θDminT[t][k]);
                @constraint(subd, cs[k,t] <= csmax[k,t]);
                @constraint(subd, cs[k,t] >= csmin[k,t]);
            else
                csmax[k,t] = 1;
                csmin[k,t] = min(cos(θDmaxT[t][k]),cos(θDminT[t][k]));
                @constraint(subd, cs[k,t] <= csmax[k,t]);
                @constraint(subd, cs[k,t] >= csmin[k,t]);
            end
            @constraint(subd, cs[k,t] >= (cos(θDmaxT[t][k]) - cos(θDminT[t][k]))/(θDmaxT[t][k] - θDminT[t][k])*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t] - θDminT[t][k]) + cos(θDminT[t][k]));
        end
        # add the strengthened bounds on ss
        for k in fData.brList
            ssmax[k,t] = sin(θDmaxT[t][k]);
            ssmin[k,t] = sin(θDminT[t][k]);
            @constraint(subd, ss[k,t] <= ssmax[k,t]);
            @constraint(subd, ss[k,t] >= ssmin[k,t]);
        end
    end

    @variable(subd,wc[k in fData.brList, t in 2:T]);
    @variable(subd,ws[k in fData.brList, t in 2:T]);
    @constraint(subd,wcEquality[k in fData.brList, t in 2:T;k[1] < k[2]], wc[k,t] == wc[(k[2],k[1],k[3]),t]);
    @constraint(subd,wsEquality[k in fData.brList, t in 2:T;k[1] < k[2]], ws[k,t] == -ws[(k[2],k[1],k[3]),t]);
    @constraint(subd,vvEquality[k in fData.brList, t in 2:T;k[1] < k[2]], vv[k,t] == vv[(k[2],k[1],k[3]),t]);
    @constraint(subd,csEquality[k in fData.brList, t in 2:T;k[1] < k[2]], cs[k,t] == cs[(k[2],k[1],k[3]),t]);
    @constraint(subd,ssEquality[k in fData.brList, t in 2:T;k[1] < k[2]], ss[k,t] == -ss[(k[2],k[1],k[3]),t]);

    M = Dict();
    for k in fData.brList
        if fData.rateA[k] < Inf
            M[k] = fData.rateA[k]^2;
        else
            M[k] = 10000;
        end
    end
    @constraint(subd,lineConstrP1[k in fData.brList, t in 2:T],p[k,t] <= (fData.g[k]*vhat[k[1],t]/(fData.τ1[k]^2) - fData.g[k]*wc[k,t] - fData.b[k]*ws[k,t]) + M[k]*(1 - x[k,t]));
    @constraint(subd,lineConstrP2[k in fData.brList, t in 2:T],p[k,t] >= (fData.g[k]*vhat[k[1],t]/(fData.τ1[k]^2) - fData.g[k]*wc[k,t] - fData.b[k]*ws[k,t]) - M[k]*(1 - x[k,t]));
    @constraint(subd,lineConstrQ1[k in fData.brList, t in 2:T],q[k,t] <= ((-fData.b[k] - fData.bc[k]/2)*vhat[k[1],t]/(fData.τ1[k]^2) + fData.b[k]*wc[k,t] - fData.g[k]*ws[k,t]) + M[k]*(1 - x[k,t]));
    @constraint(subd,lineConstrQ2[k in fData.brList, t in 2:T],q[k,t] >= ((-fData.b[k] - fData.bc[k]/2)*vhat[k[1],t]/(fData.τ1[k]^2) + fData.b[k]*wc[k,t] - fData.g[k]*ws[k,t]) - M[k]*(1 - x[k,t]));
    socList1 = Dict();
    socList2 = Dict();
    socList3 = Dict();
    socList4 = Dict();
    socList5 = Dict();
    for t in 2:T
        for k in fData.brList
            if fData.rateA[k] < Inf
                socList1[k,t] = [p[k,t],q[k,t]];
            end
            socList2[k,t] = [vv[k,t],vhat[k[1],t]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2],t]/((fData.τ2[k]^2)*sqrt(2))];
            socList3[k,t] = [tAux2[k,t],tAux3[k,t]];
            socList5[k,t] = [wc[k,t],ws[k,t],vhat[k[1],t]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2],t]/((fData.τ2[k]^2)*sqrt(2))];
        end
        for i in fData.IDList
            socList4[i,t] = [v[i,t],tAux5[i,t]];
        end
    end

    @constraint(subd,socConstraint1[k in fData.brList, t in 2:T; fData.rateA[k] < Inf],[fData.rateA[k] * x[k,t]; socList1[k,t]] in SecondOrderCone());
    @constraint(subd,socConstraint2[k in fData.brList, t in 2:T], [(vhat[k[1],t]/(fData.τ1[k]^2) + vhat[k[2],t]/(fData.τ2[k]^2))/sqrt(2); socList2[k,t]] in SecondOrderCone());
    @constraint(subd,socConstraint3[k in fData.brList, t in 2:T],[tAux4[k,t]; socList3[k,t]] in SecondOrderCone());
    @constraint(subd,socConstraint4[i in fData.IDList, t in 2:T],[tAux6[i,t]; socList4[i,t]] in SecondOrderCone());
    @constraint(subd,socConstraint5[k in fData.brList, t in 2:T],[(vhat[k[1],t]/(fData.τ1[k]^2) + vhat[k[2],t]/(fData.τ2[k]^2))/sqrt(2); socList5[k,t]] in SecondOrderCone());

    @constraint(subd,auxConstr2[k in fData.brList, t in 2:T],sqrt((1-cos(θu[t][k]))/(θu[t][k])^2)*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t]) == tAux2[k,t]);
    @constraint(subd,auxConstr3[k in fData.brList, t in 2:T],cs[k,t] - 3/4 == tAux3[k,t]);
    @constraint(subd,auxConstr4[k in fData.brList, t in 2:T],5/4 - cs[k,t] == tAux4[k,t]);
    @constraint(subd,sinMock1[k in fData.brList, t in 2:T],
                ss[k,t] <= cos(θu[t][k]/2)*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t] - θu[t][k]/2) + sin(θu[t][k]/2));
    @constraint(subd,sinMock2[k in fData.brList, t in 2:T],
                ss[k,t] >= cos(θu[t][k]/2)*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t] + θu[t][k]/2) - sin(θu[t][k]/2));
    @constraint(subd,angleDiff1[k in fData.brList, t in 2:T],(θ[k[1],t] - fData.σ[k]) - θ[k[2],t] <= θDmaxT[t][k] + 2*pi*(1 - x[k,t]));
    @constraint(subd,angleDiff2[k in fData.brList, t in 2:T],(θ[k[1],t] - fData.σ[k]) - θ[k[2],t] >= θDminT[t][k] - 2*pi*(1 - x[k,t]));
    @constraint(subd,auxConstr5[i in fData.IDList, t in 2:T],tAux5[i,t] == vhat[i,t] - 1/4);
    @constraint(subd,auxConstr6[i in fData.IDList, t in 2:T],tAux6[i,t] == vhat[i,t] + 1/4);
    @constraint(subd,v2Mock2[i in fData.IDList, t in 2:T], vhat[i,t] - (vmaxT[t][i] + vminT[t][i])*v[i,t] <= -vmaxT[t][i]*vminT[t][i]);
    @constraint(subd,vvMock1[k in fData.brList, t in 2:T],
                vv[k,t] >= vminT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vminT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));
    @constraint(subd,vvMock2[k in fData.brList, t in 2:T],
                vv[k,t] >= vmaxT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vmaxT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));
    @constraint(subd,vvMock3[k in fData.brList, t in 2:T],
                vv[k,t] <= vminT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vmaxT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vminT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));
    @constraint(subd,vvMock4[k in fData.brList, t in 2:T],
                vv[k,t] <= vmaxT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vminT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));

    @constraint(subd,wcMock1[k in fData.brList, t in 2:T],
                wc[k,t] >= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + csmin[k,t]*vv[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k,t]);
    @constraint(subd,wcMock2[k in fData.brList, t in 2:T],
                wc[k,t] >= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + csmax[k,t]*vv[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k,t]);
    @constraint(subd,wcMock3[k in fData.brList, t in 2:T],
                wc[k,t] <= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + vv[k,t]*csmax[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k,t]);
    @constraint(subd,wcMock4[k in fData.brList, t in 2:T],
                wc[k,t] <= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + vv[k,t]*csmin[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k,t]);

    @constraint(subd,wsMock1[k in fData.brList, t in 2:T],
                ws[k,t] >= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + ssmin[k,t]*vv[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k,t]);
    @constraint(subd,wsMock2[k in fData.brList, t in 2:T],
                ws[k,t] >= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + ssmax[k,t]*vv[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k,t]);
    @constraint(subd,wsMock3[k in fData.brList, t in 2:T],
                ws[k,t] <= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + vv[k,t]*ssmax[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k,t]);
    @constraint(subd,wsMock4[k in fData.brList, t in 2:T],
                ws[k,t] <= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + vv[k,t]*ssmin[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k,t]);

    @constraint(subd,tanConstr1[k in fData.brList, t in 2:T],
                    ws[k,t] - tan(θDmaxT[t][k])*wc[k,t] <= 0);
    @constraint(subd,tanConstr2[k in fData.brList, t in 2:T],
                    ws[k,t] - tan(θDminT[t][k])*wc[k,t] >= 0);

    vδ = Dict();
    θϕ = Dict();
    θδ = Dict();
    for t in 2:T
        for i in fData.IDList
            vδ[i,t] = vmaxT[t][i]+vminT[t][i];
        end
        for k in fData.brList
            θϕ[k,t] = (θDmaxT[t][k] + θDminT[t][k])/2;
            θδ[k,t] = (θDmaxT[t][k] - θDminT[t][k])/2;
        end
    end

    @constraint(subd,lncConstr1[k in fData.brList, t in 2:T],
                vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*(wc[k,t]*cos(θϕ[k,t]) + ws[k,t]*sin(θϕ[k,t])) -
                vmaxT[t][k[2]]/fData.τ2[k]*cos(θδ[k,t])*vδ[k[2],t]/fData.τ2[k]*vhat[k[1],t]/(fData.τ1[k]^2) -
                vmaxT[t][k[1]]/fData.τ1[k]*cos(θδ[k,t])*vδ[k[1],t]/fData.τ1[k]*vhat[k[2],t]/(fData.τ2[k]^2) >=
                vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]) -
                vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])));
    @constraint(subd,lncConstr2[k in fData.brList, t in 2:T],
                vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*(wc[k,t]*cos(θϕ[k,t]) + ws[k,t]*sin(θϕ[k,t])) -
                vminT[t][k[2]]/fData.τ2[k]*cos(θδ[k,t])*vδ[k[2],t]/fData.τ2[k]*vhat[k[1],t]/(fData.τ1[k]^2) -
                vminT[t][k[1]]/fData.τ1[k]*cos(θδ[k,t])*vδ[k[1],t]/fData.τ1[k]*vhat[k[2],t]/(fData.τ2[k]^2) >=
                -vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]) -
                vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])));

    # battery dynamics
    @variable(subd, ep_var[i in fData.IDList, t in 2:T]);
    @variable(subd, eq_var[i in fData.IDList, t in 2:T]);
    @variable(subd, f_var[i in fData.IDList, t in 2:T]);
    @variable(subd, I_var[i in fData.IDList, t in 1:T] >= 0);
    @constraint(subd, soc_limit[i in fData.IDList, t in 1:T], I_var[i,t] <= bData.cap[i]*y[i]);
    @constraint(subd, soc_trans[i in fData.IDList, t in 2:T], I_var[i,t] == I_var[i,t-1] - fData.Δt * f_var[i,t]);
    @constraint(subd, bat_out[i in fData.IDList, t in 2:T], [bData.uCap[i]*y[i]; [ep_var[i,t], eq_var[i,t]]] in SecondOrderCone());
    @constraint(subd, bat_eff[i in fData.IDList, l in 1:length(bData.ηα[i]), t in 2:T], ep_var[i,t] - bData.ηα[i][l]*f_var[i,t] <= bData.ηβ[i][l]);

    # power flow balance
    @variable(subd, dp_var[i in fData.IDList, t in 2:T]);
    @variable(subd, dq_var[i in fData.IDList, t in 2:T]);
    @variable(subd, h_var[i in fData.IDList, t in 2:T] >= 0);
    @variable(subd, zetap_var[i in hData.hList, t in 2:T; i in hData.hList] >= 0);
    @variable(subd, zetam_var[i in hData.hList, t in 2:T; i in hData.hList] >= 0);
    @constraint(subd,totalP[i in fData.IDList, t in 2:T], sum(p[k,t] for k in branchDict1[i]) + vhat[i,t]*fData.gs[i]
                + lpplus[i,t] - lpminus[i,t] - h_var[i,t] + dp_var[i,t] - ep_var[i,t] == sphatsum[i,t]);
    @constraint(subd,totalQ[i in fData.IDList, t in 2:T], sum(q[k,t] for k in branchDict1[i]) - vhat[i,t]*fData.bs[i]
                + lqplus[i,t] - lqminus[i,t]  + dq_var[i,t] - eq_var[i,t] == sqhatsum[i,t]);

    # bind the previous solutions
    @constraint(subd,xIni[k in fData.brList, t in 2:T],x[k,t] == xhat[k,t]);
    @constraint(subd,yIni[i in fData.IDList],y[i] == yhat[i]);
    @constraint(subd,zIni[i in hData.hList],z[i] == zhat[i]);
    u_dp = uhat["u_dp"];
    u_dm = uhat["u_dm"];
    u_hp = uhat["u_hp"];
    u_hm = uhat["u_hm"];

    # set up the variables for bilinear terms in dual formulation
    @variable(subd, pihp[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList]);
    @variable(subd, pihm[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList]);
    @variable(subd, pizetap2[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList]);
    @variable(subd, pizetam2[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList]);
    @variable(subd, pizetap3[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList]);
    @variable(subd, pizetam3[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList]);
    @variable(subd, pipdp[i in fData.IDList, t in 2:T, j in 1:4]);
    @variable(subd, pipdm[i in fData.IDList, t in 2:T, j in 1:4]);
    @variable(subd, piqdp[i in fData.IDList, t in 2:T, j in 1:4]);
    @variable(subd, piqdm[i in fData.IDList, t in 2:T, j in 1:4]);

    signDict = Dict(1 => -1, 2 => 1, 3 => -1, 4 => 1);
    u_hp_signDict = Dict(1 => Dict(), 2 => Dict(), 3 => Dict(), 4 => Dict());
    u_hm_signDict = Dict(1 => Dict(), 2 => Dict(), 3 => Dict(), 4 => Dict());
    u_dp_signDict = Dict(1 => Dict(), 2 => Dict(), 3 => Dict(), 4 => Dict());
    u_dm_signDict = Dict(1 => Dict(), 2 => Dict(), 3 => Dict(), 4 => Dict());
    for t in 2:T
        for i in hData.hList
            u_hp_signDict[1][i,t] = 1 - u_hp[i,t];
            u_hp_signDict[2][i,t] = 1 - u_hp[i,t];
            u_hp_signDict[3][i,t] = u_hp[i,t];
            u_hp_signDict[4][i,t] = u_hp[i,t];
            u_hm_signDict[1][i,t] = 1 - u_hm[i,t];
            u_hm_signDict[2][i,t] = 1 - u_hm[i,t];
            u_hm_signDict[3][i,t] = u_hm[i,t];
            u_hm_signDict[4][i,t] = u_hm[i,t];
        end
        for m in groupList
            u_dp_signDict[1][m,t] = 1 - u_dp[m,t];
            u_dp_signDict[2][m,t] = 1 - u_dp[m,t];
            u_dp_signDict[3][m,t] = u_dp[m,t];
            u_dp_signDict[4][m,t] = u_dp[m,t];
            u_dm_signDict[1][m,t] = 1 - u_dm[m,t];
            u_dm_signDict[2][m,t] = 1 - u_dm[m,t];
            u_dm_signDict[3][m,t] = u_dm[m,t];
            u_dm_signDict[4][m,t] = u_dm[m,t];
        end
    end

    @constraint(subd, pihp_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihp[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pihp_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihp[i,t,j] * signDict[j] <= abs(uData[i].RESP0[t] * uData[i].RESPmax) * u_hp_signDict[j][i,t]);
    @constraint(subd, pihm_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihm[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pihm_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihm[i,t,j] * signDict[j] <= abs(uData[i].RESP0[t] * uData[i].RESPmin) * u_hm_signDict[j][i,t]);
    @constraint(subd, pizetap2_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap2[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetap2_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap2[i,t,j] * signDict[j] <= u_hp_signDict[j][i,t]);
    @constraint(subd, pizetam2_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam2[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetam2_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam2[i,t,j] * signDict[j] <= u_hm_signDict[j][i,t]);
    @constraint(subd, pizetap3_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap3[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetap3_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap3[i,t,j] * signDict[j] <= u_hp_signDict[j][i,t]);
    @constraint(subd, pizetam3_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam3[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetam3_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam3[i,t,j] * signDict[j] <= u_hm_signDict[j][i,t]);

    @constraint(subd, pipdp_sign[i in fData.IDList, t in 2:T, j in 1:4], pipdp[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pipdm_sign[i in fData.IDList, t in 2:T, j in 1:4], pipdm[i,t,j] * signDict[j] >= 0);
    @constraint(subd, piqdp_sign[i in fData.IDList, t in 2:T, j in 1:4], piqdp[i,t,j] * signDict[j] >= 0);
    @constraint(subd, piqdm_sign[i in fData.IDList, t in 2:T, j in 1:4], piqdm[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pipdp_bound[i in fData.IDList, t in 2:T, j in 1:4], pipdp[i,t,j] * signDict[j] <= abs(uData[i].DPmax[t] - uData[i].DP0[t]) * u_dp_signDict[j][group_rev[i],t]);
    @constraint(subd, pipdm_bound[i in fData.IDList, t in 2:T, j in 1:4], pipdm[i,t,j] * signDict[j] <= abs(-uData[i].DPmin[t] + uData[i].DP0[t]) * u_dm_signDict[j][group_rev[i],t]);
    @constraint(subd, piqdp_bound[i in fData.IDList, t in 2:T, j in 1:4], piqdp[i,t,j] * signDict[j] <= abs(uData[i].DQmax[t] - uData[i].DQ0[t]) * u_dp_signDict[j][group_rev[i],t]);
    @constraint(subd, piqdm_bound[i in fData.IDList, t in 2:T, j in 1:4], piqdm[i,t,j] * signDict[j] <= abs(-uData[i].DQmin[t] + uData[i].DQ0[t]) * u_dm_signDict[j][group_rev[i],t]);

    # d_var constraints: demand uncertainty
    @constraint(subd, dp_cal[i in fData.IDList, t in 2:T], dp_var[i,t] == pipdp[i,t,3] + pipdp[i,t,4] + pipdm[i,t,3] + pipdm[i,t,4] + uData[i].DP0[t]);
    @constraint(subd, dq_cal[i in fData.IDList, t in 2:T], dq_var[i,t] == piqdp[i,t,3] + piqdp[i,t,4] + piqdm[i,t,3] + piqdm[i,t,4] + uData[i].DQ0[t]);

    # h_var constraints: renewable uncertainty
    @constraint(subd, hp_cal1[i in fData.IDList, t in 2:T; i in hData.hList], h_var[i,t] <= uData[i].RESP0[t]*(1 + expansion_factor * z[i] + 
            uData[i].RESPmax * expansion_factor * zetap_var[i,t] - uData[i].RESPmin * expansion_factor * zetam_var[i,t]) + 
            pihp[i,t,3] + pihp[i,t,4] + pihm[i,t,3] + pihm[i,t,4]);
    @constraint(subd, hp_cal2[i in fData.IDList, t in 2:T; !(i in hData.hList)], h_var[i,t] == 0.0); 
    @constraint(subd, zetap_mc1[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], zetap_var[i,t] <= z[i]);
    @constraint(subd, zetap_mc2[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], zetap_var[i,t] <= pizetap2[i,t,3] + pizetap2[i,t,4]);
    @constraint(subd, zetap_mc3[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], zetap_var[i,t] >= pizetap3[i,t,3] + pizetap3[i,t,4] + z[i] - 1);
    @constraint(subd, zetam_mc1[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], zetam_var[i,t] <= z[i]);
    @constraint(subd, zetam_mc2[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], zetam_var[i,t] <= pizetam2[i,t,3] + pizetam2[i,t,4]);
    @constraint(subd, zetam_mc3[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], zetam_var[i,t] >= pizetam3[i,t,3] + pizetam3[i,t,4] + z[i] - 1);

    # pi constraints
    @constraint(subd, pihp_constr[i in fData.IDList, t in 2:T; i in hData.hList], sum(pihp[i,t,j] for j in 1:4) == uData[i].RESP0[t] * uData[i].RESPmax);
    @constraint(subd, pihm_constr[i in fData.IDList, t in 2:T; i in hData.hList], sum(pihm[i,t,j] for j in 1:4) == -uData[i].RESP0[t] * uData[i].RESPmin);
    @constraint(subd, pizetap2_constr[i in fData.IDList, t in 2:T; i in hData.hList], sum(pizetap2[i,t,j] for j in 1:4) == 1);
    @constraint(subd, pizetam2_constr[i in fData.IDList, t in 2:T; i in hData.hList], sum(pizetam2[i,t,j] for j in 1:4) == 1);
    @constraint(subd, pizetap3_constr[i in fData.IDList, t in 2:T; i in hData.hList], sum(pizetap3[i,t,j] for j in 1:4) == 1);
    @constraint(subd, pizetam3_constr[i in fData.IDList, t in 2:T; i in hData.hList], sum(pizetam3[i,t,j] for j in 1:4) == 1);
    @constraint(subd, pipdp_constr[i in fData.IDList, t in 2:T], sum(pipdp[i,t,j] for j in 1:4) == uData[i].DPmax[t] - uData[i].DP0[t]);
    @constraint(subd, pipdm_constr[i in fData.IDList, t in 2:T], sum(pipdm[i,t,j] for j in 1:4) == uData[i].DPmin[t] - uData[i].DP0[t]);
    @constraint(subd, piqdp_constr[i in fData.IDList, t in 2:T], sum(piqdp[i,t,j] for j in 1:4) == uData[i].DQmax[t] - uData[i].DQ0[t]);
    @constraint(subd, piqdm_constr[i in fData.IDList, t in 2:T], sum(piqdm[i,t,j] for j in 1:4) == uData[i].DQmin[t] - uData[i].DQ0[t]);

    # objective value
    @objective(subd, Min, fData.cz*sum(lpplus[i,t] + lpminus[i,t] + lqplus[i,t] + lqminus[i,t] for i in fData.IDList for t in 2:T) + 
        sum(fData.cp[i].params[2]*sp[i,t] for i in fData.genIDList for t in 2:T));

    return subd;
end

function solve_dual_sub_uhat(fData,uData,hData,T,groupDict,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uhat)
    # construct and solve the dual subproblem to suit the parallelization
    subd = dual_sub_uhat(fData,uData,hData,T,groupDict,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uhat);
    optimize!(subd);

    x_coeff = shadow_price.(subd[:xIni]);
    y_coeff = shadow_price.(subd[:yIni]);
    z_coeff = shadow_price.(subd[:zIni]);

    return objective_value(subd), x_coeff, y_coeff, z_coeff;
end
