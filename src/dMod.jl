# create the dual formulation of the inner min of the min-max-min problem
function createSecond_dual(fData,uData,hData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat)
    # first-stage model without any scenarios
    θu = Dict();
    for t in 1:T
        θu[t] = Dict();
        for k in fData.brList
            θu[t][k] = max(abs(θDmaxT[t][k]),abs(θDminT[t][k]));
        end
    end
    # scaling issue exists for this primal problem.
    sprob = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "NumericFocus" => 3, "BarConvTol" => 1e-6, "MIPGap" => 1e-6, "BarQCPConvTol" => 1e-6,
        "OptimalityTol" => 1e-6, "IntFeasTol" => 1e-6, "FeasibilityTol" => 1e-6, "OutputFlag" => 1, "Threads" => 30));

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
            M[k] = fData.rateA[k]^2;
        else
            M[k] = 10000;
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
    @variable(sprob, rhppz1[i in hData.hList, t in 2:T] >= 0);
    @variable(sprob, rhppz2[i in hData.hList, t in 2:T] <= 0);
    @variable(sprob, rhpmz1[i in hData.hList, t in 2:T] >= 0);
    @variable(sprob, rhpmz2[i in hData.hList, t in 2:T] <= 0);

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
    @constraint(sprob, rdppCons3[i in fData.IDList, t in 2:T], rdpp[i,t] <= λdp[i,t] + 2 * (1 - u_dp[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdppCons4[i in fData.IDList, t in 2:T], rdpp[i,t] >= λdp[i,t] - 2 * (1 - u_dp[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rdpmCons1[i in fData.IDList, t in 2:T], rdpm[i,t] <= fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdpmCons2[i in fData.IDList, t in 2:T], rdpm[i,t] >= -fData.cz * u_dm[group_rev[i],t]);
    @constraint(sprob, rdpmCons3[i in fData.IDList, t in 2:T], rdpm[i,t] <= λdp[i,t] + 2 * (1 - u_dm[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdpmCons4[i in fData.IDList, t in 2:T], rdpm[i,t] >= λdp[i,t] - 2 * (1 - u_dm[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rdqpCons1[i in fData.IDList, t in 2:T], rdqp[i,t] <= fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdqpCons2[i in fData.IDList, t in 2:T], rdqp[i,t] >= -fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdqpCons3[i in fData.IDList, t in 2:T], rdqp[i,t] <= λdq[i,t] + 2 * (1 - u_dp[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdqpCons4[i in fData.IDList, t in 2:T], rdqp[i,t] >= λdq[i,t] - 2 * (1 - u_dp[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rdqmCons1[i in fData.IDList, t in 2:T], rdqm[i,t] <= fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdqmCons2[i in fData.IDList, t in 2:T], rdqm[i,t] >= -fData.cz * u_dp[group_rev[i],t]);
    @constraint(sprob, rdqmCons3[i in fData.IDList, t in 2:T], rdqm[i,t] <= λdq[i,t] + 2 * (1 - u_dm[group_rev[i],t]) * fData.cz);
    @constraint(sprob, rdqmCons4[i in fData.IDList, t in 2:T], rdqm[i,t] >= λdq[i,t] - 2 * (1 - u_dm[group_rev[i],t]) * fData.cz);

    @constraint(sprob, rhppCons1[i in hData.hList, t in 2:T], rhpp[i,t] <= fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppCons2[i in hData.hList, t in 2:T], rhpp[i,t] >= -fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppCons3[i in hData.hList, t in 2:T], rhpp[i,t] <= λhp[i,t] + 2 * (1 - u_hp[i,t]) * fData.cz);
    @constraint(sprob, rhppCons4[i in hData.hList, t in 2:T], rhpp[i,t] >= λhp[i,t] - 2 * (1 - u_hp[i,t]) * fData.cz);

    @constraint(sprob, rhpmCons1[i in hData.hList, t in 2:T], rhpm[i,t] <= fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmCons2[i in hData.hList, t in 2:T], rhpm[i,t] >= -fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmCons3[i in hData.hList, t in 2:T], rhpm[i,t] <= λhp[i,t] + 2 * (1 - u_hm[i,t]) * fData.cz);
    @constraint(sprob, rhpmCons4[i in hData.hList, t in 2:T], rhpm[i,t] >= λhp[i,t] - 2 * (1 - u_hm[i,t]) * fData.cz);

    @constraint(sprob, rhppzCons11[i in hData.hList, t in 2:T], rhppz1[i,t] <= fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons12[i in hData.hList, t in 2:T], rhppz1[i,t] >= -fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons13[i in hData.hList, t in 2:T], rhppz1[i,t] <= λhzp1[i,t] + 2 * (1 - u_hp[i,t]) * fData.cz);
    @constraint(sprob, rhppzCons14[i in hData.hList, t in 2:T], rhppz1[i,t] >= λhzp1[i,t] - 2 * (1 - u_hp[i,t]) * fData.cz);

    @constraint(sprob, rhpmzCons11[i in hData.hList, t in 2:T], rhpmz1[i,t] <= fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons12[i in hData.hList, t in 2:T], rhpmz1[i,t] >= -fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons13[i in hData.hList, t in 2:T], rhpmz1[i,t] <= λhzm1[i,t] + 2 * (1 - u_hm[i,t]) * fData.cz);
    @constraint(sprob, rhpmzCons14[i in hData.hList, t in 2:T], rhpmz1[i,t] >= λhzm1[i,t] - 2 * (1 - u_hm[i,t]) * fData.cz);

    @constraint(sprob, rhppzCons21[i in hData.hList, t in 2:T], rhppz2[i,t] <= fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons22[i in hData.hList, t in 2:T], rhppz2[i,t] >= -fData.cz * u_hp[i,t]);
    @constraint(sprob, rhppzCons23[i in hData.hList, t in 2:T], rhppz2[i,t] <= λhzp3[i,t] + 2 * (1 - u_hp[i,t]) * fData.cz);
    @constraint(sprob, rhppzCons24[i in hData.hList, t in 2:T], rhppz2[i,t] >= λhzp3[i,t] - 2 * (1 - u_hp[i,t]) * fData.cz);

    @constraint(sprob, rhpmzCons21[i in hData.hList, t in 2:T], rhpmz2[i,t] <= fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons22[i in hData.hList, t in 2:T], rhpmz2[i,t] >= -fData.cz * u_hm[i,t]);
    @constraint(sprob, rhpmzCons23[i in hData.hList, t in 2:T], rhpmz2[i,t] <= λhzm3[i,t] + 2 * (1 - u_hm[i,t]) * fData.cz);
    @constraint(sprob, rhpmzCons24[i in hData.hList, t in 2:T], rhpmz2[i,t] >= λhzm3[i,t] - 2 * (1 - u_hm[i,t]) * fData.cz);

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
    @constraint(sprob, I_varConstr[i in fData.IDList, t in 2:(T-1)], λbat_lim[i,t] + λbat_trans[i,t] - λbat_trans[i,t+1] == 0);
    @constraint(sprob, I_varConstrT[i in fData.IDList], λbat_lim[i,T] + λbat_trans[i,T] == 0);
    @constraint(sprob, I_varConstr1[i in fData.IDList], λbat_lim[i,1] - λbat_trans[i,2] == 0);

    # set up constraints for demand/renewable variables
    @constraint(sprob, dp_varConstr[i in fData.IDList, t in 2:T], λpi[i,t] + λdp[i,t] == 0);
    @constraint(sprob, dq_varConstr[i in fData.IDList, t in 2:T], λqi[i,t] + λdq[i,t] == 0);
    @constraint(sprob, h_varConstr[i in fData.IDList, t in 2:T], -λpi[i,t] + λhp[i,t] >= 0);
    @constraint(sprob, uzhpConstr[i in hData.hList, t in 2:T], λhzp1[i,t] + λhzp2[i,t] + λhzp3[i,t] - expansion_factor * uData[i].RESPmax * λhp[i,t] >= 0);
    @constraint(sprob, uzhmConstr[i in hData.hList, t in 2:T], λhzm1[i,t] + λhzm2[i,t] + λhzm3[i,t] + expansion_factor * uData[i].RESPmin * λhp[i,t] >= 0);

    # set up constraints for x/y/z
    @constraint(sprob, xConstr1[k in fData.brList, t in 2:T; fData.rateA[k] < Inf], M[k] * (λptrans1[k,t] - λptrans2[k,t] + λqtrans1[k,t] - λqtrans2[k,t])
                 - fData.rateA[k] * ν1[k,t] + (λθ1[k,t] - λθ2[k,t])*2*pi + λxIni[k,t] == 0);
    @constraint(sprob, xConstr2[k in fData.brList, t in 2:T; fData.rateA[k] == Inf], M[k] * (λptrans1[k,t] - λptrans2[k,t] + λqtrans1[k,t] - λqtrans2[k,t])
                 + (λθ1[k,t] - λθ2[k,t])*2*pi + λxIni[k,t] == 0);
    @constraint(sprob, yConstr[i in fData.IDList], -sum(λbat_lim[i,t] + νbat_out[i,t]*bData.uCap[i] for t in 2:T) - λbat_lim[i,1] + λyIni[i] == 0);
    @constraint(sprob, zConstr[i in hData.hList], sum(-λhzp2[i,t] - λhzp3[i,t] - λhzm2[i,t] - λhzm3[i,t] for t in 2:T) + λzIni[i] - 
                λhp[i,t] * uData[i].RESP0[t] * expansion_factor == 0);

    # set up constraints for sp/sq
    @constraint(sprob, spConstr[i in fData.genIDList, t in 2:(T-1)], λspu[i,t] + λspl[i,t] + λrampu[i,t] - λrampu[i,t+1] + λrampl[i,t] - λrampl[i,t+1] - λpi[fData.Loc[i],t] ==
                -fData.cp[i].params[2]);
    @constraint(sprob, sqConstr[i in fData.genIDList, t in 2:T], λsqu[i,t] + λsql[i,t] - λqi[fData.Loc[i],t] == 0);
    @constraint(sprob, spConstr1[i in fData.genIDList], λspu[i,1] + λspl[i,1] + λspIni[i] - λrampu[i,2] - λrampl[i,2] == 0);
    @constraint(sprob, sqConstr1[i in fData.genIDList], λsqu[i,t] + λsql[i,t] + λsqIni[i] == 0);
    @constraint(sprob, spConstrT[i in fData.genIDList], λspu[i,T] + λspl[i,T] + λrampu[i,T] + λrampl[i,T] - λpi[fData.Loc[i],T] == -fData.cp[i].params[2]);

    # set up the objective function
    @objective(sprob, Max, -sum(sum(μ4[i,2,t]/4 + ν4[i,t]/4 + λvu[i,t]*vmaxT[t][i] + λvl[i,t]*vminT[t][i] - vmaxT[t][i]*vminT[t][i]*λv[i,t] 
        + sum(bData.ηβ[i][l] * λbat_eff[i,l,t] for l in 1:length(bData.ηα[i])) + λdp[i,t]*uData[i].DP0[t] + λdq[i,t]*uData[i].DQ0[t]
        + (uData[i].DPmax[t] - uData[i].DP0[t])*rdpp[i,t] + (uData[i].DPmin[t] - uData[i].DP0[t])*rdpm[i,t] 
        + (uData[i].DQmax[t] - uData[i].DQ0[t])*rdqp[i,t] + (uData[i].DQmin[t] - uData[i].DQ0[t])*rdqm[i,t] 
            for t in 2:T) + yhat[i]*λyIni[i] for i in fData.IDList) - 
        sum(sum(fData.Pmax[i] * λspu[i,t] + fData.Pmin[i] * λspl[i,t] + fData.Qmax[i] * λsqu[i,t] + fData.Qmin[i] * λsql[i,t] 
        + λrampu[i,t] * fData.RU[i] + λrampl[i,t] * fData.RD[i]
            for t in 2:T) + sphat[i] * λspIni[i] + sqhat[i] * λsqIni[i] for i in fData.genIDList) -
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
        sum(sum(λhp[i,t]*uData[i].RESP0[t] - λhzp3[i,t] - λhzm3[i,t] + rhpp[i,t]*uData[i].RESPmax - rhpm[i,t]*uData[i].RESPmin
        + rhppz1[i,t] + rhppz2[i,t] + rhpmz1[i,t] + rhpmz2[i,t] for t in 2:T) + zhat[i]*λzIni[i] for i in hData.hList)
            );

    return sprob;
end