# create the first stage problem
function createFirst(fData,uData,hData,bData,T,vmaxT,vminT,θDmaxT,θDminT,xLimit,no_threads = 1)
    # first-stage model without any scenarios
    θu = Dict();
    for t in 1:T
        θu[t] = Dict();
        for k in fData.brList
            θu[t][k] = max(abs(θDmaxT[t][k]),abs(θDminT[t][k]));
        end
    end
    # scaling issue exists for this primal problem.
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 1, "Threads" => no_threads));
    set_string_names_on_creation(mp, false);

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

    # set up the first stage formulation
    # violation slack variables
    @variable(mp,lpplus[i in fData.IDList] >= 0);
    @variable(mp,lpminus[i in fData.IDList] >= 0);
    @variable(mp,lqplus[i in fData.IDList] >= 0);
    @variable(mp,lqminus[i in fData.IDList] >= 0);

    # OPF variables
    @variable(mp,p[k in fData.brList]);
    @variable(mp,q[k in fData.brList]);
    @variable(mp,vminT[1][i] <= v[i in fData.IDList] <= vmaxT[1][i]);
    @variable(mp,vhat[i in fData.IDList]);
    @variable(mp,θ[i in fData.IDList]);

    # switching, battery setup, renewable expansion binary variables
    @variable(mp,x[k in fData.brList,t in 1:T],Bin);
    @variable(mp,y[i in fData.IDList],Bin);
    @variable(mp,z[i in hData.hList],Bin);
    @variable(mp,V >= 0);

    # set up the reference bus angle = 0
    for i in fData.IDList
        if fData.bType[i] == 3
        @constraint(mp,θ[i] == 0);
        end
    end

    # set up the variables: active/reactive power injection
    @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList] <= fData.Pmax[i]);
    @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList] <= fData.Qmax[i]);

    sphatsum = Dict();
    for i in fData.IDList
        sphatsum[i] = @expression(mp,0.0);
        if i in keys(fData.LocRev)
            for j in fData.LocRev[i]
            sphatsum[i] += sp[j];
            end
        end
    end

    sqhatsum = Dict();
    for i in fData.IDList
        sqhatsum[i] = @expression(mp,0.0);
        if i in keys(fData.LocRev)
            for j in fData.LocRev[i]
            sqhatsum[i] += sq[j];
            end
        end
    end

    @variable(mp,tAux2[k in fData.brList]);
    @variable(mp,tAux3[k in fData.brList]);
    @variable(mp,tAux4[k in fData.brList] >= 0);
    @variable(mp,tAux5[i in fData.IDList]);
    @variable(mp,tAux6[i in fData.IDList] >= 0);

    @variable(mp,vv[k in fData.brList]);
    @variable(mp,cos(θu[1][k]) <= cs[k in fData.brList] <= 1);
    @variable(mp,-sin(θu[1][k]) <= ss[k in fData.brList] <= sin(θu[1][k]));
    # add the strengthened bounds on cs
    csmax = Dict();
    csmin = Dict();
    ssmax = Dict();
    ssmin = Dict();
    for k in fData.brList
        if (θDmaxT[1][k] >= 0)&(θDminT[1][k] >= 0)
            csmax[k] = cos(θDminT[1][k]);
            csmin[k] = cos(θDmaxT[1][k]);
            @constraint(mp, cs[k] <= csmax[k]);
            @constraint(mp, cs[k] >= csmin[k]);
        elseif (θDmaxT[1][k] < 0)&(θDminT[1][k] < 0)
            csmax[k] = cos(θDmaxT[1][k]);
            csmin[k] = cos(θDminT[1][k]);
            @constraint(mp, cs[k] <= csmax[k]);
            @constraint(mp, cs[k] >= csmin[k]);
        else
            csmax[k] = 1;
            csmin[k] = min(cos(θDmaxT[1][k]),cos(θDminT[1][k]));
            @constraint(mp, cs[k] <= csmax[k]);
            @constraint(mp, cs[k] >= csmin[k]);
        end
        @constraint(mp, cs[k] >= (cos(θDmaxT[1][k]) - cos(θDminT[1][k]))/(θDmaxT[1][k] - θDminT[1][k])*((θ[k[1]] - fData.σ[k]) - θ[k[2]] - θDminT[1][k]) + cos(θDminT[1][k]));
    end
    # add the strengthened bounds on ss
    for k in fData.brList
        ssmax[k] = sin(θDmaxT[1][k]);
        ssmin[k] = sin(θDminT[1][k]);
        @constraint(mp, ss[k] <= ssmax[k]);
        @constraint(mp, ss[k] >= ssmin[k]);
    end

    @variable(mp,wc[k in fData.brList]);
    @variable(mp,ws[k in fData.brList]);
    @constraint(mp,wcEquality[k in fData.brList;k[1] < k[2]], wc[k] == wc[(k[2],k[1],k[3])]);
    @constraint(mp,wsEquality[k in fData.brList;k[1] < k[2]], ws[k] == -ws[(k[2],k[1],k[3])]);
    @constraint(mp,vvEquality[k in fData.brList;k[1] < k[2]], vv[k] == vv[(k[2],k[1],k[3])]);
    @constraint(mp,csEquality[k in fData.brList;k[1] < k[2]], cs[k] == cs[(k[2],k[1],k[3])]);
    @constraint(mp,ssEquality[k in fData.brList;k[1] < k[2]], ss[k] == -ss[(k[2],k[1],k[3])]);

    @constraint(mp,lineConstrP[k in fData.brList],p[k] == (fData.g[k]*vhat[k[1]]/(fData.τ1[k]^2) - fData.g[k]*wc[k] - fData.b[k]*ws[k]));
    @constraint(mp,lineConstrQ[k in fData.brList],q[k] == ((-fData.b[k] - fData.bc[k]/2)*vhat[k[1]]/(fData.τ1[k]^2) + fData.b[k]*wc[k] - fData.g[k]*ws[k]));
    socList1 = Dict();
    socList2 = Dict();
    socList3 = Dict();
    socList4 = Dict();
    socList5 = Dict();
    for k in fData.brList
        if fData.rateA[k] < Inf
            socList1[k] = [p[k],q[k]];
        end
            socList2[k] = [vv[k],vhat[k[1]]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2]]/((fData.τ2[k]^2)*sqrt(2))];
            socList3[k] = [tAux2[k],tAux3[k]];
            socList5[k] = [wc[k],ws[k],vhat[k[1]]/((fData.τ1[k]^2)*sqrt(2)),vhat[k[2]]/((fData.τ2[k]^2)*sqrt(2))];
    end
    for i in fData.IDList
        socList4[i] = [v[i],tAux5[i]];
    end

    @constraint(mp,socConstraint1[k in fData.brList; fData.rateA[k] < Inf],[fData.rateA[k]; socList1[k]] in SecondOrderCone());
    @constraint(mp,socConstraint2[k in fData.brList], [(vhat[k[1]]/(fData.τ1[k]^2) + vhat[k[2]]/(fData.τ2[k]^2))/sqrt(2); socList2[k]] in SecondOrderCone());
    @constraint(mp,socConstraint3[k in fData.brList],[tAux4[k]; socList3[k]] in SecondOrderCone());
    @constraint(mp,socConstraint4[i in fData.IDList],[tAux6[i]; socList4[i]] in SecondOrderCone());
    @constraint(mp,socConstraint5[k in fData.brList],[(vhat[k[1]]/(fData.τ1[k]^2) + vhat[k[2]]/(fData.τ2[k]^2))/sqrt(2); socList5[k]] in SecondOrderCone());

    @constraint(mp,auxConstr2[k in fData.brList],sqrt((1-cos(θu[1][k]))/(θu[1][k])^2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]]) == tAux2[k]);
    @constraint(mp,auxConstr3[k in fData.brList],cs[k] - 3/4 == tAux3[k]);
    @constraint(mp,auxConstr4[k in fData.brList],5/4 - cs[k] == tAux4[k]);
    @constraint(mp,sinMock1[k in fData.brList],
                ss[k] <= cos(θu[1][k]/2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]] - θu[1][k]/2) + sin(θu[1][k]/2));
    @constraint(mp,sinMock2[k in fData.brList],
                ss[k] >= cos(θu[1][k]/2)*((θ[k[1]] - fData.σ[k]) - θ[k[2]] + θu[1][k]/2) - sin(θu[1][k]/2));
    @constraint(mp,angleDiff1[k in fData.brList],(θ[k[1]] - fData.σ[k]) - θ[k[2]] <= θDmaxT[1][k]);
    @constraint(mp,angleDiff2[k in fData.brList],(θ[k[1]] - fData.σ[k]) - θ[k[2]] >= θDminT[1][k]);
    @constraint(mp,auxConstr5[i in fData.IDList],tAux5[i] == vhat[i] - 1/4);
    @constraint(mp,auxConstr6[i in fData.IDList],tAux6[i] == vhat[i] + 1/4);
    @constraint(mp,v2Mock2[i in fData.IDList], vhat[i] - (vmaxT[1][i] + vminT[1][i])*v[i] <= -vmaxT[1][i]*vminT[1][i]);
    @constraint(mp,vvMock1[k in fData.brList],
                vv[k] >= vminT[1][k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vminT[1][k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k]));
    @constraint(mp,vvMock2[k in fData.brList],
                vv[k] >= vmaxT[1][k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmaxT[1][k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k]));
    @constraint(mp,vvMock3[k in fData.brList],
                vv[k] <= vminT[1][k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vmaxT[1][k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vminT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k]));
    @constraint(mp,vvMock4[k in fData.brList],
                vv[k] <= vmaxT[1][k[1]]*v[k[2]]/(fData.τ1[k]*fData.τ2[k]) + vminT[1][k[2]]*v[k[1]]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k]));

    @constraint(mp,wcMock1[k in fData.brList],
                wc[k] >= vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + csmin[k]*vv[k] - vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);
    @constraint(mp,wcMock2[k in fData.brList],
                wc[k] >= vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + csmax[k]*vv[k] - vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
    @constraint(mp,wcMock3[k in fData.brList],
                wc[k] <= vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + vv[k]*csmax[k] - vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k]);
    @constraint(mp,wcMock4[k in fData.brList],
                wc[k] <= vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k] + vv[k]*csmin[k] - vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k]);

    @constraint(mp,wsMock1[k in fData.brList],
                ws[k] >= vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + ssmin[k]*vv[k] - vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);
    @constraint(mp,wsMock2[k in fData.brList],
                ws[k] >= vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + ssmax[k]*vv[k] - vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
    @constraint(mp,wsMock3[k in fData.brList],
                ws[k] <= vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + vv[k]*ssmax[k] - vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k]);
    @constraint(mp,wsMock4[k in fData.brList],
                ws[k] <= vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k] + vv[k]*ssmin[k] - vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k]);

    @constraint(mp,tanConstr1[k in fData.brList],
                    ws[k] - tan(θDmaxT[1][k])*wc[k] <= 0);
    @constraint(mp,tanConstr2[k in fData.brList],
                    ws[k] - tan(θDminT[1][k])*wc[k] >= 0);

    vδ = Dict();
    θϕ = Dict();
    θδ = Dict();
    for i in fData.IDList
        vδ[i] = vmaxT[1][i]+vminT[1][i];
    end
    for k in fData.brList
        θϕ[k] = (θDmaxT[1][k] + θDminT[1][k])/2;
        θδ[k] = (θDmaxT[1][k] - θDminT[1][k])/2;
    end

    @constraint(mp,lncConstr1[k in fData.brList],
                vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k]*cos(θϕ[k]) + ws[k]*sin(θϕ[k])) -
                vmaxT[1][k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1]]/(fData.τ1[k]^2) -
                vmaxT[1][k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2]]/(fData.τ2[k]^2) >=
                vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k]) -
                vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])));
    @constraint(mp,lncConstr2[k in fData.brList],
                vδ[k[1]]*vδ[k[2]]/(fData.τ1[k]*fData.τ2[k])*(wc[k]*cos(θϕ[k]) + ws[k]*sin(θϕ[k])) -
                vminT[1][k[2]]/fData.τ2[k]*cos(θδ[k])*vδ[k[2]]/fData.τ2[k]*vhat[k[1]]/(fData.τ1[k]^2) -
                vminT[1][k[1]]/fData.τ1[k]*cos(θδ[k])*vδ[k[1]]/fData.τ1[k]*vhat[k[2]]/(fData.τ2[k]^2) >=
                -vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k])*(vminT[1][k[1]]*vminT[1][k[2]]/(fData.τ1[k]*fData.τ2[k]) -
                vmaxT[1][k[1]]*vmaxT[1][k[2]]/(fData.τ1[k]*fData.τ2[k])));

    # power flow balance
    @constraint(mp,totalP[i in fData.IDList], sum(p[k] for k in branchDict1[i]) + vhat[i]*fData.gs[i]
                + lpplus[i] - lpminus[i] == sphatsum[i] + uData[i].RESP0[1] - uData[i].DP0[1]);
    @constraint(mp,totalQ[i in fData.IDList], sum(q[k] for k in branchDict1[i]) - vhat[i]*fData.bs[i]
                + lqplus[i] - lqminus[i] == sqhatsum[i] - uData[i].DQ0[1]);

    # add the battery/switching/renewable expansion
    @constraint(mp,switching[k in fData.brList,t in 2:T], x[k,t] <= x[k,t-1]);
    @constraint(mp,switching0[k in fData.brList], x[k,1] == 1);
    @constraint(mp,switching_rev[k in fData.brList,t in 1:T;k[1] < k[2]], x[k,t] == x[(k[2],k[1],k[3]),t]);
    @constraint(mp,switching_lim, sum(x[k,T] for k in fData.brList if k[1] < k[2]) >= length(fData.brList)/2 - xLimit);

    @objective(mp, Min, fData.cz*sum(lpplus[i] + lpminus[i] + lqplus[i] + lqminus[i] for i in fData.IDList) + 
        sum(fData.cp[i].params[2]*sp[i] for i in fData.genIDList) + V + sum(bData.cost[i]*y[i] for i in fData.IDList) + 
        sum(hData.cost[i] * z[i] for i in hData.hList));

    return mp;
end

# create the second-stage problem
function createSecond(fData,uData,hData,bData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uDict)
    # first-stage model without any scenarios
    θu = Dict();
    for t in 1:T
        θu[t] = Dict();
        for k in fData.brList
            θu[t][k] = max(abs(θDmaxT[t][k]),abs(θDminT[t][k]));
        end
    end
    # scaling issue exists for this primal problem.
    subp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 0, "Threads" => 1));
    # subp = Model(COPT.Optimizer);

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

    # set up the first stage formulation
    # violation slack variables
    @variable(subp,lpplus[i in fData.IDList, t in 2:T] >= 0);
    @variable(subp,lpminus[i in fData.IDList, t in 2:T] >= 0);
    @variable(subp,lqplus[i in fData.IDList, t in 2:T] >= 0);
    @variable(subp,lqminus[i in fData.IDList, t in 2:T] >= 0);

    # OPF variables
    @variable(subp,p[k in fData.brList, t in 2:T]);
    @variable(subp,q[k in fData.brList, t in 2:T]);
    @variable(subp,vminT[t][i] <= v[i in fData.IDList, t in 2:T] <= vmaxT[t][i]);
    @variable(subp,vhat[i in fData.IDList, t in 2:T]);
    @variable(subp,θ[i in fData.IDList, t in 2:T]);

    # switching, battery setup, renewable expansion binary variables
    @variable(subp,x[k in fData.brList,t in 2:T]);
    @variable(subp,y[i in fData.IDList]);
    @variable(subp,z[i in hData.hList]);
    
    # set up the reference bus angle = 0
    for i in fData.IDList
        for t in 2:T
            if fData.bType[i] == 3
                @constraint(subp,θ[i,t] == 0);
            end
        end
    end

    # set up the variables: active/reactive power injection
    @variable(subp, fData.Pmin[i] <= sp[i in fData.genIDList, t in 1:T] <= fData.Pmax[i]);
    @variable(subp, fData.Qmin[i] <= sq[i in fData.genIDList, t in 1:T] <= fData.Qmax[i]);
    @constraint(subp, spIni[i in fData.genIDList], sp[i,1] == sphat[i]);
    @constraint(subp, sqIni[i in fData.genIDList], sq[i,1] == sqhat[i]);
    # ramping constraints
    @constraint(subp,ramp_up[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t-1] <= fData.RU[i]);
    @constraint(subp,ramp_dn[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t-1] >= fData.RD[i]);

    sphatsum = Dict();
    sqhatsum = Dict();
    for t in 2:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(subp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end

        for i in fData.IDList
            sqhatsum[i,t] = @expression(subp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(subp,tAux2[k in fData.brList, t in 2:T]);
    @variable(subp,tAux3[k in fData.brList, t in 2:T]);
    @variable(subp,tAux4[k in fData.brList, t in 2:T] >= 0);
    @variable(subp,tAux5[i in fData.IDList, t in 2:T]);
    @variable(subp,tAux6[i in fData.IDList, t in 2:T] >= 0);

    @variable(subp,vv[k in fData.brList, t in 2:T]);
    @variable(subp,cos(θu[t][k]) <= cs[k in fData.brList, t in 2:T] <= 1);
    @variable(subp,-sin(θu[t][k]) <= ss[k in fData.brList, t in 2:T] <= sin(θu[t][k]));
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
                @constraint(subp, cs[k,t] <= csmax[k,t]);
                @constraint(subp, cs[k,t] >= csmin[k,t]);
            elseif (θDmaxT[t][k] < 0)&(θDminT[t][k] < 0)
                csmax[k,t] = cos(θDmaxT[t][k]);
                csmin[k,t] = cos(θDminT[t][k]);
                @constraint(subp, cs[k,t] <= csmax[k,t]);
                @constraint(subp, cs[k,t] >= csmin[k,t]);
            else
                csmax[k,t] = 1;
                csmin[k,t] = min(cos(θDmaxT[t][k]),cos(θDminT[t][k]));
                @constraint(subp, cs[k,t] <= csmax[k,t]);
                @constraint(subp, cs[k,t] >= csmin[k,t]);
            end
            @constraint(subp, cs[k,t] >= (cos(θDmaxT[t][k]) - cos(θDminT[t][k]))/(θDmaxT[t][k] - θDminT[t][k])*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t] - θDminT[t][k]) + cos(θDminT[t][k]));
        end
        # add the strengthened bounds on ss
        for k in fData.brList
            ssmax[k,t] = sin(θDmaxT[t][k]);
            ssmin[k,t] = sin(θDminT[t][k]);
            @constraint(subp, ss[k,t] <= ssmax[k,t]);
            @constraint(subp, ss[k,t] >= ssmin[k,t]);
        end
    end

    @variable(subp,wc[k in fData.brList, t in 2:T]);
    @variable(subp,ws[k in fData.brList, t in 2:T]);
    @constraint(subp,wcEquality[k in fData.brList, t in 2:T;k[1] < k[2]], wc[k,t] == wc[(k[2],k[1],k[3]),t]);
    @constraint(subp,wsEquality[k in fData.brList, t in 2:T;k[1] < k[2]], ws[k,t] == -ws[(k[2],k[1],k[3]),t]);
    @constraint(subp,vvEquality[k in fData.brList, t in 2:T;k[1] < k[2]], vv[k,t] == vv[(k[2],k[1],k[3]),t]);
    @constraint(subp,csEquality[k in fData.brList, t in 2:T;k[1] < k[2]], cs[k,t] == cs[(k[2],k[1],k[3]),t]);
    @constraint(subp,ssEquality[k in fData.brList, t in 2:T;k[1] < k[2]], ss[k,t] == -ss[(k[2],k[1],k[3]),t]);

    M = Dict();
    for k in fData.brList
        if fData.rateA[k] < Inf
            M[k] = fData.rateA[k]^2;
        else
            M[k] = 10000;
        end
    end
    @constraint(subp,lineConstrP1[k in fData.brList, t in 2:T],p[k,t] <= (fData.g[k]*vhat[k[1],t]/(fData.τ1[k]^2) - fData.g[k]*wc[k,t] - fData.b[k]*ws[k,t]) + M[k]*(1 - x[k,t]));
    @constraint(subp,lineConstrP2[k in fData.brList, t in 2:T],p[k,t] >= (fData.g[k]*vhat[k[1],t]/(fData.τ1[k]^2) - fData.g[k]*wc[k,t] - fData.b[k]*ws[k,t]) - M[k]*(1 - x[k,t]));
    @constraint(subp,lineConstrQ1[k in fData.brList, t in 2:T],q[k,t] <= ((-fData.b[k] - fData.bc[k]/2)*vhat[k[1],t]/(fData.τ1[k]^2) + fData.b[k]*wc[k,t] - fData.g[k]*ws[k,t]) + M[k]*(1 - x[k,t]));
    @constraint(subp,lineConstrQ2[k in fData.brList, t in 2:T],q[k,t] >= ((-fData.b[k] - fData.bc[k]/2)*vhat[k[1],t]/(fData.τ1[k]^2) + fData.b[k]*wc[k,t] - fData.g[k]*ws[k,t]) - M[k]*(1 - x[k,t]));
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

    @constraint(subp,socConstraint1[k in fData.brList, t in 2:T; fData.rateA[k] < Inf],[fData.rateA[k] * x[k,t]; socList1[k,t]] in SecondOrderCone());
    @constraint(subp,socConstraint2[k in fData.brList, t in 2:T], [(vhat[k[1],t]/(fData.τ1[k]^2) + vhat[k[2],t]/(fData.τ2[k]^2))/sqrt(2); socList2[k,t]] in SecondOrderCone());
    @constraint(subp,socConstraint3[k in fData.brList, t in 2:T],[tAux4[k,t]; socList3[k,t]] in SecondOrderCone());
    @constraint(subp,socConstraint4[i in fData.IDList, t in 2:T],[tAux6[i,t]; socList4[i,t]] in SecondOrderCone());
    @constraint(subp,socConstraint5[k in fData.brList, t in 2:T],[(vhat[k[1],t]/(fData.τ1[k]^2) + vhat[k[2],t]/(fData.τ2[k]^2))/sqrt(2); socList5[k,t]] in SecondOrderCone());

    @constraint(subp,auxConstr2[k in fData.brList, t in 2:T],sqrt((1-cos(θu[t][k]))/(θu[t][k])^2)*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t]) == tAux2[k,t]);
    @constraint(subp,auxConstr3[k in fData.brList, t in 2:T],cs[k,t] - 3/4 == tAux3[k,t]);
    @constraint(subp,auxConstr4[k in fData.brList, t in 2:T],5/4 - cs[k,t] == tAux4[k,t]);
    @constraint(subp,sinMock1[k in fData.brList, t in 2:T],
                ss[k,t] <= cos(θu[t][k]/2)*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t] - θu[t][k]/2) + sin(θu[t][k]/2));
    @constraint(subp,sinMock2[k in fData.brList, t in 2:T],
                ss[k,t] >= cos(θu[t][k]/2)*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t] + θu[t][k]/2) - sin(θu[t][k]/2));
    @constraint(subp,angleDiff1[k in fData.brList, t in 2:T],(θ[k[1],t] - fData.σ[k]) - θ[k[2],t] <= θDmaxT[t][k] + 2*pi*(1 - x[k,t]));
    @constraint(subp,angleDiff2[k in fData.brList, t in 2:T],(θ[k[1],t] - fData.σ[k]) - θ[k[2],t] >= θDminT[t][k] - 2*pi*(1 - x[k,t]));
    @constraint(subp,auxConstr5[i in fData.IDList, t in 2:T],tAux5[i,t] == vhat[i,t] - 1/4);
    @constraint(subp,auxConstr6[i in fData.IDList, t in 2:T],tAux6[i,t] == vhat[i,t] + 1/4);
    @constraint(subp,v2Mock2[i in fData.IDList, t in 2:T], vhat[i,t] - (vmaxT[t][i] + vminT[t][i])*v[i,t] <= -vmaxT[t][i]*vminT[t][i]);
    @constraint(subp,vvMock1[k in fData.brList, t in 2:T],
                vv[k,t] >= vminT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vminT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));
    @constraint(subp,vvMock2[k in fData.brList, t in 2:T],
                vv[k,t] >= vmaxT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vmaxT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));
    @constraint(subp,vvMock3[k in fData.brList, t in 2:T],
                vv[k,t] <= vminT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vmaxT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vminT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));
    @constraint(subp,vvMock4[k in fData.brList, t in 2:T],
                vv[k,t] <= vmaxT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vminT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));

    @constraint(subp,wcMock1[k in fData.brList, t in 2:T],
                wc[k,t] >= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + csmin[k,t]*vv[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k,t]);
    @constraint(subp,wcMock2[k in fData.brList, t in 2:T],
                wc[k,t] >= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + csmax[k,t]*vv[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k,t]);
    @constraint(subp,wcMock3[k in fData.brList, t in 2:T],
                wc[k,t] <= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + vv[k,t]*csmax[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k,t]);
    @constraint(subp,wcMock4[k in fData.brList, t in 2:T],
                wc[k,t] <= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + vv[k,t]*csmin[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k,t]);

    @constraint(subp,wsMock1[k in fData.brList, t in 2:T],
                ws[k,t] >= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + ssmin[k,t]*vv[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k,t]);
    @constraint(subp,wsMock2[k in fData.brList, t in 2:T],
                ws[k,t] >= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + ssmax[k,t]*vv[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k,t]);
    @constraint(subp,wsMock3[k in fData.brList, t in 2:T],
                ws[k,t] <= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + vv[k,t]*ssmax[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k,t]);
    @constraint(subp,wsMock4[k in fData.brList, t in 2:T],
                ws[k,t] <= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + vv[k,t]*ssmin[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k,t]);

    @constraint(subp,tanConstr1[k in fData.brList, t in 2:T],
                    ws[k,t] - tan(θDmaxT[t][k])*wc[k,t] <= 0);
    @constraint(subp,tanConstr2[k in fData.brList, t in 2:T],
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

    @constraint(subp,lncConstr1[k in fData.brList, t in 2:T],
                vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*(wc[k,t]*cos(θϕ[k,t]) + ws[k,t]*sin(θϕ[k,t])) -
                vmaxT[t][k[2]]/fData.τ2[k]*cos(θδ[k,t])*vδ[k[2],t]/fData.τ2[k]*vhat[k[1],t]/(fData.τ1[k]^2) -
                vmaxT[t][k[1]]/fData.τ1[k]*cos(θδ[k,t])*vδ[k[1],t]/fData.τ1[k]*vhat[k[2],t]/(fData.τ2[k]^2) >=
                vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]) -
                vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])));
    @constraint(subp,lncConstr2[k in fData.brList, t in 2:T],
                vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*(wc[k,t]*cos(θϕ[k,t]) + ws[k,t]*sin(θϕ[k,t])) -
                vminT[t][k[2]]/fData.τ2[k]*cos(θδ[k,t])*vδ[k[2],t]/fData.τ2[k]*vhat[k[1],t]/(fData.τ1[k]^2) -
                vminT[t][k[1]]/fData.τ1[k]*cos(θδ[k,t])*vδ[k[1],t]/fData.τ1[k]*vhat[k[2],t]/(fData.τ2[k]^2) >=
                -vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]) -
                vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])));

    # battery dynamics
    @variable(subp, ep_var[i in fData.IDList, t in 2:T]);
    @variable(subp, eq_var[i in fData.IDList, t in 2:T]);
    @variable(subp, f_var[i in fData.IDList, t in 2:T]);
    @variable(subp, I_var[i in fData.IDList, t in 1:T] >= 0);
    @constraint(subp, soc_limit[i in fData.IDList, t in 1:T], I_var[i,t] <= bData.cap[i]*y[i]);
    @constraint(subp, soc_trans[i in fData.IDList, t in 2:T], I_var[i,t] == I_var[i,t-1] - fData.Δt * f_var[i,t]);
    @constraint(subp, bat_out[i in fData.IDList, t in 2:T], [bData.uCap[i]*y[i]; [ep_var[i,t], eq_var[i,t]]] in SecondOrderCone());
    @constraint(subp, bat_eff[i in fData.IDList, l in 1:length(bData.ηα[i]), t in 2:T], ep_var[i,t] - bData.ηα[i][l]*f_var[i,t] <= bData.ηβ[i][l]);

    # power flow balance
    @variable(subp, dp_var[i in fData.IDList, t in 2:T]);
    @variable(subp, dq_var[i in fData.IDList, t in 2:T]);
    @variable(subp, h_var[i in fData.IDList, t in 2:T] >= 0);
    @constraint(subp,totalP[i in fData.IDList, t in 2:T], sum(p[k,t] for k in branchDict1[i]) + vhat[i,t]*fData.gs[i]
                + lpplus[i,t] - lpminus[i,t] - h_var[i,t] + dp_var[i,t] - ep_var[i,t] == sphatsum[i,t]);
    @constraint(subp,totalQ[i in fData.IDList, t in 2:T], sum(q[k,t] for k in branchDict1[i]) - vhat[i,t]*fData.bs[i]
                + lqplus[i,t] - lqminus[i,t]  + dq_var[i,t] - eq_var[i,t] == sqhatsum[i,t]);

    # d_var constraints: demand uncertainty
    groupList = 1:length(groupDict[2]);
    group_rev = Dict();
    for i in fData.IDList
        for m in groupList
            if i in groupDict[2][m]
                group_rev[i] = m;
            end
        end
    end
    # @variable(subp, u_dp[m in groupList, t in 2:T], Bin);
    # @variable(subp, u_dm[m in groupList, t in 2:T], Bin);
    # @constraint(subp, uncertain_budget_d, sum(u_dp[m,t] + u_dm[m,t] for m in groupList for t in 2:T) <= Γ["d"]);
    # @constraint(subp, one_extreme_pt_d[m in groupList, t in 2:T],u_dp[m,t] + u_dm[m,t] <= 1);
    u_dp = uDict["u_dp"];
    u_dm = uDict["u_dm"];
    @constraint(subp, dp_cal[i in fData.IDList, t in 2:T], dp_var[i,t] == (uData[i].DPmax[t] - uData[i].DP0[t])*u_dp[group_rev[i],t] + 
        (uData[i].DPmin[t] - uData[i].DP0[t])*u_dm[group_rev[i],t] + uData[i].DP0[t]);
    @constraint(subp, dq_cal[i in fData.IDList, t in 2:T], dq_var[i,t] == (uData[i].DQmax[t] - uData[i].DQ0[t])*u_dp[group_rev[i],t] + 
        (uData[i].DQmin[t] - uData[i].DQ0[t])*u_dm[group_rev[i],t] + uData[i].DQ0[t]);

    # h_var constraints: renewable uncertainty
    # @variable(subp, u_hp[i in hData.hList, t in 2:T], Bin);
    # @variable(subp, u_hm[i in hData.hList, t in 2:T], Bin);
    # @variable(subp, 0 <= uz_hp[i in hData.hList, t in 2:T] <= 1);
    # @variable(subp, 0 <= uz_hm[i in hData.hList, t in 2:T] <= 1);
    # @constraint(subp, uncertain_budget_h, sum(u_hp[i,t] + u_hm[i,t] for i in hData.hList for t in 2:T) <= Γ["h"]);
    # @constraint(subp, one_extreme_pt_h[i in hData.hList, t in 2:T],u_hp[i,t] + u_hm[i,t] <= 1);
    # @constraint(subp, uzp_linear1[i in hData.hList, t in 2:T], uz_hp[i,t] <= u_hp[i,t]);
    # @constraint(subp, uzp_linear2[i in hData.hList, t in 2:T], uz_hp[i,t] <= z[i]);
    # @constraint(subp, uzp_linear3[i in hData.hList, t in 2:T], uz_hp[i,t] >= u_hp[i,t] + z[i] - 1);
    # @constraint(subp, uzm_linear1[i in hData.hList, t in 2:T], uz_hm[i,t] <= u_hm[i,t]);
    # @constraint(subp, uzm_linear2[i in hData.hList, t in 2:T], uz_hm[i,t] <= z[i]);
    # @constraint(subp, uzm_linear3[i in hData.hList, t in 2:T], uz_hm[i,t] >= u_hm[i,t] + z[i] - 1);
    u_hp = uDict["u_hp"];
    u_hm = uDict["u_hm"];
    uz_hp = Dict();
    uz_hm = Dict();
    for i in hData.hList
        for t in 2:T
            uz_hp[i,t] = z[i] * u_hp[i,t];
            uz_hm[i,t] = z[i] * u_hm[i,t];
        end
    end
    @constraint(subp, hp_cal1[i in fData.IDList, t in 2:T; i in hData.hList], h_var[i,t] <= uData[i].RESP0[t]*((1+expansion_factor * z[i]) + 
            uData[i].RESPmax*(u_hp[i,t]+expansion_factor*uz_hp[i,t]) - uData[i].RESPmin*(u_hm[i,t]+expansion_factor*uz_hm[i,t])));
    @constraint(subp, hp_cal2[i in fData.IDList, t in 2:T; !(i in hData.hList)], h_var[i,t] == 0.0);    

    # bind the previous solutions
    @constraint(subp,xIni[k in fData.brList, t in 2:T],x[k,t] == xhat[k,t]);
    @constraint(subp,yIni[i in fData.IDList],y[i] == yhat[i]);
    @constraint(subp,zIni[i in hData.hList],z[i] == zhat[i]);

    # objective value
    @objective(subp, Min, fData.cz*sum(lpplus[i,t] + lpminus[i,t] + lqplus[i,t] + lqminus[i,t] for i in fData.IDList for t in 2:T) + 
        sum(fData.cp[i].params[2]*sp[i,t] for i in fData.genIDList for t in 2:T));

    return subp;
end

# change the second-stage formulation with updated first-stage solution
function changeSecond(fData,hData,subp,xhat,yhat,zhat,sphat,sqhat)
    for t in 1:T
        if t == 1
            for i in fData.genIDList
                set_normalized_rhs(subp[:spIni][i], sphat[i]);
                set_normalized_rhs(subp[:sqIni][i], sqhat[i]);
            end
        else
            for k in fData.brList
                set_normalized_rhs(subp[:xIni][k,t], xhat[k,t]);
            end
            for i in fData.IDList
                set_normalized_rhs(subp[:yIni][i,t], yhat[i,t]);
            end
            for i in hData.hList
                set_normalized_rhs(subp[:zIni][i,t], zhat[i,t]);
            end
        end
    end
    return subp;
end

# append the scenarios to the first-stage problem, given the generated extreme point
function appendScen(mp,fData,uData,hData,T,groupDict,vmaxT,vminT,θDmaxT,θDminT,expansion_factor,uList)
    # first-stage model without any scenarios
    θu = Dict();
    for t in 1:T
        θu[t] = Dict();
        for k in fData.brList
            θu[t][k] = max(abs(θDmaxT[t][k]),abs(θDminT[t][k]));
        end
    end

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

    # obtain the big-M
    M = Dict();
    for k in fData.brList
        if fData.rateA[k] < Inf
            M[k] = fData.rateA[k] * 10;
        else
            M[k] = 10000;
        end
    end

    # obtain the demand groups
    groupList = 1:length(groupDict[2]);
    group_rev = Dict();
    for i in fData.IDList
        for m in groupList
            if i in groupDict[2][m]
                group_rev[i] = m;
            end
        end
    end

    for ω in eachindex(uList)
        # item: contains 4 elements (u_dp, u_dm, u_hp, u_hm)
        # add the scenarios

        # violation slack variables
        lpplus = @variable(mp,[i in fData.IDList, t in 2:T], lower_bound = 0);
        lpminus = @variable(mp,[i in fData.IDList, t in 2:T], lower_bound = 0);
        lqplus = @variable(mp,[i in fData.IDList, t in 2:T], lower_bound = 0);
        lqminus = @variable(mp,[i in fData.IDList, t in 2:T], lower_bound = 0);

        # OPF variables
        p = @variable(mp,[k in fData.brList, t in 2:T]);
        q = @variable(mp,[k in fData.brList, t in 2:T]);
        v = @variable(mp,[i in fData.IDList, t in 2:T], lower_bound = vminT[t][i], upper_bound = vmaxT[t][i]);
        vhat = @variable(mp,[i in fData.IDList, t in 2:T]);
        θ = @variable(mp,[i in fData.IDList, t in 2:T]);
    
        # set up the reference bus angle = 0
        for i in fData.IDList
            for t in 2:T
                if fData.bType[i] == 3
                    @constraint(mp,θ[i,t] == 0);
                end
            end
        end

        # set up the variables: active/reactive power injection
        sp = @variable(mp, [i in fData.genIDList, t in 1:T], lower_bound = fData.Pmin[i], upper_bound = fData.Pmax[i]);
        sq = @variable(mp, [i in fData.genIDList, t in 1:T], lower_bound = fData.Qmin[i], upper_bound = fData.Qmax[i]);
        @constraint(mp, [i in fData.genIDList], sp[i,1] == mp[:sp][i]);
        @constraint(mp, [i in fData.genIDList], sq[i,1] == mp[:sq][i]);

        sphatsum = Dict();
        sqhatsum = Dict();
        for t in 2:T
            for i in fData.IDList
                sphatsum[i,t] = @expression(mp,0.0);
                if i in keys(fData.LocRev)
                    for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                    end
                end
            end

            for i in fData.IDList
                sqhatsum[i,t] = @expression(mp,0.0);
                if i in keys(fData.LocRev)
                    for j in fData.LocRev[i]
                        sqhatsum[i,t] += sq[j,t];
                    end
                end
            end
        end

        tAux2 = @variable(mp,[k in fData.brList, t in 2:T]);
        tAux3 = @variable(mp,[k in fData.brList, t in 2:T]);
        tAux4 = @variable(mp,[k in fData.brList, t in 2:T], lower_bound = 0);
        tAux5 = @variable(mp,[i in fData.IDList, t in 2:T]);
        tAux6 = @variable(mp,[i in fData.IDList, t in 2:T], lower_bound = 0);

        vv = @variable(mp,[k in fData.brList, t in 2:T]);
        cs = @variable(mp,[k in fData.brList, t in 2:T], lower_bound = cos(θu[t][k]), upper_bound = 1);
        ss = @variable(mp,[k in fData.brList, t in 2:T], lower_bound = -sin(θu[t][k]), upper_bound = sin(θu[t][k]));
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
                    @constraint(mp, cs[k,t] <= csmax[k,t]);
                    @constraint(mp, cs[k,t] >= csmin[k,t]);
                elseif (θDmaxT[t][k] < 0)&(θDminT[t][k] < 0)
                    csmax[k,t] = cos(θDmaxT[t][k]);
                    csmin[k,t] = cos(θDminT[t][k]);
                    @constraint(mp, cs[k,t] <= csmax[k,t]);
                    @constraint(mp, cs[k,t] >= csmin[k,t]);
                else
                    csmax[k,t] = 1;
                    csmin[k,t] = min(cos(θDmaxT[t][k]),cos(θDminT[t][k]));
                    @constraint(mp, cs[k,t] <= csmax[k,t]);
                    @constraint(mp, cs[k,t] >= csmin[k,t]);
                end
                @constraint(mp, cs[k,t] >= (cos(θDmaxT[t][k]) - cos(θDminT[t][k]))/(θDmaxT[t][k] - θDminT[t][k])*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t] - θDminT[t][k]) + cos(θDminT[t][k]));
            end
            # add the strengthened bounds on ss
            for k in fData.brList
                ssmax[k,t] = sin(θDmaxT[t][k]);
                ssmin[k,t] = sin(θDminT[t][k]);
                @constraint(mp, ss[k,t] <= ssmax[k,t]);
                @constraint(mp, ss[k,t] >= ssmin[k,t]);
            end
        end

        wc = @variable(mp,[k in fData.brList, t in 2:T]);
        ws = @variable(mp,[k in fData.brList, t in 2:T]);
        @constraint(mp,[k in fData.brList, t in 2:T;k[1] < k[2]], wc[k,t] == wc[(k[2],k[1],k[3]),t]);
        @constraint(mp,[k in fData.brList, t in 2:T;k[1] < k[2]], ws[k,t] == -ws[(k[2],k[1],k[3]),t]);
        @constraint(mp,[k in fData.brList, t in 2:T;k[1] < k[2]], vv[k,t] == vv[(k[2],k[1],k[3]),t]);
        @constraint(mp,[k in fData.brList, t in 2:T;k[1] < k[2]], cs[k,t] == cs[(k[2],k[1],k[3]),t]);
        @constraint(mp,[k in fData.brList, t in 2:T;k[1] < k[2]], ss[k,t] == -ss[(k[2],k[1],k[3]),t]);

        @constraint(mp,[k in fData.brList, t in 2:T],p[k,t] <= (fData.g[k]*vhat[k[1],t]/(fData.τ1[k]^2) - fData.g[k]*wc[k,t] - fData.b[k]*ws[k,t]) + M[k]*(1 - mp[:x][k,t]));
        @constraint(mp,[k in fData.brList, t in 2:T],p[k,t] >= (fData.g[k]*vhat[k[1],t]/(fData.τ1[k]^2) - fData.g[k]*wc[k,t] - fData.b[k]*ws[k,t]) - M[k]*(1 - mp[:x][k,t]));
        @constraint(mp,[k in fData.brList, t in 2:T],q[k,t] <= ((-fData.b[k] - fData.bc[k]/2)*vhat[k[1],t]/(fData.τ1[k]^2) + fData.b[k]*wc[k,t] - fData.g[k]*ws[k,t]) + M[k]*(1 - mp[:x][k,t]));
        @constraint(mp,[k in fData.brList, t in 2:T],q[k,t] >= ((-fData.b[k] - fData.bc[k]/2)*vhat[k[1],t]/(fData.τ1[k]^2) + fData.b[k]*wc[k,t] - fData.g[k]*ws[k,t]) - M[k]*(1 - mp[:x][k,t]));
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

        @constraint(mp,[k in fData.brList, t in 2:T; fData.rateA[k] < Inf],[fData.rateA[k] * mp[:x][k,t]; socList1[k,t]] in SecondOrderCone());
        @constraint(mp,[k in fData.brList, t in 2:T], [(vhat[k[1],t]/(fData.τ1[k]^2) + vhat[k[2],t]/(fData.τ2[k]^2))/sqrt(2); socList2[k,t]] in SecondOrderCone());
        @constraint(mp,[k in fData.brList, t in 2:T],[tAux4[k,t]; socList3[k,t]] in SecondOrderCone());
        @constraint(mp,[i in fData.IDList, t in 2:T],[tAux6[i,t]; socList4[i,t]] in SecondOrderCone());
        @constraint(mp,[k in fData.brList, t in 2:T],[(vhat[k[1],t]/(fData.τ1[k]^2) + vhat[k[2],t]/(fData.τ2[k]^2))/sqrt(2); socList5[k,t]] in SecondOrderCone());

        @constraint(mp,[k in fData.brList, t in 2:T],sqrt((1-cos(θu[t][k]))/(θu[t][k])^2)*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t]) == tAux2[k,t]);
        @constraint(mp,[k in fData.brList, t in 2:T],cs[k,t] - 3/4 == tAux3[k,t]);
        @constraint(mp,[k in fData.brList, t in 2:T],5/4 - cs[k,t] == tAux4[k,t]);
        @constraint(mp,[k in fData.brList, t in 2:T],
                    ss[k,t] <= cos(θu[t][k]/2)*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t] - θu[t][k]/2) + sin(θu[t][k]/2));
        @constraint(mp,[k in fData.brList, t in 2:T],
                    ss[k,t] >= cos(θu[t][k]/2)*((θ[k[1],t] - fData.σ[k]) - θ[k[2],t] + θu[t][k]/2) - sin(θu[t][k]/2));
        @constraint(mp,[k in fData.brList, t in 2:T],(θ[k[1],t] - fData.σ[k]) - θ[k[2],t] <= θDmaxT[t][k] + 2*pi*(1 - mp[:x][k,t]));
        @constraint(mp,[k in fData.brList, t in 2:T],(θ[k[1],t] - fData.σ[k]) - θ[k[2],t] >= θDminT[t][k] - 2*pi*(1 - mp[:x][k,t]));
        @constraint(mp,[i in fData.IDList, t in 2:T],tAux5[i,t] == vhat[i,t] - 1/4);
        @constraint(mp,[i in fData.IDList, t in 2:T],tAux6[i,t] == vhat[i,t] + 1/4);
        @constraint(mp,[i in fData.IDList, t in 2:T], vhat[i,t] - (vmaxT[t][i] + vminT[t][i])*v[i,t] <= -vmaxT[t][i]*vminT[t][i]);
        @constraint(mp,[k in fData.brList, t in 2:T],
                    vv[k,t] >= vminT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vminT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));
        @constraint(mp,[k in fData.brList, t in 2:T],
                    vv[k,t] >= vmaxT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vmaxT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));
        @constraint(mp,[k in fData.brList, t in 2:T],
                    vv[k,t] <= vminT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vmaxT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vminT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));
        @constraint(mp,[k in fData.brList, t in 2:T],
                    vv[k,t] <= vmaxT[t][k[1]]*v[k[2],t]/(fData.τ1[k]*fData.τ2[k]) + vminT[t][k[2]]*v[k[1],t]/(fData.τ1[k]*fData.τ2[k]) - vmaxT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]));

        @constraint(mp,[k in fData.brList, t in 2:T],
                    wc[k,t] >= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + csmin[k,t]*vv[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k,t]);
        @constraint(mp,[k in fData.brList, t in 2:T],
                    wc[k,t] >= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + csmax[k,t]*vv[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k,t]);
        @constraint(mp,[k in fData.brList, t in 2:T],
                    wc[k,t] <= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + vv[k,t]*csmax[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmax[k,t]);
        @constraint(mp,[k in fData.brList, t in 2:T],
                    wc[k,t] <= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cs[k,t] + vv[k,t]*csmin[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*csmin[k,t]);

        @constraint(mp,[k in fData.brList, t in 2:T],
                    ws[k,t] >= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + ssmin[k,t]*vv[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k,t]);
        @constraint(mp,[k in fData.brList, t in 2:T],
                    ws[k,t] >= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + ssmax[k,t]*vv[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k,t]);
        @constraint(mp,[k in fData.brList, t in 2:T],
                    ws[k,t] <= vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + vv[k,t]*ssmax[k,t] - vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmax[k,t]);
        @constraint(mp,[k in fData.brList, t in 2:T],
                    ws[k,t] <= vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ss[k,t] + vv[k,t]*ssmin[k,t] - vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*ssmin[k,t]);

        @constraint(mp,[k in fData.brList, t in 2:T],ws[k,t] - tan(θDmaxT[t][k])*wc[k,t] <= 0);
        @constraint(mp,[k in fData.brList, t in 2:T],ws[k,t] - tan(θDminT[t][k])*wc[k,t] >= 0);

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

        @constraint(mp,[k in fData.brList, t in 2:T],
                    vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*(wc[k,t]*cos(θϕ[k,t]) + ws[k,t]*sin(θϕ[k,t])) -
                    vmaxT[t][k[2]]/fData.τ2[k]*cos(θδ[k,t])*vδ[k[2],t]/fData.τ2[k]*vhat[k[1],t]/(fData.τ1[k]^2) -
                    vmaxT[t][k[1]]/fData.τ1[k]*cos(θδ[k,t])*vδ[k[1],t]/fData.τ1[k]*vhat[k[2],t]/(fData.τ2[k]^2) >=
                    vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]) -
                    vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])));
        @constraint(mp,[k in fData.brList, t in 2:T],
                    vδ[k[1],t]*vδ[k[2],t]/(fData.τ1[k]*fData.τ2[k])*(wc[k,t]*cos(θϕ[k,t]) + ws[k,t]*sin(θϕ[k,t])) -
                    vminT[t][k[2]]/fData.τ2[k]*cos(θδ[k,t])*vδ[k[2],t]/fData.τ2[k]*vhat[k[1],t]/(fData.τ1[k]^2) -
                    vminT[t][k[1]]/fData.τ1[k]*cos(θδ[k,t])*vδ[k[1],t]/fData.τ1[k]*vhat[k[2],t]/(fData.τ2[k]^2) >=
                    -vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])*cos(θδ[k,t])*(vminT[t][k[1]]*vminT[t][k[2]]/(fData.τ1[k]*fData.τ2[k]) -
                    vmaxT[t][k[1]]*vmaxT[t][k[2]]/(fData.τ1[k]*fData.τ2[k])));

        # battery dynamics
        ep_var = @variable(mp, [i in fData.IDList, t in 2:T]);
        eq_var = @variable(mp, [i in fData.IDList, t in 2:T]);
        f_var = @variable(mp, [i in fData.IDList, t in 2:T]);
        I_var = @variable(mp, [i in fData.IDList, t in 1:T], lower_bound = 0);
        @constraint(mp, [i in fData.IDList, t in 1:T], I_var[i,t] <= mp[:y][i] * bData.cap[i]);
        @constraint(mp, [i in fData.IDList, t in 2:T], I_var[i,t] == I_var[i,t-1] - fData.Δt * f_var[i,t]);
        @constraint(mp, [i in fData.IDList, t in 2:T], [bData.uCap[i]*mp[:y][i]; [ep_var[i,t], eq_var[i,t]]] in SecondOrderCone());
        @constraint(mp, [i in fData.IDList, l in 1:length(bData.ηα[i]), t in 2:T], ep_var[i,t] <= bData.ηα[i][l]*f_var[i,t] + bData.ηβ[i][l]);

        # power flow balance
        dp_var = @variable(mp, [i in fData.IDList, t in 2:T]);
        dq_var = @variable(mp, [i in fData.IDList, t in 2:T]);
        h_var = @variable(mp, [i in fData.IDList, t in 2:T], lower_bound = 0);
        @constraint(mp,[i in fData.IDList, t in 2:T], sum(p[k,t] for k in branchDict1[i]) + vhat[i,t]*fData.gs[i]
                    + lpplus[i,t] - lpminus[i,t] == sphatsum[i,t] + h_var[i,t] - dp_var[i,t] + ep_var[i,t]);
        @constraint(mp,[i in fData.IDList, t in 2:T], sum(q[k,t] for k in branchDict1[i]) - vhat[i,t]*fData.bs[i]
                    + lqplus[i,t] - lqminus[i,t] == sqhatsum[i,t] - dq_var[i,t] + eq_var[i,t]);

        # ramping constraints
        @constraint(mp,[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t-1] <= fData.RU[i]);
        @constraint(mp,[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t-1] >= fData.RD[i]);

        # d_var constraints: demand uncertainty
        @constraint(mp, [i in fData.IDList, t in 2:T], dp_var[i,t] == (uData[i].DPmax[t] - uData[i].DP0[t])*uList[ω]["u_dp"][group_rev[i],t] + 
            (uData[i].DPmin[t] - uData[i].DP0[t])*uList[ω]["u_dm"][group_rev[i],t] + uData[i].DP0[t]);
        @constraint(mp, [i in fData.IDList, t in 2:T], dq_var[i,t] == (uData[i].DQmax[t] - uData[i].DQ0[t])*uList[ω]["u_dp"][group_rev[i],t] + 
            (uData[i].DQmin[t] - uData[i].DQ0[t])*uList[ω]["u_dm"][group_rev[i],t] + uData[i].DQ0[t]);

        # h_var constraints: renewable uncertainty
        @constraint(mp, [i in fData.IDList, t in 2:T; i in hData.hList], h_var[i,t] <= uData[i].RESP0[t]*((1 + expansion_factor*mp[:z][i]) + 
                uData[i].RESPmax*(uList[ω]["u_hp"][i,t]+expansion_factor*mp[:z][i]*uList[ω]["u_hp"][i,t]) - 
                uData[i].RESPmin*(uList[ω]["u_hm"][i,t]+expansion_factor*mp[:z][i]*uList[ω]["u_hm"][i,t])));
        @constraint(mp, [i in fData.IDList, t in 2:T; !(i in hData.hList)], h_var[i,t] == 0.0);    

        # objective value
        @constraint(mp, mp[:V] >= fData.cz*sum(lpplus[i,t] + lpminus[i,t] + lqplus[i,t] + lqminus[i,t] for i in fData.IDList for t in 2:T) + 
            sum(fData.cp[i].params[2]*sp[i,t] for i in fData.genIDList for t in 2:T));
    end
    return mp;
end

function solve_first(mp,fData,hData)
    # solve first-stage problem to obtain the solution
    optimize!(mp);
    objhat = objective_value(mp);

    sphat = Dict();
    sqhat = Dict();
    xhat = Dict();
    yhat = Dict();
    zhat = Dict();
    for i in fData.genIDList
        sphat[i] = value(mp[:sp][i]);
        sqhat[i] = value(mp[:sq][i]);
    end
    for t in 1:T
        for k in fData.brList
            xhat[k,t] = value(mp[:x][k,t]);
        end
    end
    for i in fData.IDList
        yhat[i] = value(mp[:y][i]);
    end
    for i in hData.hList
        zhat[i] = value(mp[:z][i]);
    end

    return objhat,sphat,sqhat,xhat,yhat,zhat;
end

function solve_second(subp, hData, groupDict)
    # solve second-stage problem to obtain the solution
    optimize!(subp);
    subp_obj = objective_value(subp);

    uDict = Dict();
    uDict["u_dp"] = Dict();
    uDict["u_dm"] = Dict();
    uDict["u_hp"] = Dict();
    uDict["u_hm"] = Dict();
    for t in 2:T
        for m in eachindex(groupDict[2])
            uDict["u_dp"][m,t] = value(subp[:u_dp][m,t]);
            uDict["u_dm"][m,t] = value(subp[:u_dm][m,t]);
        end
        for i in hData.hList
            uDict["u_hp"][i,t] = value(subp[:u_hp][i,t]);
            uDict["u_hm"][i,t] = value(subp[:u_hm][i,t]);
        end
    end
    return uDict, subp_obj;
end

function first_stage_cut(mp,fData,uData,hData,T,groupDict,vmaxT,vminT,θDmaxT,θDminT,expansion_factor,uList,iter_limit = 500, parallel_option = true)
    # for each scenario, solve the subproblem and generate cuts for the master problem
    cut_iter_bool = true;
    counter = 0;
    while (cut_iter_bool)&(counter <= iter_limit)
        # solve the first stage problem and obtain the solution
        optimize!(mp);
        objhat = objective_value(mp);

        Vhat = value(mp[:V]);
        sphat = Dict();
        sqhat = Dict();
        xhat = Dict();
        yhat = Dict();
        zhat = Dict();
        for i in fData.genIDList
            sphat[i] = value(mp[:sp][i]);
            sqhat[i] = value(mp[:sq][i]);
        end
        for t in 1:T
            for k in fData.brList
                xhat[k,t] = value(mp[:x][k,t]);
            end
        end
        for i in fData.IDList
            yhat[i] = value(mp[:y][i]);
        end
        for i in hData.hList
            zhat[i] = value(mp[:z][i]);
        end

        # based on this solution, solve the subproblem
        cut_iter_bool = false;
        if parallel_option
            obj_list, x_coeff_list, y_coeff_list, z_coeff_list = pmap(ω -> solve_dual_sub_uhat(fData,uData,hData,T,groupDict,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uList[ω]), 1:length(uList));
            for ω in eachindex(uList)
                x_coeff = x_coeff_list[ω];
                y_coeff = y_coeff_list[ω];
                z_coeff = z_coeff_list[ω];
                objV = obj_list[ω];

                if Vhat < objV - 1e-4
                    cut_iter_bool = true;
                    @constraint(mp, mp[:V] >= objV + sum(sum(x_coeff[k,t] * (mp[:x][k,t] - xhat[k,t]) for k in fData.brList) for t in 2:T) + 
                                        sum(y_coeff[i] * (mp[:y][i] - yhat[i]) for i in fData.IDList) + sum(z_coeff[i] * (mp[:z][i] - zhat[i]) for i in hData.hList));
                end
            end
        else
            for ω in eachindex(uList)
                subd = dual_sub_uhat(fData,uData,hData,T,groupDict,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uList[ω]);
                optimize!(subd);

                x_coeff = shadow_price.(subd[:xIni]);
                y_coeff = shadow_price.(subd[:yIni]);
                z_coeff = shadow_price.(subd[:zIni]);

                if Vhat < objective_value(subd) - 1e-4
                    cut_iter_bool = true;
                    @constraint(mp, mp[:V] >= objective_value(subd) + sum(sum(x_coeff[k,t] * (mp[:x][k,t] - xhat[k,t]) for k in fData.brList) for t in 2:T) + 
                                        sum(y_coeff[i] * (mp[:y][i] - yhat[i]) for i in fData.IDList) + sum(z_coeff[i] * (mp[:z][i] - zhat[i]) for i in hData.hList));
                end
            end
        end
        counter += 1;
    end

    # return the master problem with generated cuts
    return mp;
end