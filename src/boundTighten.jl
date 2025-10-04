function run_bt(fData,uData,hData,T,vmaxT,vminT,θDmaxT,θDminT,xLimit)
    # create the first-stage model
    # first-stage model without any scenarios
    θu = Dict();
    for t in 1:T
        θu[t] = Dict();
        for k in fData.brList
            θu[t][k] = max(abs(θDmaxT[t][k]),abs(θDminT[t][k]));
        end
    end
    # scaling issue exists for this primal problem.
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "NumericFocus" => 3, "BarConvTol" => 1e-6, "MIPGap" => 1e-6, "BarQCPConvTol" => 1e-6,
        "OptimalityTol" => 1e-6, "IntFeasTol" => 1e-6, "FeasibilityTol" => 1e-6, "OutputFlag" => 1, "Threads" => 1));
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

    # OPF variables
    @variable(mp,p[k in fData.brList]);
    @variable(mp,q[k in fData.brList]);
    @variable(mp,vminT[1][i] <= v[i in fData.IDList] <= vmaxT[1][i]);
    @variable(mp,vhat[i in fData.IDList]);
    @variable(mp,θ[i in fData.IDList]);

    # switching, battery setup, renewable expansion binary variables
    @variable(mp,x[k in fData.brList],Bin);
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

    @constraint(mp,lineConstrP1[k in fData.brList],p[k] <= (fData.g[k]*vhat[k[1]]/(fData.τ1[k]^2) - fData.g[k]*wc[k] - fData.b[k]*ws[k]) + 10 * fData.rateA[k] * (1 - x[k]));
    @constraint(mp,lineConstrP2[k in fData.brList],p[k] >= (fData.g[k]*vhat[k[1]]/(fData.τ1[k]^2) - fData.g[k]*wc[k] - fData.b[k]*ws[k]) - 10 * fData.rateA[k] * (1 - x[k]));
    @constraint(mp,lineConstrQ1[k in fData.brList],q[k] <= ((-fData.b[k] - fData.bc[k]/2)*vhat[k[1]]/(fData.τ1[k]^2) + fData.b[k]*wc[k] - fData.g[k]*ws[k]) + 10 * fData.rateA[k] * (1 - x[k]));
    @constraint(mp,lineConstrQ2[k in fData.brList],q[k] >= ((-fData.b[k] - fData.bc[k]/2)*vhat[k[1]]/(fData.τ1[k]^2) + fData.b[k]*wc[k] - fData.g[k]*ws[k]) - 10 * fData.rateA[k] * (1 - x[k]));
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

    @constraint(mp,socConstraint1[k in fData.brList; fData.rateA[k] < Inf],[fData.rateA[k] * x[k]; socList1[k]] in SecondOrderCone());
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
    @variable(mp, minimum(uData[i].DPmin) <= d_var_P[i in fData.IDList] <= maximum(uData[i].DPmax));
    @variable(mp, minimum(uData[i].DQmin) <= d_var_Q[i in fData.IDList] <= maximum(uData[i].DQmax));
    @variable(mp, h_var[i in fData.IDList] >= 0.0);
    @constraint(mp, h_bound_0[i in fData.IDList;!(i in hData.hList)], h_var[i] <= 0.0);
    @constraint(mp, h_bound_1[i in hData.hList], h_var[i] <= maximum(uData[i].RESP0) * (1 + expansion_factor * z[i]) * (1 + uData[i].RESPmax));
    @variable(mp,e_var_P[i in fData.IDList] >= 0);
    @variable(mp,e_var_Q[i in fData.IDList] >= 0);
    @constraint(mp, bat_out[i in fData.IDList], [bData.uCap[i]*y[i]; [e_var_P[i], e_var_Q[i]]] in SecondOrderCone());

    @constraint(mp,totalP[i in fData.IDList], sum(p[k] for k in branchDict1[i]) + vhat[i]*fData.gs[i] == sphatsum[i] - d_var_P[i] + h_var[i] + e_var_P[i]);
    @constraint(mp,totalQ[i in fData.IDList], sum(q[k] for k in branchDict1[i]) - vhat[i]*fData.bs[i] == sqhatsum[i] - d_var_Q[i] + e_var_Q[i]);

    # add the battery/switching/renewable expansion
    @constraint(mp,switching_rev[k in fData.brList;k[1] < k[2]], x[k] == x[(k[2],k[1],k[3])]);
    @constraint(mp,switching_lim, sum(x[k] for k in fData.brList if k[1] < k[2]) >= length(fData.brList)/2 - xLimit);

    # set objective function and iteratively obtain the bounds
    vmax_new = Dict();
    vmin_new = Dict();
    θDmax_new = Dict();
    θDmin_new = Dict();
    for i in fData.IDList
        @objective(mp, Max, v[i]);
        optimize!(mp);
        vmax_new[i] = objective_value(mp);
        @objective(mp, Min, v[i]);
        optimize!(mp);
        vmin_new[i] = objective_value(mp);
    end
    for k in fData.brList
        @objective(mp, Max, θ[k[1]] - θ[k[2]]);
        optimize!(mp);
        θDmax_new[k] = objective_value(mp);
        @objective(mp, Min, θ[k[1]] - θ[k[2]]);
        optimize!(mp);
        θDmin_new[k] = objective_value(mp);
    end

    return vmax_new, vmin_new, θDmax_new, θDmin_new;
end