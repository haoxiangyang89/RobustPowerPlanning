# construct the subproblem for the cutting-plane algorithm for the separation
function dual_sub(fData,uData,hData,T,groupDict,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uhat)
    # create the QC subdroblem in the MISOCP
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
    @constraint(subd, pihp_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihp[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pihp_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihp[i,t,j] * signDict[j] <= abs(uData[i].RESP0[t] * uData[i].RESPmax));
    @constraint(subd, pihm_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihm[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pihm_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihm[i,t,j] * signDict[j] <= abs(uData[i].RESP0[t] * uData[i].RESPmin));
    @constraint(subd, pizetap2_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap2[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetap2_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap2[i,t,j] * signDict[j] <= 1);
    @constraint(subd, pizetam2_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam2[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetam2_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam2[i,t,j] * signDict[j] <= 1);
    @constraint(subd, pizetap3_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap3[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetap3_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap3[i,t,j] * signDict[j] <= 1);
    @constraint(subd, pizetam3_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam3[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetam3_bound[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam3[i,t,j] * signDict[j] <= 1);

    @constraint(subd, pipdp_sign[i in fData.IDList, t in 2:T, j in 1:4], pipdp[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pipdm_sign[i in fData.IDList, t in 2:T, j in 1:4], pipdm[i,t,j] * signDict[j] >= 0);
    @constraint(subd, piqdp_sign[i in fData.IDList, t in 2:T, j in 1:4], piqdp[i,t,j] * signDict[j] >= 0);
    @constraint(subd, piqdm_sign[i in fData.IDList, t in 2:T, j in 1:4], piqdm[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pipdp_bound[i in fData.IDList, t in 2:T, j in 1:4], pipdp[i,t,j] * signDict[j] <= abs(uData[i].DPmax[t] - uData[i].DP0[t]));
    @constraint(subd, pipdm_bound[i in fData.IDList, t in 2:T, j in 1:4], pipdm[i,t,j] * signDict[j] <= abs(-uData[i].DPmin[t] + uData[i].DP0[t]));
    @constraint(subd, piqdp_bound[i in fData.IDList, t in 2:T, j in 1:4], piqdp[i,t,j] * signDict[j] <= abs(uData[i].DQmax[t] - uData[i].DQ0[t]));
    @constraint(subd, piqdm_bound[i in fData.IDList, t in 2:T, j in 1:4], piqdm[i,t,j] * signDict[j] <= abs(-uData[i].DQmin[t] + uData[i].DQ0[t]));

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
    u_dp = uhat["u_dp"];
    u_dm = uhat["u_dm"];
    u_hp = uhat["u_hp"];
    u_hm = uhat["u_hm"];
    lambdaMu = fData.cz;
    lambdaMl = -fData.cz;
    obj_term = @expression(subd, fData.cz*sum(lpplus[i,t] + lpminus[i,t] + lqplus[i,t] + lqminus[i,t] for i in fData.IDList for t in 2:T) + 
                sum(fData.cp[i].params[2]*sp[i,t] for i in fData.genIDList for t in 2:T));
    @expression(subd, obj_rest_term, sum(sum(((pihp[i,t,1] + pizetap2[i,t,1]) * lambdaMl + pizetap3[i,t,2] * lambdaMu) * u_hp[i,t] + 
                ((pihp[i,t,4] + pizetap2[i,t,4]) * lambdaMl + pizetap3[i,t,3] * lambdaMu) * (u_hp[i,t] - 1) +
                ((pihm[i,t,1] + pizetam2[i,t,1]) * lambdaMl + pizetam3[i,t,2] * lambdaMu) * u_hm[i,t] + 
                ((pihm[i,t,4] + pizetam2[i,t,4]) * lambdaMl + pizetam3[i,t,3] * lambdaMu) * (u_hm[i,t] - 1) for i in fData.IDList if i in hData.hList) + 
            sum(((pipdp[i,t,1] + piqdp[i,t,1]) * lambdaMl + (pipdp[i,t,2] + piqdp[i,t,2]) * lambdaMu) * u_dp[group_rev[i],t] + 
                ((pipdp[i,t,4] + piqdp[i,t,4]) * lambdaMl + (pipdp[i,t,3] + piqdp[i,t,3]) * lambdaMu) * (u_dp[group_rev[i],t] - 1) + 
                ((pipdm[i,t,1] + piqdm[i,t,1]) * lambdaMl + (pipdm[i,t,2] + piqdm[i,t,2]) * lambdaMu) * u_dm[group_rev[i],t] + 
                ((pipdm[i,t,4] + piqdm[i,t,4]) * lambdaMl + (pipdm[i,t,3] + piqdm[i,t,3]) * lambdaMu) * (u_dm[group_rev[i],t] - 1) for i in fData.IDList) for t in 2:T));
    @objective(subd, Min, fData.cz*sum(lpplus[i,t] + lpminus[i,t] + lqplus[i,t] + lqminus[i,t] for i in fData.IDList for t in 2:T) + 
        sum(fData.cp[i].params[2]*sp[i,t] for i in fData.genIDList for t in 2:T) + 
        sum(sum(((pihp[i,t,1] + pizetap2[i,t,1]) * lambdaMl + pizetap3[i,t,2] * lambdaMu) * u_hp[i,t] + 
                ((pihp[i,t,4] + pizetap2[i,t,4]) * lambdaMl + pizetap3[i,t,3] * lambdaMu) * (u_hp[i,t] - 1) +
                ((pihm[i,t,1] + pizetam2[i,t,1]) * lambdaMl + pizetam3[i,t,2] * lambdaMu) * u_hm[i,t] + 
                ((pihm[i,t,4] + pizetam2[i,t,4]) * lambdaMl + pizetam3[i,t,3] * lambdaMu) * (u_hm[i,t] - 1) for i in fData.IDList if i in hData.hList) + 
            sum(((pipdp[i,t,1] + piqdp[i,t,1]) * lambdaMl + (pipdp[i,t,2] + piqdp[i,t,2]) * lambdaMu) * u_dp[group_rev[i],t] + 
                ((pipdp[i,t,4] + piqdp[i,t,4]) * lambdaMl + (pipdp[i,t,3] + piqdp[i,t,3]) * lambdaMu) * (u_dp[group_rev[i],t] - 1) + 
                ((pipdm[i,t,1] + piqdm[i,t,1]) * lambdaMl + (pipdm[i,t,2] + piqdm[i,t,2]) * lambdaMu) * u_dm[group_rev[i],t] + 
                ((pipdm[i,t,4] + piqdm[i,t,4]) * lambdaMl + (pipdm[i,t,3] + piqdm[i,t,3]) * lambdaMu) * (u_dm[group_rev[i],t] - 1) for i in fData.IDList) for t in 2:T));
    # @objective(subd, Min, fData.cz*sum(lpplus[i,t] + lpminus[i,t] + lqplus[i,t] + lqminus[i,t] for i in fData.IDList for t in 2:T) + 
    #     sum(fData.cp[i].params[2]*sp[i,t] for i in fData.genIDList for t in 2:T));
    @expression(subd, u_hp_coeff[i in fData.IDList, t in 2:T; i in hData.hList], ((pihp[i,t,1] + pizetap2[i,t,1]) * lambdaMl + pizetap3[i,t,2] * lambdaMu) +
                        ((pihp[i,t,4] + pizetap2[i,t,4]) * lambdaMl + pizetap3[i,t,3] * lambdaMu));
    @expression(subd, u_hm_coeff[i in fData.IDList, t in 2:T; i in hData.hList], ((pihm[i,t,1] + pizetam2[i,t,1]) * lambdaMl + pizetam3[i,t,2] * lambdaMu) +
                        ((pihm[i,t,4] + pizetam2[i,t,4]) * lambdaMl + pizetam3[i,t,3] * lambdaMu));
    @expression(subd, u_dp_coeff[ig in groupList, t in 2:T], sum(((pipdp[i,t,1] + piqdp[i,t,1]) * lambdaMl + (pipdp[i,t,2] + piqdp[i,t,2]) * lambdaMu) +
                        ((pipdp[i,t,4] + piqdp[i,t,4]) * lambdaMl + (pipdp[i,t,3] + piqdp[i,t,3]) * lambdaMu) for i in groupDict[2][ig]));
    @expression(subd, u_dm_coeff[ig in groupList, t in 2:T], sum(((pipdm[i,t,1] + piqdm[i,t,1]) * lambdaMl + (pipdm[i,t,2] + piqdm[i,t,2]) * lambdaMu) +
                        ((pipdm[i,t,4] + piqdm[i,t,4]) * lambdaMl + (pipdm[i,t,3] + piqdm[i,t,3]) * lambdaMu) for i in groupDict[2][ig]));

    return subd, obj_term;
end

function dual_sub_modified(fData,uData,hData,T,groupDict,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uhat, regular_const = 0.01)
    # create the QC subdroblem in the MISOCP
    # dual's dual should be similar to the primal
    θu = Dict();
    for t in 1:T
        θu[t] = Dict();
        for k in fData.brList
            θu[t][k] = max(abs(θDmaxT[t][k]),abs(θDminT[t][k]));
        end
    end
    # scaling issue exists for this primal problem.
    #subd = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag" => 1, "Threads" => 1));
    subd = Model(COPT.Optimizer);

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
    @constraint(subd, pihp_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihp[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pihm_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pihm[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetap2_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap2[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetam2_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam2[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetap3_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetap3[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pizetam3_sign[i in fData.IDList, t in 2:T, j in 1:4; i in hData.hList], pizetam3[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pipdp_sign[i in fData.IDList, t in 2:T, j in 1:4], pipdp[i,t,j] * signDict[j] >= 0);
    @constraint(subd, pipdm_sign[i in fData.IDList, t in 2:T, j in 1:4], pipdm[i,t,j] * signDict[j] >= 0);
    @constraint(subd, piqdp_sign[i in fData.IDList, t in 2:T, j in 1:4], piqdp[i,t,j] * signDict[j] >= 0);
    @constraint(subd, piqdm_sign[i in fData.IDList, t in 2:T, j in 1:4], piqdm[i,t,j] * signDict[j] >= 0);

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
    u_dp = uhat["u_dp"];
    u_dm = uhat["u_dm"];
    u_hp = uhat["u_hp"];
    u_hm = uhat["u_hm"];
    lambdaMu = fData.cz;
    lambdaMl = -fData.cz;
    obj_term = @expression(subd, obj_term, fData.cz*sum(lpplus[i,t] + lpminus[i,t] + lqplus[i,t] + lqminus[i,t] for i in fData.IDList for t in 2:T) + 
                sum(fData.cp[i].params[2]*sp[i,t] for i in fData.genIDList for t in 2:T) + sum(sum(((pihp[i,t,1] + pizetap2[i,t,1]) * lambdaMl + pizetap3[i,t,2] * lambdaMu) * u_hp[i,t] + 
                ((pihp[i,t,4] + pizetap2[i,t,4]) * lambdaMl + pizetap3[i,t,3] * lambdaMu) * (u_hp[i,t] - 1) +
                ((pihm[i,t,1] + pizetam2[i,t,1]) * lambdaMl + pizetam3[i,t,2] * lambdaMu) * u_hm[i,t] + 
                ((pihm[i,t,4] + pizetam2[i,t,4]) * lambdaMl + pizetam3[i,t,3] * lambdaMu) * (u_hm[i,t] - 1) for i in fData.IDList if i in hData.hList) + 
            sum(((pipdp[i,t,1] + piqdp[i,t,1]) * lambdaMl + (pipdp[i,t,2] + piqdp[i,t,2]) * lambdaMu) * u_dp[group_rev[i],t] + 
                ((pipdp[i,t,4] + piqdp[i,t,4]) * lambdaMl + (pipdp[i,t,3] + piqdp[i,t,3]) * lambdaMu) * (u_dp[group_rev[i],t] - 1) + 
                ((pipdm[i,t,1] + piqdm[i,t,1]) * lambdaMl + (pipdm[i,t,2] + piqdm[i,t,2]) * lambdaMu) * u_dm[group_rev[i],t] + 
                ((pipdm[i,t,4] + piqdm[i,t,4]) * lambdaMl + (pipdm[i,t,3] + piqdm[i,t,3]) * lambdaMu) * (u_dm[group_rev[i],t] - 1) for i in fData.IDList) for t in 2:T));
    regularize_term = @expression(subd, regularize_term, sum(sum(sum((pihp[i,t,j] + pihm[i,t,j] + pizetap2[i,t,j] + pizetam2[i,t,j] + pizetap3[i,t,j] + pizetam3[i,t,j]) *
             signDict[j] for j in 1:4) for i in hData.hList) + 
            sum(sum((pipdp[i,t,j] + pipdm[i,t,j] + piqdp[i,t,j] + piqdm[i,t,j]) * 
             signDict[j] for j in 1:4) for i in fData.IDList) for t in 2:T));
    @objective(subd, Min, regular_const * regularize_term + fData.cz*sum(lpplus[i,t] + lpminus[i,t] + lqplus[i,t] + lqminus[i,t] for i in fData.IDList for t in 2:T) + 
        sum(fData.cp[i].params[2]*sp[i,t] for i in fData.genIDList for t in 2:T) + 
        sum(sum(((pihp[i,t,1] + pizetap2[i,t,1]) * lambdaMl + pizetap3[i,t,2] * lambdaMu) * u_hp[i,t] + 
                ((pihp[i,t,4] + pizetap2[i,t,4]) * lambdaMl + pizetap3[i,t,3] * lambdaMu) * (u_hp[i,t] - 1) +
                ((pihm[i,t,1] + pizetam2[i,t,1]) * lambdaMl + pizetam3[i,t,2] * lambdaMu) * u_hm[i,t] + 
                ((pihm[i,t,4] + pizetam2[i,t,4]) * lambdaMl + pizetam3[i,t,3] * lambdaMu) * (u_hm[i,t] - 1) for i in fData.IDList if i in hData.hList) + 
            sum(((pipdp[i,t,1] + piqdp[i,t,1]) * lambdaMl + (pipdp[i,t,2] + piqdp[i,t,2]) * lambdaMu) * u_dp[group_rev[i],t] + 
                ((pipdp[i,t,4] + piqdp[i,t,4]) * lambdaMl + (pipdp[i,t,3] + piqdp[i,t,3]) * lambdaMu) * (u_dp[group_rev[i],t] - 1) + 
                ((pipdm[i,t,1] + piqdm[i,t,1]) * lambdaMl + (pipdm[i,t,2] + piqdm[i,t,2]) * lambdaMu) * u_dm[group_rev[i],t] + 
                ((pipdm[i,t,4] + piqdm[i,t,4]) * lambdaMl + (pipdm[i,t,3] + piqdm[i,t,3]) * lambdaMu) * (u_dm[group_rev[i],t] - 1) for i in fData.IDList) for t in 2:T));
    # @objective(subd, Min, fData.cz*sum(lpplus[i,t] + lpminus[i,t] + lqplus[i,t] + lqminus[i,t] for i in fData.IDList for t in 2:T) + 
    #     sum(fData.cp[i].params[2]*sp[i,t] for i in fData.genIDList for t in 2:T));
    @expression(subd, u_hp_coeff[i in fData.IDList, t in 2:T; i in hData.hList], ((pihp[i,t,1] + pizetap2[i,t,1]) * lambdaMl + pizetap3[i,t,2] * lambdaMu) +
                        ((pihp[i,t,4] + pizetap2[i,t,4]) * lambdaMl + pizetap3[i,t,3] * lambdaMu));
    @expression(subd, u_hm_coeff[i in fData.IDList, t in 2:T; i in hData.hList], ((pihm[i,t,1] + pizetam2[i,t,1]) * lambdaMl + pizetam3[i,t,2] * lambdaMu) +
                        ((pihm[i,t,4] + pizetam2[i,t,4]) * lambdaMl + pizetam3[i,t,3] * lambdaMu));
    @expression(subd, u_dp_coeff[ig in groupList, t in 2:T], sum(((pipdp[i,t,1] + piqdp[i,t,1]) * lambdaMl + (pipdp[i,t,2] + piqdp[i,t,2]) * lambdaMu) +
                        ((pipdp[i,t,4] + piqdp[i,t,4]) * lambdaMl + (pipdp[i,t,3] + piqdp[i,t,3]) * lambdaMu) for i in groupDict[2][ig]));
    @expression(subd, u_dm_coeff[ig in groupList, t in 2:T], sum(((pipdm[i,t,1] + piqdm[i,t,1]) * lambdaMl + (pipdm[i,t,2] + piqdm[i,t,2]) * lambdaMu) +
                        ((pipdm[i,t,4] + piqdm[i,t,4]) * lambdaMl + (pipdm[i,t,3] + piqdm[i,t,3]) * lambdaMu) for i in groupDict[2][ig]));

    return subd, obj_term;
end

function dual_master(fData,hData,T,groupDict)
    Γ = Dict("d" => Γd, "h" => Γh);
    groupList = 1:length(groupDict[2]);
    group_rev = Dict();
    for i in fData.IDList
        for m in groupList
            if i in groupDict[2][m]
                group_rev[i] = m;
            end
        end
    end
    # obtain the master problem for the dual
    dp_m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "NumericFocus" => 3, "BarConvTol" => 1e-6, "MIPGap" => 1e-6, "BarQCPConvTol" => 1e-6,
        "OptimalityTol" => 1e-6, "IntFeasTol" => 1e-6, "FeasibilityTol" => 1e-6, "OutputFlag" => 0, "Threads" => 1));
    
    @variable(dp_m, u_hp[i in hData.hList, t in 2:T], Bin);
    @variable(dp_m, u_hm[i in hData.hList, t in 2:T], Bin);
    @variable(dp_m, u_dp[m in groupList, t in 2:T], Bin);
    @variable(dp_m, u_dm[m in groupList, t in 2:T], Bin);

    @constraint(dp_m, uncertain_budget_d, sum(u_dp[m,t] + u_dm[m,t] for m in groupList for t in 2:T) <= Γ["d"]);
    @constraint(dp_m, one_extreme_pt_d[m in groupList, t in 2:T], u_dp[m,t] + u_dm[m,t] <= 1);
    @constraint(dp_m, uncertain_budget_h, sum(u_hp[i,t] + u_hm[i,t] for i in hData.hList for t in 2:T) <= Γ["h"]);
    @constraint(dp_m, one_extreme_pt_h[i in hData.hList, t in 2:T], u_hp[i,t] + u_hm[i,t] <= 1);

    # provide an upper bound on the objective value
    @variable(dp_m, V <= fData.cz * sum(abs(uData[i].DP0[t]) + abs(uData[i].DQ0[t]) for t in 2:T for i in fData.IDList));

    @objective(dp_m, Max, V);

    return dp_m;
end

# construct the master problem for the cutting-plane algorithm for the separation
function dual_master_cb(fData,uData,hData,T,groupDict,Γd,Γh,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat)
    Γ = Dict("d" => Γd, "h" => Γh);
    groupList = 1:length(groupDict[2]);
    group_rev = Dict();
    for i in fData.IDList
        for m in groupList
            if i in groupDict[2][m]
                group_rev[i] = m;
            end
        end
    end
    # obtain the master problem for the dual
    dp_m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV), "NumericFocus" => 3, "BarConvTol" => 1e-6, "MIPGap" => 1e-6, "BarQCPConvTol" => 1e-6,
        "OptimalityTol" => 1e-6, "IntFeasTol" => 1e-6, "FeasibilityTol" => 1e-6, "OutputFlag" => 1, "Threads" => 1));

    # set up binary variables to represent the extreme points of the uncertainty set
    @variable(dp_m, u_hp[i in hData.hList, t in 2:T], Bin);
    @variable(dp_m, u_hm[i in hData.hList, t in 2:T], Bin);
    @variable(dp_m, u_dp[m in groupList, t in 2:T], Bin);
    @variable(dp_m, u_dm[m in groupList, t in 2:T], Bin);

    # provide an upper bound on the objective value
    @variable(dp_m, V <= fData.cz * sum(abs(uData[i].DP0[t]) + abs(uData[i].DQ0[t]) for t in 2:T for i in fData.IDList));

    # set up the constraints for the uncertainty set
    @constraint(dp_m, uncertain_budget_d, sum(u_dp[m,t] + u_dm[m,t] for m in groupList for t in 2:T) <= Γ["d"]);
    @constraint(dp_m, one_extreme_pt_d[m in groupList, t in 2:T], u_dp[m,t] + u_dm[m,t] <= 1);
    @constraint(dp_m, uncertain_budget_h, sum(u_hp[i,t] + u_hm[i,t] for i in hData.hList for t in 2:T) <= Γ["h"]);
    @constraint(dp_m, one_extreme_pt_h[i in hData.hList, t in 2:T], u_hp[i,t] + u_hm[i,t] <= 1);

    @objective(dp_m, Max, V);

    # set up the callback function to add the lazy constraints
    function QC_callback(cb_data)
        status = callback_node_status(cb_data, dp_m);
        uhat = Dict("u_hp" => Dict(), "u_hm" => Dict(), "u_dp" => Dict(), "u_dm" => Dict());
        if status == MOI.CALLBACK_MODE_STATUS_INTEGER
            # obtain the integer solutions
            for m in groupList
                for t in 2:T
                    uhat["u_dp"][m,t] = callback_value(cb_data, dp_m[:u_dp][m,t]);
                    uhat["u_dm"][m,t] = callback_value(cb_data, dp_m[:u_dm][m,t]);
                end
            end
            for i in hData.hList
                for t in 2:T
                    uhat["u_hp"][i,t] = callback_value(cb_data, dp_m[:u_hp][i,t]);
                    uhat["u_hm"][i,t] = callback_value(cb_data, dp_m[:u_hm][i,t]);
                end
            end
            # solve the subproblem
            subd = dual_subproblem(fData,uData,hData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uhat);
            optimize!(subd);
            subd_status = termination_status(subd);
            if subd_status == MOI.OPTIMAL
                obj_subd = objective_value(subd);
                V_val = callback_value(cb_data, dp_m[:V]);
                if obj_subd > V_val + 1e-4
                    # obtain the dual variables
                    πhat = Dict("u_dp" => Dict(), "u_dm" => Dict(), "u_hp" => Dict(), "u_hm" => Dict());
                    lambdaMu = 5 * fData.cz;
                    lambdaMl = -5 * fData.cz;
                    for m in groupList
                        for t in 2:T
                            if uhat["u_dp"][m,t] == 0
                                πhat["u_dp"][m,t] = sum((subd[:pipdp][i,t,1] + subd[:piqdp][i,t,1]) * lambdaMl + 
                                                    (subd[:pipdp][i,t,2] + subd[:piqdp][i,t,2]) * lambdaMu for i in groupDict[m]);
                            else
                                πhat["u_dp"][m,t] = sum((subd[:pipdp][i,t,4] + subd[:piqdp][i,t,4]) * lambdaMl + 
                                                    (subd[:pipdp][i,t,3] + subd[:piqdp][i,t,3]) * lambdaMu for i in groupDict[m]);
                            end
                            if uhat["u_dm"][m,t] == 0
                                πhat["u_dm"][m,t] = sum((subd[:pipdm][i,t,1] + subd[:piqdm][i,t,1]) * lambdaMl + 
                                                    (subd[:pipdm][i,t,2] + subd[:piqdm][i,t,2]) * lambdaMu for i in groupDict[m]);
                            else
                                πhat["u_dm"][m,t] = sum((subd[:pipdm][i,t,4] + subd[:piqdm][i,t,4]) * lambdaMl + 
                                                    (subd[:pipdm][i,t,3] + subd[:piqdm][i,t,3]) * lambdaMu for i in groupDict[m]);
                            end
                        end
                    end
                    for i in hData.hList
                        if uhat["u_hp"][i,t] == 0
                            πhat["u_hp"][i,t] = sum((subd[:pihp][i,t,1] + subd[:pizetap2][i,t,1]) * lambdaMl + 
                                                subd[:pizetap3][i,t,2] * lambdaMu);
                        else
                            πhat["u_hp"][i,t] = sum((subd[:pihp][i,t,4] + subd[:pizetap2][i,t,4]) * lambdaMl + 
                                                subd[:pizetap3][i,t,3] * lambdaMu);
                        end
                        if uhat["u_hm"][i,t] == 0
                            πhat["u_hm"][i,t] = sum((subd[:pihm][i,t,1] + subd[:pizetam2][i,t,1]) * lambdaMl + 
                                                subd[:pizetam3][i,t,2] * lambdaMu);
                        else
                            πhat["u_hm"][i,t] = sum((subd[:pihm][i,t,4] + subd[:pizetam2][i,t,4]) * lambdaMl + 
                                                subd[:pizetam3][i,t,3] * lambdaMu);
                        end
                    end

                    # generate the lazy constraint
                    lazy_constr = @build_constraint(dp_m[:V] <= obj_subd + 
                        sum(πhat["u_dp"][m,t]*(dp_m[:u_dp][m,t] - uhat["u_dp"][m,t]) + πhat["u_dm"][m,t]*(dp_m[:u_dm][m,t] - uhat["u_dm"][m,t]) for m in groupList for t in 2:T) + 
                        sum(πhat["u_hp"][i,t]*(dp_m[:u_hp][i,t] - uhat["u_hp"][i,t]) + πhat["u_hm"][i,t]*(dp_m[:u_hm][i,t] - uhat["u_hm"][i,t]) for i in hData.hList for t in 2:T));
                    MOI.submit(dp_m, MOI.LazyConstraint(cb_data), lazy_constr);
                end
            end
        end
    end

    set_attribute(dp_m, MOI.LazyConstraintCallback(), QC_callback);

    return dp_m;
end