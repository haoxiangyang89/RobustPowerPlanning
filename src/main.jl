function solve_process(fData, uData, hData, bData, T, vmaxT, vminT, θDmaxT, θDminT, xLimit, groupDict, Γd, Γh, expansion_factor, no_threads = 1)
    Γ = Dict("d" => Γd, "h" => Γh);

    # create the first-stage model
    mp = createFirst(fData,uData,hData,bData,T,vmaxT,vminT,θDmaxT,θDminT,xLimit,no_threads);

    # append the all-zero nominal scenario
    uDict = Dict();
    uDict["u_dp"] = Dict();
    uDict["u_dm"] = Dict();
    uDict["u_hp"] = Dict();
    uDict["u_hm"] = Dict();
    for t in 2:T
        for m in eachindex(groupDict[2])
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
        mp = appendScen(mp,fData,uData,hData,T,groupDict,vmaxT,vminT,θDmaxT,θDminT,expansion_factor,[uDict]);

        # solve the first stage solution to obtain the initial solution
        objhat,sphat,sqhat,xhat,yhat,zhat = solve_first(mp,fData,hData);
        tFirst_end = time();
        tFirst_elapsed = tFirst_end - tFirst_start;

        # feed the first-stage solution to the second stage problem
        # subp = createSecond(fData,uData,hData,T,groupDict[ci],Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uDict);
        tSecond_start = time();
        subp_dual = createSecond_dual(fData,uData,hData,bData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,no_threads);

        # solve the second-stage problem and obtain the worst-case scenario
        uDict, subp_obj = solve_second(subp_dual, hData, groupDict);
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
    return objhat,sphat,sqhat,xhat,yhat,zhat,uList,solveInfo;
end

function solve_process_para(fData, uData, hData, bData, T, vmaxT, vminT, θDmaxT, θDminT, xLimit, groupDict, Γd, Γh, expansion_factor, k_switch, t_switch = 0, no_threads = 1)
    Γ = Dict("d" => Γd, "h" => Γh);

    # create the first-stage model
    mp = createFirst(fData,uData,hData,bData,T,vmaxT,vminT,θDmaxT,θDminT,xLimit);
    if t_switch != 0
        @constraint(mp, mp[:x][k_switch,t_switch] == 0);
    end

    # append the all-zero nominal scenario
    uDict = Dict();
    uDict["u_dp"] = Dict();
    uDict["u_dm"] = Dict();
    uDict["u_hp"] = Dict();
    uDict["u_hm"] = Dict();
    for t in 2:T
        for m in eachindex(groupDict[2])
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
    iter_no = 0;

    while keepIter
        iter_no += 1;
        # append the scenario to the first stage
        tFirst_start = time();
        push!(uList, uDict);
        mp = appendScen(mp,fData,uData,hData,T,groupDict,vmaxT,vminT,θDmaxT,θDminT,expansion_factor,[uDict]);

        # solve the first stage solution to obtain the initial solution
        objhat,sphat,sqhat,xhat,yhat,zhat = solve_first(mp,fData,hData);
        tFirst_end = time();
        tFirst_elapsed = tFirst_end - tFirst_start;

        # feed the first-stage solution to the second stage problem
        # subp = createSecond(fData,uData,hData,T,groupDict[ci],Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uDict);
        tSecond_start = time();
        subp_dual = createSecond_dual(fData,uData,hData,bData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat);

        # solve the second-stage problem and obtain the worst-case scenario
        uDict, subp_obj = solve_second(subp_dual,hData, groupDict);
        tSecond_end = time();
        tSecond_elapsed = tSecond_end - tSecond_start;

        # tell when to stop: if the scenario has been generated
        if uDict in uList
            keepIter = false;
        end

        push!(solveInfo["iter"], [objhat, tFirst_elapsed, subp_obj, tSecond_elapsed]);
        println("-------------------------- Iteration $(iter_no), relaxation lb = $(objhat), worst-case ub = $(subp_obj) --------------------------");
    end

    # obtain the optimal value and optimal solution
    objhat,sphat,sqhat,xhat,yhat,zhat = solve_first(mp,fData,hData);
    return objhat,sphat,sqhat,xhat,yhat,zhat,uList,solveInfo;
end

function solve_process_cut(fData, uData, hData, bData, T, vmaxT, vminT, θDmaxT, θDminT, xLimit, groupDict, Γd, Γh, expansion_factor, no_threads = 1)
    Γ = Dict("d" => Γd, "h" => Γh);

    # create the first-stage model
    mp = createFirst(fData,uData,hData,bData,T,vmaxT,vminT,θDmaxT,θDminT,xLimit,no_threads);

    # append the all-zero nominal scenario
    uDict = Dict();
    uDict["u_dp"] = Dict();
    uDict["u_dm"] = Dict();
    uDict["u_hp"] = Dict();
    uDict["u_hm"] = Dict();
    for t in 2:T
        for m in eachindex(groupDict[2])
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
        # generate cuts for each appended scenario
        tFirst_start = time();
        push!(uList, uDict);
        mp = first_stage_cut(mp,fData,uData,hData,bData,T,groupDict,vmaxT,vminT,θDmaxT,θDminT,expansion_factor,uList,500);

        # solve the first stage solution to obtain the initial solution
        objhat,sphat,sqhat,xhat,yhat,zhat = solve_first(mp,fData,hData);
        tFirst_end = time();
        tFirst_elapsed = tFirst_end - tFirst_start;

        # feed the first-stage solution to the second stage problem
        # subp = createSecond(fData,uData,hData,T,groupDict[ci],Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uDict);
        tSecond_start = time();
        subp_dual = createSecond_dual(fData,uData,hData,bData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,no_threads);

        # solve the second-stage problem and obtain the worst-case scenario
        uDict, subp_obj = solve_second(subp_dual, hData, groupDict);
        tSecond_end = time();
        tSecond_elapsed = tSecond_end - tSecond_start;

        # tell when to stop: if the scenario has been generated
        if uDict in uList
            keepIter = false;
        end

        push!(solveInfo["iter"], [objhat, tFirst_elapsed, subp_obj, tSecond_elapsed]);
    end    

    objhat,sphat,sqhat,xhat,yhat,zhat = solve_first(mp,fData,hData);
    return objhat,sphat,sqhat,xhat,yhat,zhat,uList,solveInfo;
end

function test_second(fData,uData,hData,bData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,gamma_results)
    obj_Dict = Dict();

    uDict = Dict();
    uDict["u_dp"] = Dict();
    uDict["u_dm"] = Dict();
    uDict["u_hp"] = Dict();
    uDict["u_hm"] = Dict();
    for t in 2:T
        for m in eachindex(groupDict[ci][2])
            uDict["u_dp"][m,t] = randn()*0.3;
            uDict["u_dm"][m,t] = randn()*0.3;
        end
        for i in hData.hList
            uDict["u_hp"][i,t] = randn()*0.3;
            uDict["u_hm"][i,t] = randn()*0.3;
        end
    end

    for gd in 1:3
        for gh in 1:3
            sphat = gamma_results[gd,gh][2];
            sqhat = gamma_results[gd,gh][3];
            xhat = gamma_results[gd,gh][4];
            yhat = gamma_results[gd,gh][5];
            zhat = gamma_results[gd,gh][6];
            subp = createSecond(fData,uData,hData,bData,T,groupDict[ci],Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uDict);
            optimize!(subp);
            obj_Dict[gd,gh] = objective_value(subp);              
        end
    end

    gd = 0;
    gh = 0;
    sphat = gamma_results[gd,gh][2];
    sqhat = gamma_results[gd,gh][3];
    xhat = gamma_results[gd,gh][4];
    yhat = gamma_results[gd,gh][5];
    zhat = gamma_results[gd,gh][6];
    subp = createSecond(fData,uData,hData,bData,T,groupDict[ci],Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat,uDict);
    optimize!(subp);
    obj_Dict[gd,gh] = objective_value(subp);              

    return obj_Dict;
end