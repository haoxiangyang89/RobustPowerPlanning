vList = [];
dp_m = dual_master_sep(fData,uData,hData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,xhat,yhat,zhat,sphat,sqhat);
optimize!(dp_m);
u_hp_val = Dict();
u_hm_val = Dict();
u_dp_val = Dict();
u_dm_val = Dict();
for t in 2:T
    for i in hData.hList
        u_hp_val[i,t] = value(dp_m[:u_hp][i,t]);
        u_hm_val[i,t] = value(dp_m[:u_hm][i,t]);
    end
    for m in groupList
        u_dp_val[m,t] = value(dp_m[:u_dp][m,t]);
        u_dm_val[m,t] = value(dp_m[:u_dm][m,t]);
    end
end
V_val = value(dp_m[:V]);


sprob = dual_sub_sep(fData,uData,hData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,
    xhat,yhat,zhat,sphat,sqhat,u_dp_val,u_dm_val,u_hp_val,u_hm_val);
optimize!(sprob);
obj_sp = objective_value(sprob);
d_val_udp = Dict();
d_val_udm = Dict();
d_val_uhp = Dict();
d_val_uhm = Dict();
for t in 2:T
    for m in groupList
        d_val_udp[m,t] = dual(sprob[:udp_constr][m,t]);
        d_val_udm[m,t] = dual(sprob[:udm_constr][m,t]);
    end
    for i in hData.hList
        d_val_uhp[i,t] = dual(sprob[:uhp_constr][i,t]);
        d_val_uhm[i,t] = dual(sprob[:uhm_constr][i,t]);
    end
end
@constraint(dp_m, dp_m[:V] <= obj_sp - sum(d_val_uhp[i,t] * (dp_m[:u_hp][i,t] - u_hp_val[i,t]) + d_val_uhm[i,t] * (dp_m[:u_hm][i,t] - u_hm_val[i,t]) for i in hData.hList for t in 2:T) 
                                            - sum(d_val_udp[m,t] * (dp_m[:u_dp][m,t] - u_dp_val[m,t]) + d_val_udm[m,t] * (dp_m[:u_dm][m,t] - u_dm_val[m,t]) for m in groupList for t in 2:T))

#------------------------------------------------------------------------------------
optimize!(dp_m);
u_hp_val = Dict();
u_hm_val = Dict();
u_dp_val = Dict();
u_dm_val = Dict();
for t in 2:T
    for i in hData.hList
        u_hp_val[i,t] = value(dp_m[:u_hp][i,t]);
        u_hm_val[i,t] = value(dp_m[:u_hm][i,t]);
    end
    for m in groupList
        u_dp_val[m,t] = value(dp_m[:u_dp][m,t]);
        u_dm_val[m,t] = value(dp_m[:u_dm][m,t]);
    end
end
V_val = value(dp_m[:V]);                              

subd = dual_sub_sep_p(fData,uData,hData,T,groupDict,Γ,expansion_factor,vmaxT,vminT,θDmaxT,θDminT,
    xhat,yhat,zhat,sphat,sqhat,u_dp_val,u_dm_val,u_hp_val,u_hm_val);
optimize!(subd);
obj_sp_p = objective_value(subd);
d_val_udp = Dict();
d_val_udm = Dict();
d_val_uhp = Dict();
d_val_uhm = Dict();
for t in 2:T
    for m in groupList
        d_val_udp[m,t] = value(subd[:λeqn_u_dp][m,t]);
        d_val_udm[m,t] = value(subd[:λeqn_u_dm][m,t]);
    end
    for i in hData.hList
        d_val_uhp[i,t] = value(subd[:λeqn_u_hp][i,t]);
        d_val_uhm[i,t] = value(subd[:λeqn_u_hm][i,t]);
    end
end
@constraint(dp_m, dp_m[:V] <= obj_sp_p - sum(d_val_uhp[i,t] * (dp_m[:u_hp][i,t] - u_hp_val[i,t]) + d_val_uhm[i,t] * (dp_m[:u_hm][i,t] - u_hm_val[i,t]) for i in hData.hList for t in 2:T) 
                                            - sum(d_val_udp[m,t] * (dp_m[:u_dp][m,t] - u_dp_val[m,t]) + d_val_udm[m,t] * (dp_m[:u_dm][m,t] - u_dm_val[m,t]) for m in groupList for t in 2:T))

push!(vList, (V_val, obj_sp_p));