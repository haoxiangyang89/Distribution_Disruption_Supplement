@everywhere function extForm(td, ωd, inheritData, baseProb, τ, Δt, T, fData, bData, dData, pDistr,coneBool = false)
    # extensive formulation could not have variable/constraint names
    # inheritData: [1]: sp, [2]: w, [3]: u (only contains the information for the linking time period)

    println("========= Disruption time $(td), scenario $(ωd) modeling =========");
    # precalculate data
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];
    Bparams = Dict();
    for t in td:T
        # create B parameters
        for k in fData.brList
            # if the line is disrupted and it is within disruption time
            if (((k[1],k[2]) == ωd)|((k[2],k[1]) == ωd))&(t <= td + τ)
                Bparams[k,t] = 0;
            else
                Bparams[k,t] = 1;
            end
        end
        for i in fData.genIDList
            if (i == ωd)&(t <= td + τ)
                Bparams[i,t] = 0;
            else
                Bparams[i,t] = 1;
            end
        end
    end

    # set up the variables
    spDict = Dict();
    sqDict = Dict();
    for i in fData.genIDList
        for t in (td - 1):T
            spDict[i,t] = @variable(mExt, lowerbound = fData.Pmin[i],  upperbound = fData.Pmax[i], basename="sp_$(td)_$(i)_$(t)");
            sqDict[i,t] = @variable(mExt, lowerbound = fData.Qmin[i],  upperbound = fData.Qmax[i], basename="sq_$(td)_$(i)_$(t)");
        end
    end
    sphatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mExt,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += spDict[j,t];
                end
            end
        end
    end
    sqhatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sqhatsum[i,t] = @expression(mExt,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sqDict[j,t];
                end
            end
        end
    end

    pDict = Dict();
    qDict = Dict();
    vDict = Dict();
    wDict = Dict();
    yDict = Dict();
    zpDict = Dict();
    zqDict = Dict();
    lppDict = Dict();
    lqpDict = Dict();
    lpmDict = Dict();
    lqmDict = Dict();
    uDict = Dict();
    fsDict = Dict();
    for k in fData.brList
        for t in td:T
            pDict[k,t] = @variable(mExt, basename="p_$(td)_$(k)_$(t)");
            qDict[k,t] = @variable(mExt, basename="q_$(td)_$(k)_$(t)");
        end
    end
    for i in fData.IDList
        for t in td:T
            vDict[i,t] = @variable(mExt, lowerbound = fData.Vmin[i]^2, upperbound = fData.Vmax[i]^2, basename="v_$td");
            lppDict[i,t] = @variable(mExt, lowerbound = 0, basename="lpp_$(td)_$(i)_$(t)");
            lqpDict[i,t] = @variable(mExt, lowerbound = 0, basename="lqp_$(td)_$(i)_$(t)");
            lpmDict[i,t] = @variable(mExt, lowerbound = 0, basename="lpm_$(td)_$(i)_$(t)");
            lqmDict[i,t] = @variable(mExt, lowerbound = 0, basename="lqm_$(td)_$(i)_$(t)");
        end
    end
    for i in bData.IDList
        wDict[i,td - 1] = @variable(mExt, lowerbound = 0, upperbound = bData.cap[i], basename="w_$(td)_$(i)_$(td - 1)");
        uDict[i] = @variable(mExt, lowerbound = 0, upperbound = bData.uCap[i], basename="u_$(td)_$(i)");
        for t in td:T
            wDict[i,t] = @variable(mExt, lowerbound = 0, upperbound = bData.cap[i], basename="w_$(td)_$(i)_$(t)");
            yDict[i,t] = @variable(mExt, basename="y_$(td)_$(i)_$(t)");
            zpDict[i,t] = @variable(mExt, basename="zp_$(td)_$(i)_$(t)");
            zqDict[i,t] = @variable(mExt, basename="zq_$(td)_$(i)_$(t)");
        end
    end

    # set up the constraints
    for i in fData.IDList
        for t in td:T
            @constraint(mExt, sum(zpDict[b,t] for b in bData.IDList if bData.Loc[b] == i) + lppDict[i,t] - lpmDict[i,t] +
                sphatsum[i,t] - dData.pd[i][t] == sum(pDict[k,t] for k in fData.branchDict1[i]));
            @constraint(mExt, sum(zqDict[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqpDict[i,t] - lqmDict[i,t] +
                sqhatsum[i,t] - dData.qd[i][t] == sum(qDict[k,t] for k in fData.branchDict1[i]));
        end
    end
    for i in fData.genIDList
        if td != 1
            @constraint(mExt, spDict[i,td - 1] == inheritData[1][i]);
        end
        for t in td:T
            if (t != 1)&(Bparams[i,t] == 1)
                @constraint(mExt, spDict[i,t] - spDict[i,t - 1] <= fData.RU[i]);
                @constraint(mExt, spDict[i,t] - spDict[i,t - 1] >= fData.RD[i]);
            end
            if Bparams[i,t] == 0
                @constraint(mExt, spDict[i,t] == 0);
                @constraint(mExt, sqDict[i,t] == 0);
            end
        end
    end
    for k in fData.brList
        for t in td:T
            @constraint(mExt, pDict[k,t] == -pDict[(k[2],k[1],k[3]),t]);
            @constraint(mExt, qDict[k,t] == -qDict[(k[2],k[1],k[3]),t]);
            if Bparams[k,t] == 1
                if coneBool
                    @constraint(mExt, norm([pDict[k,t],qDict[k,t]]) <= fData.rateA[k]);
                else
                    @constraint(mExt, pDict[k,t]^2 + qDict[k,t]^2 <= fData.rateA[k]^2);
                end
                @constraint(mExt, vDict[k[2],t] == vDict[k[1],t] - 2*(Rdict[k]*pDict[k,t] + Xdict[k]*qDict[k,t]));
            else
                @constraint(mExt, pDict[k,t] == 0);
                @constraint(mExt, qDict[k,t] == 0);
            end
        end
    end
    for i in bData.IDList
        if td != 1
            @constraint(mExt, uDict[i] == inheritData[3][i]);
        end
        @constraint(mExt, wDict[i,td - 1] == inheritData[2][i]);
        for t in td:T
            @constraint(mExt, wDict[i,t] == wDict[i,t - 1] - yDict[i,t]*Δt);
            if coneBool
                @constraint(mExt, norm([zpDict[i,t],zqDict[i,t]]) <= uDict[i]);
            else
                @constraint(mExt, zpDict[i,t]^2 + zqDict[i,t]^2 <= uDict[i]^2);
            end
            for l in 1:length(bData.ηα[i])
                @constraint(mExt, zpDict[i,t] <= bData.ηα[i][l]*yDict[i,t] + bData.ηβ[i][l]);
            end
            @constraint(mExt, wDict[i,t] <= bData.cap[i]);
        end
    end

    # recursion through the possible scenario
    tList = sort([t for t in keys(pDistr.tDistrn)]);
    objExpr = getobjective(mExt);
    if td == 1
        objExpr += sum(bData.cost[i]*uDict[i] for i in bData.IDList);
        # only the first pass will execute this
    end
    for tp in tList
        if td == 1
            if td + tp > T
                # if the next disruption is over the time horizon
                dExpr = fData.cz*sum(sum(lppDict[i,t] + lqpDict[i,t] + lpmDict[i,t] + lqmDict[i,t] for i in fData.IDList) for t in td:T);
                for t in td:T
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            if coneBool
                                fsDict[i,t] = @variable(mExt, lowerbound = 0, basename = "fs_$(td)_$(i)_$(t)");
                                dExpr += fData.cp[i].params[1]*fsDict[i,t] + fData.cp[i].params[2]*spDict[i,t];
                                @constraint(mExt, norm([spDict[i,t],fsDict[i,t] - 1/4]) <= fsDict[i,t] + 1/4);
                            else
                                dExpr += fData.cp[i].params[1]*(spDict[i,t]^2) + fData.cp[i].params[2]*spDict[i,t];
                            end
                        elseif fData.cp[i].n == 2
                            dExpr += fData.cp[i].params[1]*spDict[i,t];
                        end
                    end
                end
                objExpr += baseProb*pDistr.tDistrn[tp]*dExpr;
                @objective(mExt, Min, objExpr);
            else
                # if not
                dExpr = fData.cz*sum(sum(lppDict[i,t] + lqpDict[i,t] + lpmDict[i,t] + lqmDict[i,t] for i in fData.IDList) for t in td:(td + tp - 1));
                for t in td:(td + tp - 1)
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            if coneBool
                                fsDict[i,t] = @variable(mExt, lowerbound = 0, basename = "fs_$(td)_$(i)_$(t)");
                                dExpr += fData.cp[i].params[1]*fsDict[i,t] + fData.cp[i].params[2]*spDict[i,t];
                                @constraint(mExt, norm([spDict[i,t],fsDict[i,t] - 1/4]) <= fsDict[i,t] + 1/4);
                            else
                                dExpr += fData.cp[i].params[1]*(spDict[i,t]^2) + fData.cp[i].params[2]*spDict[i,t];
                            end
                        elseif fData.cp[i].n == 2
                            dExpr += fData.cp[i].params[1]*spDict[i,t];
                        end
                    end
                end
                for ω in Ω
                    objExpr += baseProb*pDistr.tDistrn[tp]*pDistr.ωDistrn[ω]*dExpr;
                    spInherit = Dict();
                    wInherit = Dict();
                    uInherit = Dict();
                    for i in fData.genIDList
                        spInherit[i] = spDict[i,td + tp - 1];
                    end
                    for i in bData.IDList
                        wInherit[i] = wDict[i,td + tp - 1];
                        uInherit[i] = uDict[i];
                    end
                    inheritData = [spInherit,wInherit,uInherit];
                    @objective(mExt, Min, objExpr);
                    # println("=========================",1," ",td," ",tp,"=========================");
                    # println(mExt.obj);
                    global mExt = extForm(td + tp, ω, inheritData, baseProb*pDistr.tDistrn[tp]*pDistr.ωDistrn[ω], τ, Δt, T, fData, bData, dData, pDistr, coneBool);
                    # println("=========================",11," ",td," ",tp,"=========================");
                    # println(mExt.obj);
                    objExpr = getobjective(mExt);
                end
            end
        else
            if td + τ + tp > T
                # if the next disruption is over the time horizon
                # add to the objective function
                dExpr = fData.cz*sum(sum(lppDict[i,t] + lqpDict[i,t] + lpmDict[i,t] + lqmDict[i,t] for i in fData.IDList) for t in td:T);
                for t in td:T
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            if coneBool
                                fsDict[i,t] = @variable(mExt, lowerbound = 0, basename = "fs_$(td)_$(i)_$(t)");
                                dExpr += fData.cp[i].params[1]*fsDict[i,t] + fData.cp[i].params[2]*spDict[i,t];
                                @constraint(mExt, norm([spDict[i,t],fsDict[i,t] - 1/4]) <= fsDict[i,t] + 1/4);
                            else
                                dExpr += fData.cp[i].params[1]*(spDict[i,t]^2) + fData.cp[i].params[2]*spDict[i,t];
                            end
                        elseif fData.cp[i].n == 2
                            dExpr += fData.cp[i].params[1]*spDict[i,t];
                        end
                    end
                end
                objExpr += baseProb*pDistr.tDistrn[tp]*dExpr;
                @objective(mExt, Min, objExpr);
            else
                # if not
                # add the current part to the objective function
                # recursion to the next disruption
                dExpr = fData.cz*sum(sum(lppDict[i,t] + lqpDict[i,t] + lpmDict[i,t] + lqmDict[i,t] for i in fData.IDList) for t in td:(td + tp + τ - 1));
                for t in td:(td + tp + τ - 1)
                    for i in fData.genIDList
                        # add generator cost
                        if fData.cp[i].n == 3
                            if coneBool
                                fsDict[i,t] = @variable(mExt, lowerbound = 0, basename = "fs_$(td)_$(i)_$(t)");
                                dExpr += fData.cp[i].params[1]*fsDict[i,t] + fData.cp[i].params[2]*spDict[i,t];
                                @constraint(mExt, norm([spDict[i,t],fsDict[i,t] - 1/4]) <= fsDict[i,t] + 1/4);
                            else
                                dExpr += fData.cp[i].params[1]*(spDict[i,t]^2) + fData.cp[i].params[2]*spDict[i,t];
                            end
                        elseif fData.cp[i].n == 2
                            dExpr += fData.cp[i].params[1]*spDict[i,t];
                        end
                    end
                end
                for ω in Ω
                    objExpr += baseProb*pDistr.tDistrn[tp]*pDistr.ωDistrn[ω]*dExpr;
                    spInherit = Dict();
                    wInherit = Dict();
                    uInherit = Dict();
                    for i in fData.genIDList
                        spInherit[i] = spDict[i,td + tp + τ - 1];
                    end
                    for i in bData.IDList
                        wInherit[i] = wDict[i,td + tp + τ - 1];
                        uInherit[i] = uDict[i];
                    end
                    inheritData = [spInherit,wInherit,uInherit];
                    @objective(mExt, Min, objExpr);
                    # println("=========================",2," ",td," ",tp,"=========================");
                    # println(mExt.obj);
                    global mExt = extForm(td + tp + τ, ω, inheritData, baseProb*pDistr.tDistrn[tp]*pDistr.ωDistrn[ω], τ, Δt, T, fData, bData, dData, pDistr, coneBool);
                    # println("=========================",22," ",td," ",tp,"=========================");
                    # println(mExt.obj);
                    objExpr = getobjective(mExt);
                end
            end
        end
    end
    return mExt;
end
