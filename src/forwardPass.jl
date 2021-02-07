# forward pass of the SDDP algorithm
function noDisruptionBuild(Δt, T, qpopt = false, solveOpt = true)
    # precalculate data
    Rdict = Dict();
    Xdict = Dict();
    for k in fData.brList
        Rdict[k] = fData.g[k]/(fData.g[k]^2 + fData.b[k]^2);
        Xdict[k] = -fData.b[k]/(fData.g[k]^2 + fData.b[k]^2);
    end
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    # construct the first stage without disruption occurring
    if qpopt
        mp = Model(solver = IpoptSolver(print_level = 0, linear_solver = "ma27"));
    else
        mp = Model(solver = GurobiSolver(GUROBI_ENV,OutputFlag = 0,Threads = 1));
    end

    # set up the variables
    @variable(mp, fData.Pmin[i] <= sp[i in fData.genIDList, t in 1:T] <= fData.Pmax[i]);
    @variable(mp, fData.Qmin[i] <= sq[i in fData.genIDList, t in 1:T] <= fData.Qmax[i]);
    sphatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end
    end
    sqhatsum = Dict();
    for t in 1:T
        for i in fData.IDList
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in 1:T]);
    @variable(mp, q[k in fData.brList, t in 1:T]);
    @variable(mp, fData.Vmin[i]^2 <= v[i in fData.IDList, t in 1:T] <= fData.Vmax[i]^2);
    @variable(mp, 0 <= w[i in bData.IDList, t in 0:T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in 1:T]);
    @variable(mp, zp[i in bData.IDList, t in 1:T]);
    @variable(mp, zq[i in bData.IDList, t in 1:T]);
    @variable(mp, lpp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in 1:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in 2:T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in 1:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in 1:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, pequal[k in fData.brList, t in 1:T], p[k,t] == -p[(k[2],k[1],k[3]),t]);
    @constraint(mp, qequal[k in fData.brList, t in 1:T], q[k,t] == -q[(k[2],k[1],k[3]),t]);
    if qpopt
        @constraint(mp, lineThermal[k in fData.brList, t in 1:T], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2);
        @constraint(mp, bThermal[i in bData.IDList, t in 1:T], zp[i,t]^2+zq[i,t]^2 <= u[i]^2);
    else
        @constraint(mp, lineThermal[k in fData.brList, t in 1:T], norm([p[k,t],q[k,t]]) <= fData.rateA[k]);
        @constraint(mp, bThermal[i in bData.IDList, t in 1:T], norm([zp[i,t],zq[i,t]]) <= u[i]);
    end
    @constraint(mp, powerflow[k in fData.brList, t in 1:T], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in 2:T], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in 1:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in 1:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,0] == bData.bInv[i]);

    # set up the cuts, here tp is the disruption time
    for tp in 2:T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in 1:T]);
    if qpopt
        @constraint(mp, genCost[i in fData.genIDList, t in 1:T], fs[i,t] ==
            sum(fData.cp[i].params[nco]*sp[i,t]^(fData.cp[i].n - nco) for nco in 1:(fData.cp[i].n - 1)));
    else
        @variable(mp,tAux1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3] >= 0);
        @variable(mp,tAux2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3]);
        @variable(mp,tAux3[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3]);

        @constraint(mp,gcAux1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux1[i,t] == fs[i,t] +
            (fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
        @constraint(mp,gcAux2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux2[i,t] == fs[i,t] +
            (-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
        @constraint(mp,gcAux3[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3], tAux3[i,t] == sqrt(fData.cp[i].params[1])*sp[i,t] +
            fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1])));
        @constraint(mp, genCost1[i in fData.genIDList, t in 1:T; fData.cp[i].n == 3],
            norm([tAux2[i,t],tAux3[i,t]]) <= tAux1[i,t]);
        @constraint(mp, genCost2[i in fData.genIDList, t in 1:T; fData.cp[i].n == 2],
            fs[i,t] == fData.cp[i].params[1]*sp[i,t]);
    end

    fDict,lDict,θDict = calDualC1(T, fData, pDistr);
    objExpr = @expression(mp, sum(bData.cost[i]*u[i] for i in bData.IDList) +
        sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in 1:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in 2:T));

    @objective(mp, Min, objExpr);

    if solveOpt
        # solve the problem
        # optimize!(mp, with_optimizer(Gurobi.Optimizer, GUROBI_ENV, OutputFlag = 0,
        #     QCPDual = 1, NumericFocus = 3, BarQCPConvTol = 1e-9, FeasibilityTol = 1e-9));
        statusMp = solve(mp);
        if statusMp != :Optimal
            mp.solver = GurobiSolver(GUROBI_ENV,OutputFlag = 0,NumericFocus = 3,Threads = 1);
            statusMp = solve(mp);
        end
        #with_optimizer(Ipopt.Optimizer, linear_solver = "ma27", acceptable_tol = 1e-8, print_level = 0, max_iter = 10000));
        mpObj = getobjectivevalue(mp);
        println("First stage, solving status $(statusMp)");
        # obtain the solutions
        solSp = Dict();
        solSq = Dict();
        solw = Dict();
        solu = Dict();
        solLp = Dict();
        solLq = Dict();
        solzp = Dict();
        for i in fData.genIDList
            for t in 1:T
                if abs(getvalue(sp[i,t])) > 1e-5
                    solSp[i,t] = getvalue(sp[i,t]);
                else
                    solSp[i,t] = 0;
                end
                if abs(getvalue(sq[i,t])) > 1e-5
                    solSq[i,t] = getvalue(sq[i,t]);
                else
                    solSq[i,t] = 0;
                end
            end
        end
        for i in bData.IDList
            if abs(getvalue(u[i])) > 1e-5
                solu[i] = getvalue(u[i]);
                for t in 1:T
                    if abs(getvalue(w[i,t])) > 1e-5
                        solw[i,t] = getvalue(w[i,t]);
                    else
                        solw[i,t] = 0;
                    end
                    if abs(getvalue(zp[i,t])) > 1e-6
                        solzp[i,t] = getvalue(zp[i,t]);
                    else
                        solzp[i,t] = 0;
                    end
                end
            else
                solu[i] = 0;
                for t in 1:T
                    solw[i,t] = 0;
                    if abs(getvalue(zp[i,t])) > 1e-6
                        solzp[i,t] = getvalue(zp[i,t]);
                    else
                        solzp[i,t] = 0;
                    end
                end
            end
        end
        for i in fData.IDList
            for t in 1:T
                if (abs(getvalue(lpp[i,t])) > 1e-8)|(abs(getvalue(lpm[i,t])) > 1e-8)
                    solLp[i,t] = getvalue(lpp[i,t]) - getvalue(lpm[i,t]);
                else
                    solLp[i,t] = 0;
                end
                if (abs(getvalue(lqp[i,t])) > 1e-8)|(abs(getvalue(lqm[i,t])) > 1e-8)
                    solLq[i,t] = getvalue(lqp[i,t]) - getvalue(lqm[i,t]);
                else
                    solLq[i,t] = 0;
                end
            end
        end

        sol = solData(solSp,solSq,solw,solu,solLp,solLq,solzp);
        return sol,mpObj;
    else
        return mp;
    end
end

function fBuild(td, ωd, currentSol, τ, Δt, T, qpopt = false, solveOpt = true, hardened = [])
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
            if (((k[1],k[2]) in ωd)|((k[2],k[1]) in ωd))&(t <= td + τ)
                if ((k[1],k[2]) in hardened) | ((k[2],k[1]) in hardened)
                    Bparams[k,t] = 1;
                else
                    Bparams[k,t] = 0;
                end
            else
                Bparams[k,t] = 1;
            end
        end
        for i in fData.genIDList
            if (i in ωd)&(t <= td + τ)
                if i in hardened
                    Bparams[i,t] = 1;
                else
                    Bparams[i,t] = 0;
                end
            else
                Bparams[i,t] = 1;
            end
        end
    end
    for i in fData.genIDList
        Bparams[i,td - 1] = 1;
    end

    if qpopt
        mp = Model(solver = IpoptSolver(print_level = 0, linear_solver = "ma27"));
    else
        mp = Model(solver = GurobiSolver(GUROBI_ENV,OutputFlag = 0,Threads = 1));
    end

    # set up the variables
    @variable(mp, fData.Pmin[i]*Bparams[i,t] <= sp[i in fData.genIDList,t in (td - 1):T] <= fData.Pmax[i]*Bparams[i,t]);
    @variable(mp, fData.Qmin[i]*Bparams[i,t] <= sq[i in fData.genIDList,t in td:T] <= fData.Qmax[i]*Bparams[i,t]);
    sphatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sphatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sphatsum[i,t] += sp[j,t];
                end
            end
        end
    end
    sqhatsum = Dict();
    for t in td:T
        for i in fData.IDList
            sqhatsum[i,t] = @expression(mp,0.0);
            if i in keys(fData.LocRev)
                for j in fData.LocRev[i]
                    sqhatsum[i,t] += sq[j,t];
                end
            end
        end
    end

    @variable(mp, p[k in fData.brList, t in td:T]);
    @variable(mp, q[k in fData.brList, t in td:T]);
    @variable(mp, fData.Vmin[i]^2 <= v[i in fData.IDList, t in td:T] <= fData.Vmax[i]^2);
    @variable(mp, 0 <= w[i in bData.IDList, t in (td - 1):T] <= bData.cap[i]);
    @variable(mp, y[i in bData.IDList, t in td:T]);
    @variable(mp, zp[i in bData.IDList, t in td:T]);
    @variable(mp, zq[i in bData.IDList, t in td:T]);
    @variable(mp, lpp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqp[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lpm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, lqm[i in fData.IDList, t in td:T] >= 0);
    @variable(mp, 0 <= u[i in bData.IDList] <= bData.uCap[i]);
    @variable(mp, θ[tp in (td + τ + 1):T, ω in Ω] >= 0);

    # set up the constraints
    @constraint(mp, pBalance[i in fData.IDList, t in td:T], sum(zp[b,t] for b in bData.IDList if bData.Loc[b] == i) + lpp[i,t] - lpm[i,t] +
        sphatsum[i,t] - dData.pd[i][t] == sum(p[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, qBalance[i in fData.IDList, t in td:T], sum(zq[b,t] for b in bData.IDList if bData.Loc[b] == i) + lqp[i,t] - lqm[i,t] +
        sqhatsum[i,t] - dData.qd[i][t] == sum(q[k,t] for k in fData.branchDict1[i]));
    @constraint(mp, pequal[k in fData.brList, t in td:T], p[k,t] == -p[(k[2],k[1],k[3]),t]);
    @constraint(mp, qequal[k in fData.brList, t in td:T], q[k,t] == -q[(k[2],k[1],k[3]),t]);
    if qpopt
        @constraint(mp, lineThermal1[k in fData.brList, t in td:T;Bparams[k,t] == 1], p[k,t]^2+q[k,t]^2 <= fData.rateA[k]^2);
        @constraint(mp, bThermal[i in bData.IDList, t in td:T], zp[i,t]^2+zq[i,t]^2 <= u[i]^2);
    else
        @constraint(mp, lineThermal1[k in fData.brList, t in td:T;Bparams[k,t] == 1], norm([p[k,t],q[k,t]]) <= fData.rateA[k]);
        @constraint(mp, bThermal[i in bData.IDList, t in td:T], norm([zp[i,t],zq[i,t]]) <= u[i]);
    end
    @constraint(mp, lineThermal2[k in fData.brList, t in td:T;Bparams[k,t] == 0], p[k,t] == 0);
    @constraint(mp, lineThermal3[k in fData.brList, t in td:T;Bparams[k,t] == 0], q[k,t] == 0);
    @constraint(mp, powerflow1[k in fData.brList, t in td:T;Bparams[k,t] == 1], v[k[2],t] == v[k[1],t] - 2*(Rdict[k]*p[k,t] + Xdict[k]*q[k,t]));
    @constraint(mp, rampUp[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] <= fData.RU[i]);
    @constraint(mp, rampDown[i in fData.genIDList, t in td:T; Bparams[i,t] == 1], sp[i,t] - sp[i,t - 1] >= fData.RD[i]);
    @constraint(mp, bInv[i in bData.IDList, t in td:T], w[i,t] == w[i,t-1] - y[i,t]*Δt);
    # @constraint(mp, bThermal1[i in bData.IDList, t in td:T], zp[i,t] <= u[i]);
    # @constraint(mp, bThermal3[i in bData.IDList, t in td:T], -zp[i,t] <= u[i]);
    # #@constraint(mp, bThermal2[i in bData.IDList, t in 1:T], zq[i,t] <= u[i]);
    # @constraint(mp, bThermal2[i in bData.IDList, t in td:T], zq[i,t] == 0);
    @constraint(mp, bEfficient[i in bData.IDList, l in 1:length(bData.ηα[i]), t in td:T], zp[i,t] <= bData.ηα[i][l]*y[i,t] + bData.ηβ[i][l]);
    @constraint(mp, bInvIni[i in bData.IDList], w[i,td - 1] == currentSol.w[i,td - 1]);
    @constraint(mp, spIni[i in fData.genIDList], sp[i,td - 1] == currentSol.sp[i,td - 1]);
    @constraint(mp, uIni[i in bData.IDList], u[i] == currentSol.u[i]);

    # set up the cuts, here tp is the disruption time
    for tp in (td + τ + 1):T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    @constraint(mp, θ[tp,ω] >= cutDict[tp,ω][l].vhat +
                        sum(cutDict[tp,ω][l].λ[i]*sp[i,tp - 1] for i in fData.genIDList) +
                        sum(cutDict[tp,ω][l].γ[i]*w[i,tp - 1] + cutDict[tp,ω][l].μ[i]*u[i] for i in bData.IDList));
                end
            end
        end
    end

    # set up the objective function
    @variable(mp,fs[i in fData.genIDList, t in td:T]);
    if qpopt
        @constraint(mp, genCost[i in fData.genIDList, t in td:T], fs[i,t] ==
            sum(fData.cp[i].params[nco]*sp[i,t]^(fData.cp[i].n - nco) for nco in 1:(fData.cp[i].n - 1)));
    else
        @variable(mp,tAux1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3] >= 0);
        @variable(mp,tAux2[i in fData.genIDList, t in td:T; fData.cp[i].n == 3]);
        @variable(mp,tAux3[i in fData.genIDList, t in td:T; fData.cp[i].n == 3]);

        @constraint(mp,gcAux1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux1[i,t] == fs[i,t] +
            (fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
        @constraint(mp,gcAux2[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux2[i,t] == fs[i,t] +
            (-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
        @constraint(mp,gcAux3[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux3[i,t] == sqrt(fData.cp[i].params[1])*sp[i,t] +
            fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1])));
        @constraint(mp, genCost1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3],
            norm([tAux2[i,t],tAux3[i,t]]) <= tAux1[i,t]);
        @constraint(mp, genCost2[i in fData.genIDList, t in td:T; fData.cp[i].n == 2],
            fs[i,t] == fData.cp[i].params[1]*sp[i,t]);
    end

    fDict,lDict,θDict = calDualC(td, τ, T, fData, pDistr);
    objExpr = @expression(mp, sum(sum(fDict[i,t]*fs[i,t] for i in fData.genIDList) +
        sum(lDict[i,t]*(lpp[i,t] + lpm[i,t] + lqp[i,t] + lqm[i,t]) for i in fData.IDList) for t in td:T) +
        sum(sum(θDict[t,ω]*θ[t,ω] for ω in Ω) for t in (td + τ + 1):T));
    @objective(mp, Min, objExpr);

    if solveOpt
        # solve the problem
        # optimize!(mp, with_optimizer(Gurobi.Optimizer, GUROBI_ENV, OutputFlag = 0,
        #     QCPDual = 1, NumericFocus = 3, BarQCPConvTol = 1e-9, FeasibilityTol = 1e-9));
        statusMp = solve(mp);
        if statusMp != :Optimal
            mp.solver = GurobiSolver(GUROBI_ENV,OutputFlag = 0,NumericFocus = 3,Threads = 1);
            statusMp = solve(mp);
        end
        #optimize!(mp, with_optimizer(Ipopt.Optimizer, linear_solver = "ma27", acceptable_tol = 1e-8, print_level = 0, max_iter = 10000));
        mpObj = getobjectivevalue(mp);
        println("Disruption time $(td), scenario $(ωd), solving status $(statusMp)");
        # obtain the solutions
        solSp = Dict();
        solSq = Dict();
        solw = Dict();
        solu = Dict();
        solLp = Dict();
        solLq = Dict();
        solzp = Dict();
        for i in fData.genIDList
            for t in td:T
                if abs(getvalue(sp[i,t])) > 1e-5
                    solSp[i,t] = getvalue(sp[i,t]);
                else
                    solSp[i,t] = 0;
                end
                if abs(getvalue(sq[i,t])) > 1e-5
                    solSq[i,t] = getvalue(sq[i,t]);
                else
                    solSq[i,t] = 0;
                end
            end
        end
        for i in bData.IDList
            if abs(getvalue(u[i])) > 1e-5
                solu[i] = getvalue(u[i]);
                for t in td:T
                    if abs(getvalue(w[i,t])) > 1e-5
                        solw[i,t] = getvalue(w[i,t]);
                    else
                        solw[i,t] = 0;
                    end
                    if abs(getvalue(zp[i,t])) > 1e-6
                        solzp[i,t] = getvalue(zp[i,t]);
                    else
                        solzp[i,t] = 0;
                    end
                end
            else
                solu[i] = 0;
                for t in td:T
                    solw[i,t] = 0;
                    if abs(getvalue(zp[i,t])) > 1e-6
                        solzp[i,t] = getvalue(zp[i,t]);
                    else
                        solzp[i,t] = 0;
                    end
                end
            end
        end
        for i in fData.IDList
            for t in td:T
                if (abs(getvalue(lpp[i,t])) > 1e-8)|(abs(getvalue(lpm[i,t])) > 1e-8)
                    solLp[i,t] = getvalue(lpp[i,t]) - getvalue(lpm[i,t]);
                else
                    solLp[i,t] = 0;
                end
                if (abs(getvalue(lqp[i,t])) > 1e-8)|(abs(getvalue(lqm[i,t]))> 1e-8)
                    solLq[i,t] = getvalue(lqp[i,t]) - getvalue(lqm[i,t]);
                else
                    solLq[i,t] = 0;
                end
            end
        end

        sol = solData(solSp,solSq,solw,solu,solLp,solLq,solzp);
        return sol,mpObj;
    else
        return mp;
    end
end

function constructForwardM(td, ωd, sol, τ, Δt, T, qpopt = false, hardenend = [])
    # construct the math program given the state variables and current stage
    if td == 1
        # if it is the no-disruption problem
        sol,objV = noDisruptionBuild(Δt, T, qpopt);
    else
        # if it is f_{ht}^ω
        sol,objV = fBuild(td, ωd, sol, τ, Δt, T, qpopt, true, hardenend);
    end
    return sol,objV;
end

function buildPath(T, Δt, qpopt = false, pathList = [], hardened = [])
    disT = 1;
    ωd = 0;
    τω = 0;
    costn = 0;
    solHist = [];
    currentLB = 0;
    currentSol = solData(Dict(),Dict(),Dict(),Dict(),Dict(),Dict(),Dict());
    iter = 1;
    while disT <= T
        # solve the current stage problem, state variables are passed
        nowT = disT;
        currentSol,objV = constructForwardM(disT, ωd, currentSol, τω, Δt, T, qpopt, hardened);
        push!(solHist,(currentSol,nowT,ωd,τω));

        # generate disruption
        if pathList == []
            tp,ωd,τω = genScenario(pDistr);
        else
            tp,ωd,τω = pathList[iter];
        end
        iter += 1;
        if nowT == 1
            currentLB = objV;
            disT += tp;
            disT = min(disT, T + 1);
            # calculate the cost of the solution until the next disruption time
            costn += sum(currentSol.u[i]*bData.cost[i] for i in bData.IDList);
            costn = calCostF(costn, currentSol, T, fData, nowT, disT);
        else
            disT += tp + τω;
            disT = min(disT, T + 1);
            # calculate the cost of the solution until the next disruption time
            costn = calCostF(costn, currentSol, T, fData, nowT, disT);
        end
    end
    println("Path Built!");
    return [solHist,currentLB,costn];
end

function exeForward(T, Δt, N, qpopt = false, pathDict = Dict(), hardened = [])
    # execution of forward pass
    # input: N: the number of trial points;
    #       cutDict: set of currently generated cuts (global in every core)
    # output: solList: a list of solution paths
    solDict = Dict();
    costDict = Dict();
    objV = 0;
    currentLB = 0;
    # for n in 1:N
    #     # for each trial path
    #     returnData = buildPath(τ, T, Δt, fData, bData, dData, pDistr, cutDict);
    #     solDict[n] = returnData[1];
    #     costDict[n] = returnData[2];
    # end
    if pathDict == Dict()
        for i in 1:N
            pathDict[i] = [];
        end
    end
    returnData = pmap(i -> buildPath(T, Δt, qpopt, pathDict[i], hardened), 1:N);
    for n in 1:N
        solDict[n] = returnData[n][1];
        costDict[n] = returnData[n][3];
    end
    currentLB = returnData[1][2];
    return solDict, currentLB, costDict;
end

function exeForward_simuOpt(T, Δt, N, iterNo, qpopt = false, hardened = [])
    # execution of forward pass, with samples selected to cover more time
    # input: N: the number of trial points;
    #       cutDict: set of currently generated cuts (global in every core)
    # output: solList: a list of solution paths
    solDict = Dict();
    costDict = Dict();
    objV = 0;
    currentLB = 0;
    genCutsCount = Dict();
    for t in 2:T
        genCutsCount[t] = 0;
    end
    pathDict,pathTimeList = pathSimu_cover(N, T, pDistr, genCutsCount, iterNo);
    for i in 1:N
        # update genCutsCount
        for j in 1:length(pathTimeList[i])
            genCutsCount[pathTimeList[i][j]] += 1;
        end
    end
    returnData = pmap(i -> buildPath(T, Δt, qpopt, pathDict[i], hardened), 1:N);
    for n in 1:N
        solDict[n] = returnData[n][1];
        costDict[n] = returnData[n][3];
    end
    currentLB = returnData[1][2];
    return solDict, currentLB, costDict;
end

function exeForward_last(T, Δt, N, qpopt = false, hardened = [])
    # execution of forward pass, with samples selected to cover more time
    # input: N: the number of trial points;
    #       cutDict: set of currently generated cuts (global in every core)
    # output: solList: a list of solution paths
    solDict = Dict();
    costDict = Dict();
    objV = 0;
    currentLB = 0;
    pathDict,pathTimeList = pathSimu_cover_last(N, T, pDistr);
    returnData = pmap(i -> buildPath(T, Δt, qpopt, pathDict[i], hardened), 1:N);
    for n in 1:N
        solDict[n] = returnData[n][1];
        costDict[n] = returnData[n][3];
    end
    currentLB = returnData[1][2];
    return solDict, currentLB, costDict;
end
