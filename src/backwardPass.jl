# backward pass of the SDDP algorithm
function fBuild_D(td, ωd, currentPath, τ, Δt, T, qpopt = false, solveOpt = true, hardened = [])
    prevtpInd = maximum([i for i in 1:length(currentPath) if currentPath[i][2] < td]);
    currentSol = currentPath[prevtpInd][1];
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
        mp = Model(solver = GurobiSolver(GUROBI_ENV, OutputFlag = 0, QCPDual = 1,Threads = 1,NumericFocus = 3));
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
        @constraint(mp, lineThermal1[k in fData.brList, t in td:T;Bparams[k,t] == 1], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2);
        @constraint(mp, bThermal[i in bData.IDList, t in td:T], zp[i,t]^2 + zq[i,t]^2 <= u[i]^2);
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
    @constraint(mp, bThermal1[i in bData.IDList, t in td:T], zp[i,t] <= u[i]);
    @constraint(mp, bThermal2[i in bData.IDList, t in td:T], zq[i,t] <= u[i]);
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
    objExpr = calObj(mp, td, τ, T, fData, bData, pDistr, sp, lpp, lqp, lpm, lqm, θ);
    @objective(mp, Min, objExpr);

    if solveOpt
        # solve the problem
        # optimize!(mp, with_optimizer(Gurobi.Optimizer, GUROBI_ENV, OutputFlag = 0,
        #     QCPDual = 1, NumericFocus = 3, BarQCPConvTol = 1e-9, FeasibilityTol = 1e-9));
        #optimize!(mp, with_optimizer(Ipopt.Optimizer, linear_solver = "ma27", print_level = 0, acceptable_tol = 1e-8, max_iter = 10000));
        statusMp = solve(mp);
        println(statusMp, " ", td, " ", ωd);
        # obtain the primal solutions & obj value
        vhat = getobjectivevalue(mp);
        # obtain the solutions
        solSp = Dict();
        solw = Dict();
        solu = Dict();
        for i in fData.genIDList
            solSp[i] = getvalue(sp[i,td - 1]);
        end
        for i in bData.IDList
            solu[i] = getvalue(u[i]);
            solw[i] = getvalue(w[i,td - 1]);
        end
        # obtain the dual solutions
        dsolλ = Dict();
        dsolγ = Dict();
        dsolμ = Dict();
        for i in fData.genIDList
            dsolλ[i] = getdual(spIni[i]);
            vhat -= dsolλ[i]*solSp[i];
            if abs(dsolλ[i]) < 1e-4
                if dsolλ[i] < 0
                    vhat += dsolλ[i]*fData.Pmax[i];
                end
                dsolλ[i] = 0;
            end
        end
        for i in bData.IDList
            dsolγ[i] = getdual(bInvIni[i]);
            vhat -= dsolγ[i]*solw[i];
            if abs(dsolγ[i]) < 1e-4
                if dsolγ[i] < 0
                    vhat += dsolγ[i]*bData.cap[i];
                end
                dsolγ[i] = 0;
            end
            dsolμ[i] = getdual(uIni[i]);
            vhat -= dsolμ[i]*solu[i];
            if abs(dsolμ[i]) < 1e-4
                if dsolμ[i] < 0
                    vhat += dsolμ[i]*bData.uCap[i];
                end
                dsolμ[i] = 0;
            end
        end
        cutTemp = cutData(statusMp,dsolλ,dsolγ,dsolμ,vhat);
        return cutTemp;
    else
        return mp;
    end
end

function dfBuild_D(td, ωd, currentPath, τ, Δt, T, qpopt = false, solveOpt = true, hardened = [])
    # trying to find the previous solution that can be used for this td
    prevtpInd = maximum([i for i in 1:length(currentPath) if currentPath[i][2] < td]);
    currentSol = currentPath[prevtpInd][1];
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
        dp = Model(solver = IpoptSolver(print_level = 0, linear_solver = "ma27"));
    else
        dp = Model(solver = GurobiSolver(GUROBI_ENV,OutputFlag = 0,Threads = 1));
    end

    # simple bounds dual variables
    @variable(dp, λspu[i in fData.genIDList,t in td:T] >= 0);
    @variable(dp, λspl[i in fData.genIDList,t in td:T] <= 0);
    @variable(dp, λsqu[i in fData.genIDList,t in td:T] >= 0);
    @variable(dp, λsql[i in fData.genIDList,t in td:T] <= 0);
    @variable(dp, λvu[i in fData.IDList,t in td:T] >= 0);
    @variable(dp, λvl[i in fData.IDList,t in td:T] <= 0);
    @variable(dp, λwu[i in bData.IDList,t in td:T] >= 0);
    @variable(dp, λuu[i in bData.IDList] >= 0);

    # cuts dual variables
    λcuts = Dict();
    spCuts = Dict();
    wCuts = Dict();
    for tp in (td-1):T
        spCuts[tp] = [];
        wCuts[tp] = [];
    end
    for tp in (td + τ + 1):T
        for ω in Ω
            if (tp,ω) in keys(cutDict)
                for l in 1:length(cutDict[tp,ω])
                    λcuts[tp,ω,l] = @variable(dp, upperbound = 0);
                    push!(spCuts[tp - 1],(tp,ω,l));
                    push!(wCuts[tp - 1],(tp,ω,l));
                end
            end
        end
    end

    # linear constraints dual variables
    @variable(dp, λpb[i in fData.IDList, t in td:T]);
    @variable(dp, λqb[i in fData.IDList, t in td:T]);
    @variable(dp, λpe[k in fData.brList, t in td:T]);
    @variable(dp, λqe[k in fData.brList, t in td:T]);
    @variable(dp, λlt2[k in fData.brList, t in td:T]);
    @variable(dp, λlt3[k in fData.brList, t in td:T]);
    @variable(dp, λpf[k in fData.brList, t in td:T]);
    @variable(dp, λru[i in fData.genIDList, t in td:T] >= 0);
    @variable(dp, λrd[i in fData.genIDList, t in td:T] <= 0);
    @variable(dp, λbInv[i in bData.IDList, t in td:T]);
    @variable(dp, λbCoeff[i in bData.IDList, l in 1:length(bData.ηα[i]), t in td:T] >= 0);
    @variable(dp, λsp[i in fData.genIDList]);
    @variable(dp, λw[i in bData.IDList]);
    @variable(dp, λu[i in bData.IDList]);
    @variable(dp, λgc2[i in fData.genIDList, t in td:T]);
    @variable(dp, λgcAux1[i in fData.genIDList, t in td:T]);
    @variable(dp, λgcAux2[i in fData.genIDList, t in td:T]);
    @variable(dp, λgcAux3[i in fData.genIDList, t in td:T]);

    # SOCP constraints dual variables
    @variable(dp,μ1[k in fData.brList, t in td:T, j in 1:2]);
    @variable(dp,μ2[i in bData.IDList, t in td:T, j in 1:2]);
    @variable(dp,μ3[i in fData.genIDList, t in td:T, j in 1:2]);
    @variable(dp,ν1[k in fData.brList, t in td:T] >= 0);
    @variable(dp,ν2[i in bData.IDList, t in td:T] >= 0);
    @variable(dp,ν3[i in fData.genIDList, t in td:T] >= 0);

    # dual constraints:
    # sp
    bool3 = Dict();
    for i in fData.genIDList
        if fData.cp[i].n == 3
            bool3[i] = 1;
        else
            bool3[i] = 0;
        end
    end
    fDict,lDict,θDict = calDualC(td, τ, T, fData, pDistr);
    @constraint(dp, spConstraint[i in fData.genIDList, t in td:(T - 1)],
        λspu[i,t] + λspl[i,t] + λpb[fData.Loc[i],t] + λru[i,t]*Bparams[i,t] -
        λru[i,t+1]*Bparams[i,t+1] + λrd[i,t]*Bparams[i,t] - λrd[i,t+1]*Bparams[i,t+1] -
        λgcAux3[i,t]*sqrt(fData.cp[i].params[1])*bool3[i] - λgc2[i,t]*fData.cp[i].params[1]*(1 - bool3[i]) -
        sum(λcuts[item]*cutDict[(item[1],item[2])][item[3]].λ[i] for item in spCuts[t] if t >= td + τ) == 0);
    @constraint(dp, spConstraintT[i in fData.genIDList],
        λspu[i,T] + λspl[i,T] + λpb[fData.Loc[i],T] + λru[i,T]*Bparams[i,T] + λrd[i,T]*Bparams[i,T] -
        λgcAux3[i,T]*sqrt(fData.cp[i].params[1])*bool3[i] - λgc2[i,T]*fData.cp[i].params[1]*(1 - bool3[i]) == 0);# for all time period after td
    @constraint(dp, spConstraint1[i in fData.genIDList], -λru[i,td]*Bparams[i,td] - λrd[i,td]*Bparams[i,td] + λsp[i] == 0); # for td-1
    # sq
    @constraint(dp, sqConstraint[i in fData.genIDList, t in td:T],
        λsqu[i,t] + λsql[i,t] + λqb[fData.Loc[i],t] == 0);
    # p
    @constraint(dp, pConstraint[k in fData.brList, t in td:T], -λpb[k[1],t] + λpe[k,t] + λpe[(k[2],k[1],k[3]),t] +
        (1 - Bparams[k,t])*λlt2[k,t] + 2*Rdict[k]*λpf[k,t]*Bparams[k,t] - μ1[k,t,1]*Bparams[k,t] == 0);
    # q
    @constraint(dp, qConstraint[k in fData.brList, t in td:T], -λqb[k[1],t] + λqe[k,t] + λqe[(k[2],k[1],k[3]),t] +
        (1 - Bparams[k,t])*λlt3[k,t] + 2*Xdict[k]*λpf[k,t]*Bparams[k,t] - μ1[k,t,2]*Bparams[k,t] == 0);
    # v
    @constraint(dp, vConstraint[i in fData.IDList, t in td:T], λvu[i,t] + λvl[i,t] + sum(λpf[k,t]*Bparams[k,t] for k in fData.branchDict2[i]) +
        sum(-λpf[k,t]*Bparams[k,t] for k in fData.branchDict1[i]) == 0);
    # lpp
    @constraint(dp, lppConstraint[i in fData.IDList, t in td:T], λpb[i,t] + lDict[i,t] >= 0);
    # lpm
    @constraint(dp, lpmConstraint[i in fData.IDList, t in td:T], -λpb[i,t] + lDict[i,t] >= 0);
    # lqp
    @constraint(dp, lqpConstraint[i in fData.IDList, t in td:T], λqb[i,t] + lDict[i,t] >= 0);
    # lqm
    @constraint(dp, lqmConstraint[i in fData.IDList, t in td:T], -λqb[i,t] + lDict[i,t] >= 0);

    # w
    @constraint(dp, wConstraint[i in bData.IDList, t in td:(T - 1)], λwu[i,t] + λbInv[i,t] - λbInv[i,t+1] -
        sum(λcuts[item]*cutDict[(item[1],item[2])][item[3]].γ[i] for item in wCuts[t] if t >= td + τ) >= 0);
    @constraint(dp, wConstraintT[i in bData.IDList], λwu[i,T] + λbInv[i,T] >= 0);# time after td
    @constraint(dp, wConstraint1[i in bData.IDList], -λbInv[i,td] + λw[i] >= 0);   # td-1
    # u
    @constraint(dp, uConstraint[i in bData.IDList], λuu[i] + λu[i] - sum(ν2[i,t] for t in td:T) -
        sum(λcuts[item]*cutDict[(item[1],item[2])][item[3]].μ[i] for item in keys(λcuts) if item[1] >= td + τ + 1) >= 0);
    # zp
    @constraint(dp, zpConstraint1[i in bData.IDList, t in td:T], λpb[bData.Loc[i],t] + sum(λbCoeff[i,l,t] for l in 1:length(bData.ηα[i])) -
        μ2[i,t,1] == 0);
    # zq
    @constraint(dp, zqConstraint1[i in bData.IDList, t in td:T], λqb[bData.Loc[i],t] - μ2[i,t,2] == 0);
    # y
    @constraint(dp, yConstraint1[i in bData.IDList, t in td:T], λbInv[i,t]*Δt - sum(bData.ηα[i][l]*λbCoeff[i,l,t] for l in 1:length(bData.ηα[i])) == 0);

    # θ
    @constraint(dp, θConstraint1[tp in (td + τ + 1):T, ω in Ω; (tp,ω) in keys(cutDict)],
        sum(λcuts[tp,ω,l] for l in 1:length(cutDict[tp,ω])) + θDict[tp,ω] >= 0);

    # f
    @constraint(dp, fConstraint3[i in fData.genIDList, t in td:T;fData.cp[i].n == 3],
        -λgcAux1[i,t] - λgcAux2[i,t] + fDict[i,t] >= 0);
    @constraint(dp, fConstraint2[i in fData.genIDList, t in td:T;fData.cp[i].n == 2], λgc2[i,t] + fDict[i,t] >= 0);

    # tAux1
    @constraint(dp, tAux1Constraint[i in fData.genIDList, t in td:T;fData.cp[i].n == 3],
        -ν3[i,t] + λgcAux1[i,t] == 0);
    # tAux2
    @constraint(dp, tAux2Constraint[i in fData.genIDList, t in td:T;fData.cp[i].n == 3],
        -μ3[i,t,1] + λgcAux2[i,t] == 0);
    # tAux3
    @constraint(dp, tAux3Constraint[i in fData.genIDList, t in td:T;fData.cp[i].n == 3],
        -μ3[i,t,2] + λgcAux3[i,t] == 0);

    # SOCP
    if qpopt
        for t in td:T
            for k in fData.brList
                if Bparams[k,t] == 1
                    @constraint(dp, μ1[k,t,1]^2 + μ1[k,t,2]^2 <= ν1[k,t]^2);
                end
            end
            for i in bData.IDList
                @constraint(dp, μ2[i,t,1]^2 + μ2[i,t,2]^2 <= ν2[i,t]^2);
            end
            for i in fData.genIDList
                if fData.cp[i].n == 3
                    @constraint(dp, μ3[i,t,1]^2 + μ3[i,t,2]^2 <= ν3[i,t]^2);
                end
            end
        end
    else
        for t in td:T
            for k in fData.brList
                if Bparams[k,t] == 1
                    μ1List = [μ1[k,t,1],μ1[k,t,2]];
                    @constraint(dp,norm(μ1List) <= ν1[k,t]);
                end
            end
            for i in bData.IDList
                μ2List = [μ2[i,t,1],μ2[i,t,2]];
                @constraint(dp,norm(μ2List) <= ν2[i,t]);
            end
            for i in fData.genIDList
                if fData.cp[i].n == 3
                    μ3List = [μ3[i,t,1],μ3[i,t,2]];
                    @constraint(dp,norm(μ3List) <= ν3[i,t]);
                end
            end
        end
    end

    # objective function
    @objective(dp, Max, sum(sum(-λspu[i,t]*fData.Pmax[i]*Bparams[i,t] - λspl[i,t]*fData.Pmin[i]*Bparams[i,t] -
            λsqu[i,t]*fData.Qmax[i]*Bparams[i,t] - λsql[i,t]*fData.Qmin[i]*Bparams[i,t] -
            λru[i,t]*fData.RU[i]*Bparams[i,t] - λrd[i,t]*fData.RD[i]*Bparams[i,t] for i in fData.genIDList) for t in td:T) +
            sum(sum(-λvu[i,t]*fData.Vmax[i]^2 - λvl[i,t]*fData.Vmin[i]^2 - λpb[i,t]*dData.pd[i][t] - λqb[i,t]*dData.qd[i][t]
            for i in fData.IDList) for t in td:T) + sum(sum(-ν1[k,t]*fData.rateA[k]*Bparams[k,t] for k in fData.brList) +
            sum(-λwu[i,t]*bData.cap[i] - sum(λbCoeff[i,l,t]*bData.ηβ[i][l] for l in 1:length(bData.ηα[i]))
            for i in bData.IDList) for t in td:T) + sum(-λsp[i]*currentSol.sp[i,td - 1] for i in fData.genIDList) +
            sum(-λuu[i]*bData.uCap[i] - λu[i]*currentSol.u[i] - λw[i]*currentSol.w[i,td - 1] for i in bData.IDList) +
            sum(sum(-λgcAux1[i,t]*((fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1])) -
            λgcAux2[i,t]*((-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1])) -
            λgcAux3[i,t]*(fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1]))) for i in fData.genIDList if fData.cp[i].n == 3) for t in td:T) -
            sum(sum(sum(λcuts[tp,ω,l]*cutDict[tp,ω][l].vhat for l in 1:length(cutDict[tp,ω])) for ω in Ω if (tp,ω) in keys(cutDict)) for tp in (td + τ + 1):T)
            );

    # return the cut or the dual problem
    if solveOpt
        statusDp = solve(dp);
        println(statusDp, " ", td, " ", ωd);
        # obtain the primal solutions & obj value
        vhat = getobjectivevalue(dp);

        solSp = Dict();
        solw = Dict();
        solu = Dict();
        for i in fData.genIDList
            solSp[i] = currentSol.sp[i,td - 1];
        end
        for i in bData.IDList
            solu[i] = currentSol.u[i];
            solw[i] = currentSol.w[i,td - 1];
        end
        # obtain the dual solutions
        dsolλ = Dict();
        dsolγ = Dict();
        dsolμ = Dict();
        for i in fData.genIDList
            dsolλ[i] = -getvalue(λsp[i]);
            vhat -= dsolλ[i]*solSp[i];
            if abs(dsolλ[i]) < 1e-4
                if dsolλ[i] < 0
                    vhat += dsolλ[i]*fData.Pmax[i];
                end
                dsolλ[i] = 0;
            end
        end
        for i in bData.IDList
            dsolγ[i] = -getvalue(λw[i]);
            vhat -= dsolγ[i]*solw[i];
            if abs(dsolγ[i]) < 1e-4
                if dsolγ[i] < 0
                    vhat += dsolγ[i]*bData.cap[i];
                end
                dsolγ[i] = 0;
            end
            dsolμ[i] = -getvalue(λu[i]);
            vhat -= dsolμ[i]*solu[i];
            if abs(dsolμ[i]) < 1e-4
                if dsolμ[i] < 0
                    vhat += dsolμ[i]*bData.uCap[i];
                end
                dsolμ[i] = 0;
            end
        end

        cutTemp = cutData(statusDp,dsolλ,dsolγ,dsolμ,vhat);
        return cutTemp;
    else
        return dp;
    end

    # SOCP constraints
    # @constraint(mp, lineThermal1[k in fData.brList, t in td:T;Bparams[k,t] == 1], p[k,t]^2 + q[k,t]^2 <= fData.rateA[k]^2);
    # @constraint(mp, bThermal[i in bData.IDList, t in td:T], zp[i,t]^2 + zq[i,t]^2 <= u[i]^2);
    # @constraint(mp,gcAux1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux1[i,t] == f[i,t] +
    #     (fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    # @constraint(mp,gcAux2[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux2[i,t] == f[i,t] +
    #     (-fData.cp[i].params[1] + fData.cp[i].params[2]^2)/(4*fData.cp[i].params[1]));
    # @constraint(mp,gcAux3[i in fData.genIDList, t in td:T; fData.cp[i].n == 3], tAux3[i,t] == sqrt(fData.cp[i].params[1])*sp[i,t] +
    #     fData.cp[i].params[2]/(2*sqrt(fData.cp[i].params[1])));
    # @constraint(mp, genCost1[i in fData.genIDList, t in td:T; fData.cp[i].n == 3],
    #     norm([tAux2[i,t],tAux3[i,t]]) <= tAux1[i,t]);
end

function constructBackwardM(td, T, Δt, trialPaths, matchedTrial, qpopt = false, hardened = [])
    # construct the math program given the state variables and current stage
    Ω = [ω for ω in keys(pDistr.ωDistrn)];
    paraSet = Iterators.product(Ω,matchedTrial);

    cutCurrentData = pmap(item -> dfBuild_D(td, pDistr.ωDict[item[1]], trialPaths[item[2]], pDistr.ωτ[item[1]], Δt, T, qpopt, true, hardened), paraSet);
    for j in procs()
        remotecall_fetch(cutUpdate,j,td,Ω,paraSet,cutCurrentData);
    end

    # for ω in Ω
    #     # solve the later stage problem
    #     cutCurrent = fBuild_D(td, ω, prevSol, τ, Δt, T, fData, bData, dData, pDistr, cutDict);
    #     if (td,ω) in keys(cutDict)
    #         push!(cutDict[td,ω],cutCurrent);
    #     else
    #         cutDict[td,ω] = [cutCurrent];
    #     end
    # end
    # return cutDict;
end

function exeBackward(T, Δt, trialPaths, qpopt = false, hardened = [])
    # execution of forward pass
    # input: trialPaths: the collection of trial points
    #        cutDict: previously generated cuts (preset in every core)
    # output: update the cutDict
    tpDict = Dict();
    for n in keys(trialPaths)
        tpDict[n] = [trialPaths[n][i][2] for i in 1:length(trialPaths[n])];
    end
    for t in T:-1:2
        matchedTrial = [];
        possiblePath = [];
        for n in keys(trialPaths)
            pathTemp = [(trialPaths[n][i][2],trialPaths[n][i][3],trialPaths[n][i][4]) for i in 2:length(trialPaths[n]) if trialPaths[n][i][2] < t];
            if t in tpDict[n]
                if !(pathTemp in possiblePath)
                    push!(matchedTrial,n);
                    push!(possiblePath,pathTemp);
                end
            end
        end
        if possiblePath != []
            constructBackwardM(t, T, Δt, trialPaths, matchedTrial, qpopt, hardened);
        end
        println("Time $(t) Passed");
    end
end

function exeBackward_last(T, Δt, trialPaths, qpopt = false, hardened = [])
    # execution of forward pass
    # input: trialPaths: the collection of trial points
    #        cutDict: previously generated cuts (preset in every core)
    # output: update the cutDict
    tpDict = Dict();
    for n in keys(trialPaths)
        tpDict[n] = [trialPaths[n][i][2] for i in 1:length(trialPaths[n])];
    end
    τMin = minimum([i for i in values(pDistr.ωτ)]);
    for t in T:-1:(T-τMin)
        matchedTrial = [];
        possiblePath = [];
        for n in keys(trialPaths)
            pathTemp = [(trialPaths[n][i][2],trialPaths[n][i][3],trialPaths[n][i][4]) for i in 2:length(trialPaths[n]) if trialPaths[n][i][2] < t];
            if t in tpDict[n]
                if !(pathTemp in possiblePath)
                    push!(matchedTrial,n);
                    push!(possiblePath,pathTemp);
                end
            end
        end
        if possiblePath != []
            constructBackwardM(t, T, Δt, trialPaths, matchedTrial, qpopt, hardened);
        end
        println("Time $(t) Passed");
    end
end

function exeBackwardAll(T, Δt, trialPaths, qpopt = false, hardened = [])
    # execution of forward pass
    # input: trialPaths: the collection of trial points
    #        cutDict: previously generated cuts (preset in every core)
    # output: update the cutDict
    tpDict = Dict();
    for n in keys(trialPaths)
        possibleTList = [t for t in 1:T];
        for i in 2:length(trialPaths[n])
            t = trialPaths[n][i][2];
            for tp in 1:trialPaths[n][i][4]
                delInd = findfirst(x -> x == t+tp, possibleTList);
                if delInd != nothing
                    deleteat!(possibleTList,delInd);
                end
            end
        end
        tpDict[n] = possibleTList;
    end
    for t in T:-1:2
        possiblePath = [];
        matchedTrial = [];
        for n in keys(trialPaths)
            if t in tpDict[n]
                pathTemp = [(trialPaths[n][i][2],trialPaths[n][i][3],trialPaths[n][i][4]) for i in 2:length(trialPaths[n]) if trialPaths[n][i][2] < t];
                if !(pathTemp in possiblePath)
                    push!(matchedTrial,n);
                    push!(possiblePath,pathTemp);
                end
            end
        end
        if possiblePath != []
            constructBackwardM(t, T, Δt, trialPaths, matchedTrial, qpopt, hardened);
        end
        println("Time $(t) Passed");
    end
end
