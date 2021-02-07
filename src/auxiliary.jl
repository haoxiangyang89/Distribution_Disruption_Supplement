# auxiliary functions

using Distributions;
import Base.rand;

struct CategoricalSamplerNew <: Sampleable{Univariate,Discrete}
    mass::Vector{Float64}
    category::Vector{Float64}
end

function rand(s::CategoricalSamplerNew)
    s.category[rand(Categorical(s.mass))]
end

function genScenario_old(pDistr)
    # generate a disruption time
    tSupport = [i for i in keys(pDistr.tDistrn)];
    tProb = [pDistr.tDistrn[i] for i in tSupport];
    tDistrObj = Categorical(tProb);
    t = rand(tDistrObj);

    # generate a disruption location
    ωSupport = [i for i in keys(pDistr.ωDistrn)];
    ωProb = [pDistr.ωDistrn[i] for i in ωSupport];
    ωDistrObj = Categorical(ωProb);
    ω = rand(ωDistrObj);
    return tSupport[t],ωSupport[ω];
end

function genScenario(pDistr)
    # generate a disruption time
    tSupport = [i for i in keys(pDistr.tDistrn)];
    tProb = [pDistr.tDistrn[i] for i in tSupport];
    tDistrObj = Categorical(tProb);
    t = rand(tDistrObj);

    # generate a disruption location
    ωSupport = [i for i in keys(pDistr.ωDistrn)];
    ωProb = [pDistr.ωDistrn[i] for i in ωSupport];
    ωDistrObj = Categorical(ωProb);
    ω = rand(ωDistrObj);
    scenInd = ωSupport[ω];
    return tSupport[t],pDistr.ωDict[scenInd],pDistr.ωτ[scenInd];
end

function modifyOmega(pDistr,hardComp)
    ωDistrNew = Dict();
    ωDictNew = Dict();
    ωτNew = Dict();
    for iHard in hardComp
        releaseProb += pDistr.ωDistrn[iHard];
    end
    avgNo = length(values(pDistr.ωDistrn)) - length(hardComp);
    for i in keys(pDistr.ωDistrn)
        if !(i in hardComp)
            # if it is not the hardened component
            ωDistrNew[i] = pDistr.ωDistrn[i] + releaseProb/avgNo;
            ωDictNew[i] = pDistr.ωDict[i];
            ωτNew[i] = pDistr.ωτ[i];
        end
    end
    pDistrNew = probDistrn(pDistr.tDistrn,ωDistrNew,ωDictNew);
    return pDistrNew;
end

function modifyT(pDistr,λD,T)
    # modify the time distribution using the given λD
    tDistrnNew = Dict();
    for t in 1:(T-1)
        tDistrnNew[t] = exp(-λD*(t - 1)) - exp(-λD*t);
    end
    # probability of the time T+1
    tDistrnNew[T] = 1 - sum([tDistrnNew[t] for t in 1:(T-1)]);
    pDistrNew = probDistrn(tDistrnNew,pDistr.ωDistrn,pDistr.ωDict,pDistr.ωτ);
    return pDistrNew;
end

function calObj(mp, td, τ, T, fData, bData, pDistr, sp, lpp, lqp, lpm, lqm, θ)
    # function to formulate the objective function for f_D
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    objExpr = @expression(mp, 0);
    for tp in 1:maximum(keys(pDistr.tDistrn))
        dExpr = @expression(mp, 0);
        if tp <= T - (td + τ)
            for t in td:(tp + td + τ - 1)
                for i in fData.genIDList
                    # add generator cost
                    if fData.cp[i].n == 3
                        dExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
                    elseif fData.cp[i].n == 2
                        dExpr += fData.cp[i].params[1]*sp[i,t];
                    end
                end
                # add load shed cost
                dExpr += fData.cz*(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*(dExpr + sum(pDistr.ωDistrn[ω]*θ[tp + td + τ,ω] for ω in Ω));
        else
            for t in td:T
                for i in fData.genIDList
                    # add generator cost
                    if fData.cp[i].n == 3
                        dExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
                    elseif fData.cp[i].n == 2
                        dExpr += fData.cp[i].params[1]*sp[i,t];
                    end
                end
                # add load shed cost
                dExpr += fData.cz*(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*dExpr;
        end
    end
    return objExpr;
end

function calObjDet(mp, T, fData, bData, sp, lpp, lqp, lpm, lqm, u)
    # function to formulate the objective function for deterministic case

    objExpr = @expression(mp, sum(bData.cost[i]*u[i] for i in bData.IDList) +
        fData.cz*sum(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList) for t in 1:T));
    for t in 1:T
        for i in fData.genIDList
            # add generator cost
            if fData.cp[i].n == 3
                objExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
            elseif fData.cp[i].n == 2
                objExpr += fData.cp[i].params[1]*sp[i,t];
            end
        end
    end
    return objExpr;
end

function calObj1(mp, T, fData, bData, pDistr, sp, lpp, lqp, lpm, lqm, θ, u)
    # function to formulate the objective function for f_1
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    objExpr = @expression(mp, sum(bData.cost[i]*u[i] for i in bData.IDList));
    for tp in 1:maximum(keys(pDistr.tDistrn))
        dExpr = @expression(mp, 0);
        if tp < T
            for t in 1:tp
                for i in fData.genIDList
                    # add generator cost
                    if fData.cp[i].n == 3
                        dExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
                    elseif fData.cp[i].n == 2
                        dExpr += fData.cp[i].params[1]*sp[i,t];
                    end
                end
                # add load shed cost
                dExpr += fData.cz*(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*(dExpr + sum(pDistr.ωDistrn[ω]*θ[tp + 1,ω] for ω in Ω));
        else
            for t in 1:T
                for i in fData.genIDList
                    # add generator cost
                    if fData.cp[i].n == 3
                        dExpr += fData.cp[i].params[1]*(sp[i,t]^2) + fData.cp[i].params[2]*sp[i,t];
                    elseif fData.cp[i].n == 2
                        dExpr += fData.cp[i].params[1]*sp[i,t];
                    end
                end
                # add load shed cost
                dExpr += fData.cz*(sum(lpp[i,t] + lqp[i,t] + lpm[i,t] + lqm[i,t] for i in fData.IDList));
            end
            objExpr += pDistr.tDistrn[tp]*dExpr;
        end
    end
    return objExpr;
end

function calDualC1(T, fData, pDistr)
    # function to formulate the objective function for f_D
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    fDict = Dict();
    lDict = Dict();
    θDict = Dict();
    for tp in 1:maximum(keys(pDistr.tDistrn))
        if tp < T
            for t in 1:tp
                for i in fData.genIDList
                    if !((i,t) in keys(fDict))
                        fDict[i,t] = pDistr.tDistrn[tp];
                    else
                        fDict[i,t] += pDistr.tDistrn[tp];
                    end
                end
                for i in fData.IDList
                    if !((i,t) in keys(lDict))
                        lDict[i,t] = pDistr.tDistrn[tp]*fData.cz;
                    else
                        lDict[i,t] += pDistr.tDistrn[tp]*fData.cz;
                    end
                end
            end
            for ω in Ω
                θDict[tp + 1,ω] = pDistr.tDistrn[tp]*pDistr.ωDistrn[ω];
            end
        else
            for t in 1:T
                for i in fData.genIDList
                    if !((i,t) in keys(fDict))
                        fDict[i,t] = pDistr.tDistrn[tp];
                    else
                        fDict[i,t] += pDistr.tDistrn[tp];
                    end
                end
                for i in fData.IDList
                    if !((i,t) in keys(lDict))
                        lDict[i,t] = pDistr.tDistrn[tp]*fData.cz;
                    else
                        lDict[i,t] += pDistr.tDistrn[tp]*fData.cz;
                    end
                end
            end
        end
    end
    return fDict,lDict,θDict;
end

function calCostF(costn, currentSol, T, fData, nowT, disT)
    # function to calculate costs for forward pass
    costn += sum(sum(fData.cz*(abs(currentSol.lp[i,t]) + abs(currentSol.lq[i,t])) for i in fData.IDList) for t in nowT:(disT - 1));
    for t in nowT:(disT - 1)
        for i in fData.genIDList
            # add generator cost
            if fData.cp[i].n == 3
                costn += fData.cp[i].params[1]*(currentSol.sp[i,t]^2) + fData.cp[i].params[2]*currentSol.sp[i,t];
            elseif fData.cp[i].n == 2
                costn += fData.cp[i].params[1]*currentSol.sp[i,t];
            end
        end
    end
    return costn;
end

function calDualC(td, τ, T, fData, pDistr)
    # function to formulate the objective function for f_D
    Ω = [ω for ω in keys(pDistr.ωDistrn)];

    fDict = Dict();
    lDict = Dict();
    θDict = Dict();
    for tp in 1:maximum(keys(pDistr.tDistrn))
        if tp <= T - (td + τ)
            for t in td:(tp + td + τ - 1)
                for i in fData.genIDList
                    if !((i,t) in keys(fDict))
                        fDict[i,t] = pDistr.tDistrn[tp];
                    else
                        fDict[i,t] += pDistr.tDistrn[tp];
                    end
                end
                for i in fData.IDList
                    if !((i,t) in keys(lDict))
                        lDict[i,t] = pDistr.tDistrn[tp]*fData.cz;
                    else
                        lDict[i,t] += pDistr.tDistrn[tp]*fData.cz;
                    end
                end
            end
            for ω in Ω
                θDict[tp + td + τ,ω] = pDistr.tDistrn[tp]*pDistr.ωDistrn[ω];
            end
        else
            for t in td:T
                for i in fData.genIDList
                    if !((i,t) in keys(fDict))
                        fDict[i,t] = pDistr.tDistrn[tp];
                    else
                        fDict[i,t] += pDistr.tDistrn[tp];
                    end
                end
                for i in fData.IDList
                    if !((i,t) in keys(lDict))
                        lDict[i,t] = pDistr.tDistrn[tp]*fData.cz;
                    else
                        lDict[i,t] += pDistr.tDistrn[tp]*fData.cz;
                    end
                end
            end
        end
    end
    return fDict,lDict,θDict;
end

function simuPath(T,pDistr)
    pathList = [];
    nowT = 1;
    while nowT <= T
        tp,ωd,τω = genScenario(pDistr);
        push!(pathList, (tp,ωd,τω));
        if nowT == 1
            nowT += tp;
            nowT = min(nowT, T + 1);
        else
            nowT += tp + τω;
            nowT = min(nowT, T + 1);
        end
    end
    return pathList;
end

function pathSimu_cover(N,T,pDistr,genCutsCount,iterNo)
    # maximize the coverage from the sample
    pathSet = [];
    pathTimeList = [];
    for i in 1:(10*N)
        pathList = simuPath(T,pDistr);
        push!(pathSet,pathList);
        pathTimes = [];
        tNow = 1;
        for j in 1:length(pathList)
            if tNow == 1
                tNow += pathList[j][1];
            else
                tNow += pathList[j][3] + pathList[j][1];
            end
            if tNow < T
                push!(pathTimes, tNow);
            end
        end
        push!(pathTimeList, pathTimes);
    end

    sp = Model(solver = GurobiSolver(OutputFlag = 0));
    @variable(sp, x[i in 1:length(pathSet)], Bin);
    @constraint(sp, sum(x[i] for i in 1:length(pathSet)) == N);
    @objective(sp, Max, sum(sum((N*iterNo - genCutsCount[t])*x[i] for t in 2:T if t in [pathTimeList[i][j] for j in 1:length(pathTimeList[i])])
        for i in 1:length(pathSet)));
    solve(sp);
    pathSetsel = Dict();
    pathTimesel = Dict();
    iCount = 0;
    for i in 1:length(pathSet)
        if getvalue(sp[:x][i]) == 1
            iCount += 1;
            pathSetsel[iCount] = pathSet[i];
            pathTimesel[iCount] = pathTimeList[i];
        end
    end
    return pathSetsel,pathTimesel;
end

function pathSimu_cover_last(N,T,pDistr)
    # maximize the coverage for the terminal stage problem from the sample
    pathSet = [];
    pathTimeList = [];
    selNo = 0;
    while selNo < N
        pathList = simuPath(T,pDistr);
        pathTimes = [];
        tNow = 1;
        τList = false;
        for j in 1:length(pathList)
            if tNow == 1
                tNow += pathList[j][1];
            else
                tNow += pathList[j][3] + pathList[j][1];
            end
            if tNow < T
                push!(pathTimes, tNow);
                if tNow >= T - pathList[j][3]
                    τList = true;
                end
            end
        end
        if τList
            push!(pathSet,pathList);
            push!(pathTimeList, pathTimes);
            selNo += 1;
        end
    end
    return pathSet,pathTimeList;
end

# initialize the cutData in every core
function cutIni(cutData)
    global cutDict = cutData;
end

# update the cutData in every core with the new cuts
function cutUpdate(td,Ω,paraSet,cutCurrentData)
    # cutCurrentData is a list
    for ω in Ω
        itemInd = 0;
        for item in paraSet
            itemInd += 1;
            if item[1] == ω
                if (cutCurrentData[itemInd].solStatus == :Optimal)
                    if (td,ω) in keys(cutDict)
                        push!(cutDict[td,ω],cutCurrentData[itemInd]);
                    else
                        cutDict[td,ω] = [cutCurrentData[itemInd]];
                    end
                end
            end
        end
    end
end

function yzCal(y,bData,loc)
    noPiece = length(bData.ηα[loc]);
    fv = Inf;
    for n in 1:noPiece
        if y*bData.ηα[loc][n]+bData.ηβ[loc][n] < fv
            fv = y*bData.ηα[loc][n]+bData.ηβ[loc][n];
        end
    end
    return fv;
end

# change the cost of generation by a multiplier
function changeCost(fDatal, bDatal, cmulti)
    fDataLocal = deepcopy(fDatal);
    bDataLocal = deepcopy(bDatal);
    for i in fDataLocal.genIDList
        for nIter in 1:fDataLocal.cp[i].n
            fDataLocal.cp[i].params[nIter] = fDataLocal.cp[i].params[nIter]*cmulti^(fDataLocal.cp[i].n - nIter);
        end
    end
    fDataLocal.cz = fDataLocal.cz*cmulti;

    for i in bDataLocal.IDList
        bDataLocal.cost[i] = bDataLocal.cost[i]*cmulti;
    end

    global fData = fDataLocal;
    global bData = bDataLocal;
end

# break a component
function breakComponent(fDatal, bItem, bType)
    fDataLocal = deepcopy(fDatal);
    if bType == "g"
        # generators
        gListTemp = copy(fDataLocal.genIDList);
        deleteat!(gListTemp,findfirst(fDataLocal.genIDList,bItem));
        fDataLocal.genIDList = gListTemp;
    elseif bType == "l"
        # lines
        lineListTemp = deepcopy(fDataLocal.brList);
        for k in lineIDList
            if ((k[1],k[2]) == bItem) | ((k[2],k[1]) == bItem)
                deleteat!(lineListTemp,findfirst(fDataLocal.brList,k));
            end
        end
        fDataLocal.brList = lineListTemp;
    end
    global fData = fDataLocal;
end

# Add τ to the last element of pathDict
function addTau(pathDict,τ)
    newpathDict = Dict();
    for i in keys(pathDict)
        newPath = [];
        for item in pathDict[i]
            push!(newPath,(item[1],item[2],τ));
        end
        newpathDict[i] = newPath;
    end
    return newpathDict;
end

# get the index of scenarios
function reverseScen(pathDict,τ,pDistr)
    newpathDict = Dict();
    for i in keys(pathDict)
        newPath = [];
        for item in pathDict[i]
            # find the index of scenario
            jInd = 0;
            for j in keys(pDistr.ωDict)
                if item[2] in pDistr.ωDict[j]
                    jInd = j;
                end
            end
            push!(newPath,(item[1],pDistr.ωDict[jInd],pDistr.ωτ[jInd]));
        end
        newpathDict[i] = newPath;
    end
    return newpathDict;
end
