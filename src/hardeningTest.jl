# test 3 hardening tests
using Distributed;
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

@everywhere function obtainStat(fData,NN,T,solSDDP)
    # obtain the load mismatch and generation
    spSDDP = Dict();
    sqSDDP = Dict();
    LSDDP = Dict();
    spList = [];
    sqList = [];
    LList = [];
    for i in 1:NN
        spSDDP[i] = Dict();
        sqSDDP[i] = Dict();
        LSDDP[i] = Dict();
        for g in fData.genIDList
            spSDDP[i][g] = Dict();
            sqSDDP[i][g] = Dict();
        end
        for ib in fData.IDList
            LSDDP[i][ib] = Dict();
        end

        dSDDPList = [solSDDP[i][j][2] for j in 1:length(solSDDP[i])];
        push!(dSDDPList,T+1);
        for j in 1:length(solSDDP[i])
            # obtain the inventory & actual injection to network until the next solution
            for t in dSDDPList[j]:(dSDDPList[j+1] - 1)
                for g in fData.genIDList
                    spSDDP[i][g][t] = solSDDP[i][j][1].sp[g,t];
                    sqSDDP[i][g][t] = solSDDP[i][j][1].sq[g,t];
                end
                for ib in fData.IDList
                    LSDDP[i][ib][t] = solSDDP[i][j][1].lp[ib,t] + solSDDP[i][j][1].lq[ib,t];
                end
            end
        end
        # obtain the statistics to show total generation and power mismatch
        push!(spList,sum(sum(spSDDP[i][g][t] for t in 1:T) for g in fData.genIDList));
        push!(sqList,sum(sum(sqSDDP[i][g][t] for t in 1:T) for g in fData.genIDList));
        push!(LList,sum(sum(LSDDP[i][ib][t] for t in 1:T) for ib in fData.IDList));
    end
    return spList,sqList,LList;
end

N = 5;

caseList = [13,33,123];
Δt = 0.25;
iterMax = 20;
T = 24;
NN = 1000;
pathTrain = load("pathHist_600.jld");
pathListDRaw = load("pathHist_1000.jld");

for ci in 1:length(caseList)
    pathDict = pathTrain["pathDict"][ci][T];
    pathDictTest = pathListDRaw["pathDict"][ci][T];

    dataList = Dict();
    τ = Int64(1/6*T);
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T,τ,1e4);
        #remotecall_fetch(readInData_old,j,T,[(2,9),(8,12),(10,13)],10000,0);
    end
    cutDictPG = preGen(T, Δt, N, iterMax, false, Dict(), []);
    # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false,
    #     20, 20, cutDictPG, false, 0, [], 0, pathDict);
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
        max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), cutDictPG, false, 0, []);
    solSDDP, LBSDDP, costSDDP = exeForward(T, Δt, NN, false, pathDictTest);
    spList,sqList,LList = obtainStat(fData,NN,T,solSDDP);

    dataList["NoD"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,LBSDDP,costSDDP,spList,sqList,LList];

    for ω in keys(pDistr.ωDistrn)
    #for ω in [(2,9),(8,12),(10,13)]
        cutDictPG = preGen(T, Δt, N, iterMax, false, Dict(), pDistr.ωDict[ω]);
        # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false,
        #     20, 20, cutDictPG, false, 0, [ω], 0, pathDict);
        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
            max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), cutDictPG, false, 0, pDistr.ωDict[ω]);
        solSDDP, LBSDDP, costSDDP = exeForward(T, Δt, NN, false, pathDictTest, pDistr.ωDict[ω]);
        spList,sqList,LList = obtainStat(fData,NN,T,solSDDP);

        dataList[ω] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,LBSDDP,costSDDP,spList,sqList,LList];
        save("hardResults_$(ci).jld","data",dataList);
    end
end
