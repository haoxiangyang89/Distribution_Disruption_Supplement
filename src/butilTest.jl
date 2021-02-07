# test the utilization rate of batteries
using Distributed;
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;
dataList = Dict();

caseList = [13,33,123];
Δt = 0.25;
iterMax = 20;
T = 24;
NN = 1000;
τ = 4;
pathListDRaw = load("pathHist_1000.jld");

for ci in 1:length(caseList)
    pathDictA = pathListDRaw["pathDict"][ci];

    dataList[ci] = Dict();
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T,τ);
    end

    # select a preset pathDict
    pathDict = pathDictA[T];
    solDet,costDet = exeDet(T, Δt, fData, bData, dData, pDistr, NN, pathDict);

    # train the stochastic strategy
    cutDictPG = preGen(T, Δt, N, iterMax);
    # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, 20, 20, cutDictPG);
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
        max(Int64(round(500/N)),20), max(Int64(round(500/N)),20),cutDictPG);
    for j in procs()
        remotecall_fetch(cutIni,j,cutDict);
    end

    solSDDP, LBSDDP, costSDDP = exeForward(T, Δt, NN, false, pathDict);

    # read the battery utilization data from solDet and solSDDP
    uDet = solDet[1][1][1].u;
    uSDDP = solSDDP[1][1][1].u;
    zDet = Dict();
    zSDDP = Dict();
    wDet = Dict();
    wSDDP = Dict();
    statD = Dict();
    statS = Dict();
    for i in 1:NN
        zDet[i] = Dict();
        zSDDP[i] = Dict();
        wDet[i] = Dict();
        wSDDP[i] = Dict();
        for b in bData.IDList
            zDet[i][b] = Dict();
            zSDDP[i][b] = Dict();
            wDet[i][b] = Dict();
            wSDDP[i][b] = Dict();
        end

        # for each scenario, obtain the deterministic solution and the stochastic solution
        dDetList = [solDet[i][j][2] for j in 1:length(solDet[i])];
        push!(dDetList,T+1);
        for j in 1:length(solDet[i])
            # obtain the inventory & actual injection to network until the next solution
            for t in dDetList[j]:(dDetList[j+1] - 1)
                for b in bData.IDList
                    zDet[i][b][t] = solDet[i][j][1].zp[b,t];
                    wDet[i][b][t] = solDet[i][j][1].w[b,t];
                end
            end
        end

        dSDDPList = [solSDDP[i][j][2] for j in 1:length(solSDDP[i])];
        push!(dSDDPList,T+1);
        for j in 1:length(solSDDP[i])
            # obtain the inventory & actual injection to network until the next solution
            for t in dSDDPList[j]:(dSDDPList[j+1] - 1)
                for b in bData.IDList
                    zSDDP[i][b][t] = solSDDP[i][j][1].zp[b,t];
                    wSDDP[i][b][t] = solSDDP[i][j][1].w[b,t];
                end
            end
        end
        # obtain the statistics to show the utilization
        for b in bData.IDList
            wListD = [wDet[i][b][t] for t in 1:T];
            meanwD = mean(wListD);
            wListS = [wSDDP[i][b][t] for t in 1:T];
            meanwS = mean(wListS);
            chargeD = 0;
            chargeS = 0;
            dischargeD = 0;
            dischargeS = 0;
            cdCoeff = [];
            wDet[i][b][0] = bData.bInv[b];
            wSDDP[i][b][0] = bData.bInv[b];
            for t in 1:T
                yD = wDet[i][b][t] - wDet[i][b][t - 1];
                zD = zDet[i][b][t];
                if zD > 0
                    dischargeD += zD;
                else
                    chargeD += zD;
                end

                yS = (wSDDP[i][b][t - 1] - wSDDP[i][b][t])/Δt;
                zS = zSDDP[i][b][t];
                if zS > 0
                    dischargeS += zS;
                    zSs = yzCal(yS,bData,b);
                    if (zS/zSs <= 1)&(zS/zSs > 0)
                        push!(cdCoeff,zS/zSs);
                    end
                elseif zS < 0
                    chargeS += zS;
                    zSs = yzCal(yS,bData,b);
                    if (zSs/zS <= 1)&(zSs/zS > 0)
                        push!(cdCoeff,zSs/zS);
                    end
                end
            end
            statD[i,b] = [meanwD,chargeD,dischargeD];
            if cdCoeff != []
                cdCoeffAvg = mean(cdCoeff);
            else
                cdCoeffAvg = -Inf;
            end
            statS[i,b] = [meanwS,chargeS,dischargeS,cdCoeffAvg,cdCoeff];
        end
    end
    save("butilResults_$(ci).jld","detOut",[uDet,statD],"stochOut",[uSDDP,statS]);
end
