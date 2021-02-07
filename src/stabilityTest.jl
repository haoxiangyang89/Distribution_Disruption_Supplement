# test of solution stability
using Distributed;
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;
dataList = Dict();

caseList = [13,33,123];
Δt = 0.25;
iterMax = 20;
trainNo = 20;
NN = 5000;

for ci in 1:length(caseList)
    TList = [24,36,48,72,96];
    pathListDRaw = load("pathHist_5000.jld");
    pathDictA = pathListDRaw["pathDict"][ci];
    dataList[ci] = Dict();

    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T,τ);
        end
        pathDict = pathDictA[T];
        dataList[ci][T] = Dict();

        # train trainNo of strategies and obtain the 95% PI
        for tn in 1:trainNo
            # train the stochastic programming strategy
            cutDictPG = preGen(T, Δt, N, iterMax);
            #cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, 30, 30, cutDictPG);
            cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
                max(Int64(round(500/N)),20), max(Int64(round(500/N)),20),cutDictPG);
            for j in procs()
                remotecall_fetch(cutIni,j,cutDict);
            end

            # test against 5000 samples
            solSDDP, LBSDDP, costSDDP = exeForward(T, Δt, NN, false, pathDict);
            listSDDP = [costSDDP[i] for i in 1:NN];
            meanSDDP = mean(listSDDP);
            sigmaSDDP = std(listSDDP);
            println(round(meanSDDP,digits = 2)," ",round(meanSDDP - 1.96*sigmaSDDP,digits = 2)," ",round(meanSDDP + 1.96*sigmaSDDP,digits = 2));
            dataList[ci][T][tn] = [LBSDDP, costSDDP, listSDDP, meanSDDP, sigmaSDDP];
            save("stabilityResults.jld","data",dataList);
        end
    end
end
