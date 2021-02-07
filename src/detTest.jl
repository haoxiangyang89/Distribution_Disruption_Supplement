# test deterministic vs. stochastic
using Distributed;
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

caseList = [13,33,123];
NN = 1000;
Δt = 0.25;
N = 5;
iterMax = 20;
pathListDRaw = load("pathHist_1000.jld");
TList = [24,36,48,72,96];

for ci in 1:length(caseList)
    pathDictA = pathListDRaw["pathDict"][ci];
    detOut = Dict();
    stochOut = Dict();
    detCostDict = Dict();
    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T,τ);
            #remotecall_fetch(readInData_old,j,T,ωSet0,10000,0);
        end

        # select a preset pathDict
        pathDict = pathDictA[T];
        solDet,costDet = exeDet(T, Δt, fData, bData, dData, pDistr, NN, pathDict);
        listDet = [costDet[i] for i in 1:NN];
        meanDet = mean(listDet);
        sigmaDet = std(listDet);
        println(round(meanDet,digits = 2)," ",
                round(meanDet - 1.96*sigmaDet,digits = 2)," ",
                round(meanDet + 1.96*sigmaDet,digits = 2));
        detOut[T] = [costDet,listDet,meanDet,sigmaDet];

        # train the stochastic programming strategy
        cutDictPG = preGen(T, Δt, N, iterMax);
        # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, 10, 10, cutDictPG);
        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
            max(Int64(round(500/N)),20), max(Int64(round(500/N)),20),cutDictPG);
        for j in procs()
            remotecall_fetch(cutIni,j,cutDict);
        end

        solSDDP, LBSDDP, costSDDP = exeForward(T, Δt, NN, false, pathDict);
        listSDDP = [costSDDP[i] for i in 1:NN];
        meanSDDP = mean(listSDDP);
        sigmaSDDP = std(listSDDP);
        println(round(meanSDDP,digits = 2)," ",
                round(meanSDDP - 1.96*sigmaSDDP,digits = 2)," ",
                round(meanSDDP + 1.96*sigmaSDDP,digits = 2));
        stochOut[T] = [LBSDDP, costSDDP, listSDDP, meanSDDP, sigmaSDDP, LBHist];

        detSol,detObj = detBuild(Δt, T, fData, bData, dData);
        detCostDict[T] = detObj;

        #println("Output: ",round(meanDet,1)," & ",round(meanSDDP,1)," & ", round(meanDet - meanSDDP,1)," & ", round((meanDet - meanSDDP)/meanDet*100,1));
        #save("detResults_old.jld", "detOut", detOut, "stochOut", stochOut, "nomOut", detCostDict);
        save("detResults_$(ci).jld", "detOut", detOut, "stochOut", stochOut, "nomOut", detCostDict);
    end
end

# for T in TList
#     println(T," & ",round(detOut[T][3],2)," & ",round(stochOut[T][4],2)," & ",round((detOut[T][3] - stochOut[T][4])/detOut[T][3]*100,2)," & ");
# end

# pathListData = pmap(i -> simuPath(τ,T,pDistr), 1:NN);
# pathDict = Dict();
# for i in 1:length(pathListData)
#     pathDict[i] = pathListData[i];
# end
# solDet,costDet = exeDet(τ, T, Δt, fData, bData, dData, pDistr, NN, pathDict);
# listDet = [costDet[i] for i in 1:NN];
# meanDet = mean(listDet);
# sigmaDet = std(listDet);
# println(round(meanDet,2)," ",round(meanDet - 1.96*sigmaDet,2)," ",round(meanDet + 1.96*sigmaDet,2));
#
# cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, 20, 20);
# solSDDP, LBSDDP, costSDDP = exeForward(τ, T, Δt, fData, bData, dData, pDistr, NN, cutDict, pathDict);
# listSDDP = [costSDDP[i] for i in 1:NN];
# meanSDDP = mean(listSDDP);
# sigmaSDDP = std(listSDDP);
# println(round(meanSDDP,2)," ",round(meanSDDP - 1.96*sigmaSDDP,2)," ",round(meanSDDP + 1.96*sigmaSDDP,2));
#
# outputData = [solDet,costDet,meanDet,sigmaDet,
#     cutDict,LBHist,UBHist,UBuHist,UBlHist,solSDDP,LBSDDP,costSDDP,meanSDDP,sigmaSDDP];
# save("detOut.jld","data",outputData);
