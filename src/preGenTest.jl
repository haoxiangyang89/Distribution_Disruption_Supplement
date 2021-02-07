# test to check if preGen helps reduce the solution time
using Distributed;
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;
dataList = Dict();

caseList = [13,33,123];
T = 96;
τ = Int64(1/6*T);
Δt = 0.25;
iterMax = 20;
pathTrain = load("pathHist_600.jld");

for ci in 1:length(caseList)
    pathDict = pathTrain["pathDict"][ci][T];
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T,τ);
    end

    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, 30, false,false, 2, 2);

    startPGT = time();
    cutDictPG = preGen(T, Δt, N, iterMax);
    preGenT = time() - startPGT;

    startT = time();
    # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, 10, 10, cutDictPG,
    #     false, 200, [], 0, pathDict);
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
        max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), cutDictPG, false, 200, [], 0, pathDict);
    elapsedT = time() - startT;
    dataList["preGen"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT,preGenT];

    save("pgResults_$(ci).jld","data",dataList);
end
