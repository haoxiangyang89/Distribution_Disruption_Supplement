# test to compare allGen and only disruption gen
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
pathTrain = load("pathHist_600.jld");

for ci in 1:length(caseList)
    pathDict = pathTrain["pathDict"][ci][T];
    for j in procs()
        remotecall_fetch(readInData,j,ci,caseList,T,τ);
    end

    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, 20, false,false, 2, 2);

    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, true, false, 20, 20, Dict(),
        false, 200, [], 0, pathDict);
    elapsedT = time() - startT;
    dataList["allGen"] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT];

    save("agResults_$(ci).jld","data",dataList);
end
