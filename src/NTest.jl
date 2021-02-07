# test 3: number of trial paths tests
using Distributed;
addprocs(20);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
#pmap(i -> importIpopt(),1:30);

NList = [1,5,10,15,20];
dataList = Dict();
iterMax = 20;
NN = 5;

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

    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, 30, false,false, 2, 2);

    startPGT = time();
    cutDictPG = preGen(T, Δt, NN, iterMax);
    preGenT = time() - startPGT;
    cutDictPGOri = deepcopy(cutDictPG);

    for N in NList
        startT = time();
        # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false,
        #     Int64(round(100/N)),Int64(round(100/N)), cutDictPG, false, 200, [], 0, pathDict);
        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
            max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), cutDictPG, false, 200, [], 0, pathDict);
        elapsedT = time() - startT;
        dataList[N] = [LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT,preGenT];

        save("NResults_$(ci).jld","NOut",dataList);
        cutDictPG = deepcopy(cutDictPGOri);
    end
end
