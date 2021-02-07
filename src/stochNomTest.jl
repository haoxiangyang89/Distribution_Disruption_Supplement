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
    stochNomOut = Dict();
    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T);
            #remotecall_fetch(readInData_old,j,T,ωSet0,10000,0);
        end
        pathDict = pathDictA[T];
        pathDict1 = Dict();
        noP = 1;
        foundBool = false;
        while !(foundBool)
            if length(pathDict[noP]) == 1
                pathDict1[1] = pathDict[noP];
                foundBool = true;
            else
                noP += 1;
            end
        end

        # train the stochastic programming strategy
        cutDictPG = preGen(τ, T, Δt, N, iterMax);
        # cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false, 10, 10, cutDictPG);
        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false,
            max(Int64(round(500/N)),20), max(Int64(round(500/N)),20),cutDictPG);
        for j in procs()
            remotecall_fetch(cutIni,j,cutDict);
        end

        solSDDP, LBSDDP, costSDDP = exeForward(τ, T, Δt, 1, false, pathDict1);
        stochNomOut[T] = [solSDDP,LBSDDP,costSDDP];

        save("stochNomResults_$(ci).jld", "stochNomOut", stochNomOut);
    end
end
