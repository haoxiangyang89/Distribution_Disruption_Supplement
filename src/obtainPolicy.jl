# obtain the policies for different time horizon
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

Δt = 0.25;
N = 5;
iterMax = 20;

caseList = [13,33,123];

policyDict = Dict();
for ci in 1:length(caseList)
    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T);
        end

        startPGT = time();
        cutDictPG = preGen(τ, T, Δt, N, iterMax);
        preGenT = time() - startPGT;

        startT = time();
        cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, max(Int64(round(200/N)),20), max(Int64(round(200/N)),20),cutDictPG);
        elapsedT = time() - startT;

        policyDict[T] = [cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist,elapsedT,preGenT];
        save("policy_trained_$(ci).jld","policyDict",policyDict);
    end
end
