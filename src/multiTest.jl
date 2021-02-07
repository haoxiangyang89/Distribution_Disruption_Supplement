# test multiple components broken
# upper bound tests
using Distributed;
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
#pmap(i -> importIpopt(),1:30);

dataList = Dict();
iterMax = 20;

caseList = [13,33,123];
T = 24;
Δt = 0.25;
N = 5;
ci = 3;
NN = 1000;

# obtain an N-1 policy
τ = Int64(1/6*T);
for j in procs()
    remotecall_fetch(readInData,j,ci,caseList,T,τ);
end

startT = time();
# cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false,
#     Int64(round(100/N)),Int64(round(100/N)), cutDictPG, false, 200, [], 0, pathDict);
cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(T, Δt, N, false, false,
    max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), Dict(), false, 500, [], 0, Dict());
elapsedT = time() - startT;

cutDictN1 = deepcopy(cutDict);

for j in procs()
    remotecall_fetch(readInData_tau,j,ci,caseList,T,1e4,0,"Multi",τ);
end

# get multiple disruption cuts
startT = time();
# cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, true, false,
#     Int64(round(100/N)),Int64(round(100/N)), cutDictPG, false, 200, [], 0, pathDict);
cutDict,LBHistm,UBHistm,UBuHistm,UBlHistm,timeHistm = solveMain(T, Δt, N, false, false,
    max(Int64(round(500/N)),20), max(Int64(round(500/N)),20), Dict(), false, 500, [], 0, Dict());
elapsedT = time() - startT;
dataList = [LBHist,UBHist,UBuHist,UBlHist,timeHist,
            LBHistm,UBHistm,UBuHistm,UBlHistm,timeHistm,elapsedT];
cutDictNM = deepcopy(cutDict);

# plug in the N-1 policy in the multi-disruption case
# sample 1000 scenarios
pathDict = pmap(i -> simuPath(T,pDistr), 1:1000);
solDet,costDet = exeDet(T, Δt, fData, bData, dData, pDistr, NN, pathDict);
for j in procs()
    remotecall_fetch(cutIni,j,cutDictN1);
end
solSDDP1, LBSDDP1, costSDDP1 = exeForward(T, Δt, NN, false, pathDict);
for j in procs()
    remotecall_fetch(cutIni,j,cutDictNM);
end
solSDDPM, LBSDDPM, costSDDPM = exeForward(T, Δt, NN, false, pathDict);

save("MultiResults.jld","MultiOut",dataList,
    #"N1Cuts",cutDictN1,"MultiCuts",cutDict,
    "Results",[solDet, costDet, solSDDP1, LBSDDP1, costSDDP1, solSDDPM, LBSDDPM, costSDDPM]);
