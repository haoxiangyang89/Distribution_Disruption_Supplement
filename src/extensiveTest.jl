# test 1: compare SDDP and extensive formulation
using Distributed;
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();
#pmap(i -> importIpopt(),1:30);

τ = 2;
T = 12;
Δt = 0.25;
N = 5;
###########################################################################
data = [];
ωSet = [5,(2,9),(8,12)];

for T in [8,12,16]
    for j in procs()
        remotecall_fetch(readInData_old,j,T,ωSet,10000);
    end

    # global mExt = Model(solver = IpoptSolver(linear_solver = "ma27"));
    global mExt = Model(solver = GurobiSolver());
    startTE = time();
    global mExt = extForm(1, 0, [[],bData.bInv,[]], 1, τ, Δt, T, fData, bData, dData, pDistr, true);
    solve(mExt);
    elapsedTE = time() - startTE;
    mExtObj = getobjectivevalue(mExt);

    startT = time();
    #cutDict,LBHist,UBHist,UBuHist,UBlHist = solveMain(τ, T, Δt, fData, pDistr, bData, dData, N, false, 20, 20);
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, 100, 100);
    elapsedT = time() - startT;
    push!(data,(T,elapsedTE,mExtObj,elapsedT,LBHist));
end
save("extensiveSDDP.jld","data",data);

τ = 2;
T = 8;
Δt = 0.25;
N = 5;
disAdd = "testProbRead_96.csv";
pDistr = readDisruption(disAdd,"csv");
ωSet0 = [i for i in keys(pDistr.ωDistrn)];
ωSetTot = [[5,(2,9),(8,12)],
        [2,3,5,(2,3),(2,9),(8,12),(9,12),(10,13)],
        ωSet0];

for i in 1:3
    for j in procs()
        remotecall_fetch(readInData_old,j,T,ωSetTot[i]);
    end

    global mExt = Model(solver = GurobiSolver(NumericFocus = 3));
    startTE = time();
    global mExt = extForm(1, 0, [[],bData.bInv,[]], 1, τ, Δt, T, fData, bData, dData, pDistr, true);
    solve(mExt);
    elapsedTE = time() - startTE;
    mExtObj = getobjectivevalue(mExt);

    startT = time();
    cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist = solveMain(τ, T, Δt, N, false, false, 100, 100);
    elapsedT = time() - startT;
    push!(data,(T,elapsedT,cutDict,LBHist,UBHist,UBuHist,UBlHist));
end
