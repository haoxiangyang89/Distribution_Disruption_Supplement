# generate samples for detTest and stabilityTest
using Distributed;
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

Δt = 0.25;
N = 5;
NN = 1000;
NN2 = 5000;
NNTest = 600;

caseList = [13,33,123];
dataDict = Dict();
dataDict2 = Dict();
dataDictTest = Dict();

for ci in 1:length(caseList)
    pathDict = Dict();
    TList = [24,36,48,72,96];
    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T,τ);
        end

        pathListData = pmap(i -> simuPath(T,pDistr), 1:NN);
        pathDict[T] = Dict();
        for i in 1:length(pathListData)
            pathDict[T][i] = pathListData[i];
        end
    end
    dataDict[ci] = pathDict;
end
# reverse the old generated sample
for ci in 1:2
    pathDict = Dict();
    TList = [24,36,48,72,96];
    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T,τ);
        end

        pathDict[T] = reverseScen(dataDict[ci][T],τ,pDistr);
    end
    dataDict[ci] = pathDict;
end
save("pathHist_1000.jld","pathDict",dataDict);

for ci in 1:length(caseList)
    pathDict = Dict();
    TList = [24,36,48,72,96];
    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T,τ);
        end

        pathListData = pmap(i -> simuPath(T,pDistr), 1:NNTest);
        pathDict[T] = Dict();
        for i in 1:length(pathListData)
            pathDict[T][i] = pathListData[i];
        end
    end
    dataDictTest[ci] = pathDict;
end
# reverse the old generated sample
for ci in 1:2
    pathDict = Dict();
    TList = [24,36,48,72,96];
    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T,τ);
        end

        pathDict[T] = reverseScen(dataDictTest[ci][T],τ,pDistr);
    end
    dataDictTest[ci] = pathDict;
end

save("pathHist_600.jld","pathDict",dataDictTest);

for ci in 1:length(caseList)
    pathDict = Dict();
    TList = [24,36,48,72,96];
    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T,τ);
        end

        pathListData = pmap(i -> simuPath(T,pDistr), 1:NN2);
        pathDict[T] = Dict();
        for i in 1:length(pathListData)
            pathDict[T][i] = pathListData[i];
        end
    end
    dataDict2[ci] = pathDict;
end
for ci in 1:2
    pathDict = Dict();
    TList = [24,36,48,72,96];
    for T in TList
        τ = Int64(1/6*T);
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T,τ);
        end

        pathDict[T] = reverseScen(dataDict2[ci][T],τ,pDistr);
    end
    dataDict2[ci] = pathDict;
end

save("pathHist_5000.jld","pathDict",dataDict2);
