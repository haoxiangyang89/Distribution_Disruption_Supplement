# main program structure construction

function solveMain(T, Δt, N, allGen = false, qpopt = false, iterMin = 100,
    iterMax = 1000, cutDict = Dict(), ubGen = false, ubM = 200, hardened = [], simuRule = 0, simuSample = Dict())
    # readin data and execute the SDDP algorithm

    UB = 9999999999;
    LB = -9999999999;
    iterNo = 0;
    keepIter = true;
    LBHist = [];
    UBHist = [];
    UBuHist = [];
    UBlHist = [];
    timeHist = [];
    # initialize the cutDict
    for j in procs()
        remotecall_fetch(cutIni,j,cutDict);
    end
    sampleCounter = 0;

    # while the termination criteria is not met
    while (keepIter)&(iterNo <= iterMax)
        iterNo += 1;
        iterStart = time();
        if simuRule == 0
            if simuSample == Dict()
                trialPaths,currentLB,currentUBDict = exeForward(T, Δt, N, qpopt, Dict(), hardened);
            else
                pathIter = Dict();
                for sNo in 1:N
                    pathIter[sNo] = simuSample[sampleCounter+sNo];
                end
                trialPaths,currentLB,currentUBDict = exeForward(T, Δt, N, qpopt, pathIter, hardened);
                sampleCounter += N;
            end
        else
            trialPaths,currentLB,currentUBDict = exeForward_simuOpt(T, Δt, N, iterNo, qpopt, hardened);
        end
        push!(LBHist,currentLB);

        # forward pass: obtain the trial paths
        if !(ubGen)
            currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
        else
            trialPathsUB,currentLBUB,currentUBDict = exeForward(T, Δt, ubM, qpopt, Dict(), hardened);
            currentUBList = [currentUBDict[ubkey] for ubkey in keys(currentUBDict)];
        end
        push!(UBHist,mean(currentUBList));
        currentUBu = mean(currentUBList) + 1.96*std(currentUBList);
        push!(UBuHist,currentUBu);
        currentUBl = mean(currentUBList) - 1.96*std(currentUBList);
        push!(UBlHist,currentUBl);

        if (currentLB <= currentUBu)&(currentLB >= currentUBl)&(iterNo >= iterMin)
            keepIter = false;
        else
            if allGen
                exeBackwardAll(T, Δt, trialPaths, qpopt, hardened);
            else
                exeBackward(T, Δt, trialPaths, qpopt, hardened);
            end
        end
        iterElapsed = time() - iterStart;
        push!(timeHist,iterElapsed);
        println("========= Iteration $(iterNo) Finished, LB = $(round(currentLB, digits = 2)),
            UB = [$(round(currentUBl,digits = 2)),$(round(currentUBu, digits = 2))], Time = $(iterElapsed) sec. =========")
    end
    return cutDict,LBHist,UBHist,UBuHist,UBlHist,timeHist;
end

function preGen(T, Δt, N, iterMax, qpopt = false, cutDict = Dict(), hardened = [])
    # pregenerate cuts

    UB = 9999999999;
    LB = -9999999999;
    iterNo = 0;
    keepIter = true;
    LBHist = [];
    UBHist = [];
    UBuHist = [];
    UBlHist = [];
    timeHist = [];
    # initialize the cutDict
    for j in procs()
        remotecall_fetch(cutIni,j,cutDict);
    end

    while iterNo <= iterMax
        iterNo += 1;
        trialPaths,currentLB,currentUBDict = exeForward_last(T, Δt, N, qpopt, hardened);
        exeBackward_last(T, Δt, trialPaths, qpopt, hardened);
    end

    return cutDict;
end
