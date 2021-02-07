# test breaking a component and examine the deterministic costs
addprocs(30);
@everywhere include("loadMod.jl");
@everywhere const GUROBI_ENV = Gurobi.Env();

N = 5;

caseList = [13,33,123];
Δt = 0.25;
iterMax = 20;
T = 24;

for ci in 1:length(caseList)
    dataList = Dict();
    τ = Int64(1/6*T);

    for ω in keys(pDistr.ωDistrn)
        for j in procs()
            remotecall_fetch(readInData,j,ci,caseList,T);
            if length(ω) == 1
                remotecall_fetch(breakComponent, j, fData, ω, "g");
            else
                remotecall_fetch(breakComponent, j, fData, ω, "l");
            end
        end

        # break the component and then calculate the deterministic cost
        detSol,detObj = detBuild(Δt, T, fData, bData, dData);
        dataList[ω] = [detSol,detOb];
    end

    save("breakResults_$(ci).jld","data",dataList);
end
