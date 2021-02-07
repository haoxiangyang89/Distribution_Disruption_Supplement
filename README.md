# Distribution_Disruption_Supplement
Code and supplement data repository for the paper "Optimal Power Flow in Distribution Networks under N-1 Disruptions: A Multi-stage Stochastic Programming Approach"

---
Main body:
 * loadMod.jl: scripts to load required functions
 * def.jl: definitions of data structure
 * readin.jl: functions to read in power network data
 * main.jl: functions to execute the SDDP algorithm
 * forwardPass.jl: execute backwardPass
 * backwardPass.jl: execute backwardPass
 * auxiliary.jl: auxiliary functions to process samples and change data structures
 * detForm.jl: formulate the deterministic policy
 * extForm.jl: formulate the extensive formulation (only for small samples)
---
Test functions:
 * allGenTest.jl/dOnlyTest.jl: generate results between DOnly and GenAll (Section 4.1.1)
 * NTest.jl: generate results comparing different N^p (Section 4.1.2)
 * detTest.jl/stochNomTest.jl/stabilityTest.jl: generate results comparing deterministic policy vs. stochastic policy (Section 4.2 Table 3/4)
 * ubTest.jl: generate results of statistical upper bound (Section 4.2 Figure 5)
 * butilTest.jl: generate results of battery capacity and utilization (Section 4.3.1)
 * hardeningTest.jl: generate results of component hardening (Section 4.3.2)
 * tauTest.jl: generate results of recovery phase duration (Section 4.3.3)
 * multiTest.jl: generate results of multiple-component disruptin (Section 4.3.4)
---
Results data:
 * dOnly_lb/dOnly_time.csv: Table 1, Figure 3
 * N_lb/N_time.csv: Table 2, Figure 4
 * detResults.jld: Table 3
 * stabilityResults.jld: Table 4
 * ub.csv/UBResults.jld: Figure 5
 * butilResults.jld: Table 5
 * hardResults.jld: Table 6
 * tauResults.jld: Table 7
 * tauResults_specific.jld: Table 8
 * multiResults.jld: Table 9
