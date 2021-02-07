# load functions and packages

using JuMP, Ipopt, Gurobi, Combinatorics, JLD, HDF5, DelimitedFiles, Statistics, LinearAlgebra;

include("def.jl");
include("readin.jl");
include("main.jl");
include("detForm.jl");
include("auxiliary.jl");
include("forwardPass.jl");
include("backwardPass.jl");
