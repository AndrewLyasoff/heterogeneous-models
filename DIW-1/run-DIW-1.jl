######################################################################################################
#
# Julia code
#
# Executes method described in Sec. 4 in the paper
#       "Another look at the distribution of income and wealth in the Macroeconomy" (DIW) by Andrew Lyasoff
#
# Prrovides an alternative solution to the example from the paper
#       Krusell, Per, and Anthony Smith. (1998). Income and wealth heterogeneity in the macroeconomy.
#                         Journal of Political Economy 106 867-896.
#                
# Supplements the paper "Another look at the distribution of income and wealth in the Macroeconomy" (DIW)
#                                        by Andrew Lyasoff (www.andrewlyasoff.com)
#
# Copyright © 2019-2022 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################

# test the number of avaiulable CPU
println(" available CPUs: ", length(Sys.cpu_info()) )


@time begin
    using SpecialFunctions
    using DataFrames, GLM
    using FastGaussQuadrature
    using Serialization, FileIO
    using LinearAlgebra
    using Polynomials
    using Dierckx
    using Roots
    using Gnuplot
    using Printf
    using ForwardDiff
    using NLsolve
    using NLopt
    import Base.Threads.@spawn
    include("functions-DIW-1.jl");
    include("ini-setup-DIW-1.jl"); # include("ini-setup-DIW-1-mod.jl") to change β from original
    gqx, gqw = gausslegendre( 100 );
    GQX, GQW = gausslegendre( 10_000 );
end;


begin
    IS=[1,2]; # idiosyncratic states
    AS=[1,2]; # productivity states
    RRR=1.0;  # risk aversion
    ISpd=[[0.96, 0.04], [0.9, 0.1]] #same as PJL' (population distribution over employment in high and low)
end;

#creating grids
begin
    xs=collect(0.5:0.05:1.0);
    ys=collect(0.8:0.02375:0.99);
    lx=length(xs);
    ly=length(ys);
    lx,ly
end


#initial step
@time out1 = period_Tm1(xs, ys, 10.5, 100.0, [0.5 for i=1:10], 1.0e-4, 1.0e-5, 1.0e-4, 1.0e-14, α, NN, XX, AS, IS, 15, 30);
#7.584958 seconds

begin
    result = out1[1];
    interp = out1[2];
    K_prev_spl = out1[3];
    sol_prev_spl = out1[4];
end;

#may be skipped; needed only to force Julia to compile the main function 
@time begin
    out0=period_T(4, 2, xs, ys, interp, K_prev_spl, 100.0, sol_prev_spl, 1.0e-4, 1.0e-5, 1.0e-5, 1.0e-4, 1.0e-12, α, NN, XX, AS, IS, 20, 30, 40);
end;
# 16.067788 seconds

#the main program
@time out2 = period_T(1000, 2, xs, ys, interp, K_prev_spl, 100.0, sol_prev_spl, 1.0e-4, 1.0e-5, 1.0e-5, 1.0e-4, 1.0e-12, α, NN, XX, AS, IS, 20, 30, 40);
#354.599245 seconds (34.45 G allocations: 752.972 GiB, 48.47% gc time)

#=
out2[1] = total number of completed iterations
out2[2] = results from the last iteration
out2[3] = results from next to last iterations
out2[2][x][i,j] = solution from the last iteration at aggregate state x and grid point [i,j]
out2[2][x][i,j][1] = number of iterations between steps 3-7  ::Int64
out2[2][x][i,j][2] = future distribution mismatch  ::Float64
out2[2][x][i,j][3] = future distribution in high and in low states ::Vector{Vector{Float64}}
out2[2][x][i,j][4=end] = ::Tuple{Int64, Float64, Float64, Tuple{Vector{Float64}, Int64, Vector{Float64}}}
out2[2][x][i,j][4][1] = number of iterations between steps 3-6
out2[2][x][i,j][4][2] = mismatch in the future distribution
out2[2][x][i,j][4][3] = installed capital
out2[2][x][i,j][4][4] = actual solution (list of 10 floats)
         out2[2][x][i,j][4][4][9:10] = portfolio lines intercepts in high and low states
         out2[2][x][i,j][4][4][1:8]  = furure consumption lines intercepts 4*(u-1)+2*(y-1)+v
                                   g(s,ξ,σ) corresponds to entry 4*(s-1)+2*(ξ-1)+σ
                        [g(1,1,1)=1, g(1,1,2)=2, g(1,2,1)=3, g(1,2,2)=4,
                                   g(2,1,1)=5, g(2,1,2)=6, g(2,2,1)=7, g(2,2,2)=8]
                                          y=1:[1,2,5,6], y=2:[3,4,7,8]   
=#


#serialize("output-DIW-1-mod.jls", out2);
serialize("output-DIW-1.jls", out2);
