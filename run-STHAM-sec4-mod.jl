######################################################################################################
#
# Julia code
#
# Executes the method described in Sec. 4 in the paper
#       "Self-Consistent Transport in Heterogeneous Agent Models" [STHAM] by Andrew Lyasoff
#
# Prrovides an alternative solution to the example from the paper
#       Krusell, Per, and Anthony Smith. (1998). Income and wealth heterogeneity in the macroeconomy.
#                         Journal of Political Economy 106 867-896.
#                
# Supplements the paper
#       "Self-Consistent Transport in Heterogeneous Agent Models" [STHAM] by Andrew Lyasoff
#
# Copyright ©2019-2025 Andrew Lyasoff <mathema@lyasoff.net>
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
    using Interpolations
    using Roots
    using Gnuplot
    using Printf
    using ForwardDiff
    using NLsolve
    using NLopt
    import Base.Threads.@spawn
    include("functions-STHAM-sec4.jl");
    include("ini-setup-STHAM-sec4-mod.jl");  # changes β in the original "ini-setup-STHAM-sec4.jl"
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
    #grid_ave=collect(0.4:0.0125:1.0);
    grid_ave=collect(0.4:0.01:1.0);
    len_grid=length(grid_ave)
end


#initial step
@time out1 = period_Tm1(grid_ave, 5.5, [0.5 for i=1:10], 1.0e-5, 1.0e-4, 1.0e-14, α, β, δ, NN, XX, AS, IS, ISpd, atpm, itpm, 15, 30);
#  4.805775 seconds (34.02 M allocations: 1.511 GiB, 6.46% gc time, 1285.51% compilation time)


begin
    result = out1[1];
    interp = out1[2];
    K_prev_spl = out1[3];
    sol_prev_spl = out1[4];
end;

#may be skipped; needed only to force Julia to compile the main function
@time out0=period_T(4, 2, grid_ave, interp, K_prev_spl, sol_prev_spl, 1.0e-5, 1.0e-5, 1.0e-4, 1.0e-12, α, β, δ, NN, XX, AS, IS, ISpd, atpm, itpm, 20, 50, 40);
#  3.689246 seconds (221.36 M allocations: 8.557 GiB, 27.15% gc time, 356.34% compilation time)

#the main program
@time out2 = period_T(1000, 2, grid_ave, interp, K_prev_spl, sol_prev_spl, 1.0e-5, 1.0e-5, 1.0e-4, 1.0e-12, α, β, δ, NN, XX, AS, IS, ISpd, atpm, itpm, 20, 50, 40);
# 131.083352 seconds (9.98 G allocations: 371.299 GiB, 32.50% gc time)

#=
EXPLANATION OF out2:

out2[1] = total number of completed iterations
out2[2] = results from the last iteration
out2[3] = results from next to last iterations
out2[2][x][i] = solution from the last iteration at aggregate state x and grid point (average) [i]
       output from 'solve_loc'
out2[2][x][i][1] = number of iterations between steps 3-7  ::Int64
out2[2][x][i][2] = future distribution mismatch  ::Float64
out2[2][x][i][3] = future distribution (as total average) in high and low states ::Vector{Float64}
              output from 'find_K02:
out2[2][x][i][4=end] = local solution::Tuple{Int64, Float64, Float64, Tuple{Bool, Vector{Float64}}}
out2[2][x][i][4][1] = number of iterations between steps 3-6
out2[2][x][i][4][2] = convergence in the search for capital
out2[2][x][i][4][3] = installed capital
                        output from 'get_cs'
out2[2][x][i][4][4] = complete solution to the system
out2[2][x][i][4][4][1] = solution exists::Bool
out2[2][x][i][4][4][2][9:10] = portfolio lines intercepts in high and low states
out2[2][x][i][4][4][2][1:8]  = furure consumption lines intercepts 4*(u-1)+2*(y-1)+v
                                   g(u,y,v) corresponds to entry 4*(u-1)+2*(y-1)+v
                        [g(1,1,1)=1, g(1,1,2)=2, g(1,2,1)=3, g(1,2,2)=4,
                                   g(2,1,1)=5, g(2,1,2)=6, g(2,2,1)=7, g(2,2,2)=8]
                                          y=1:[1,2,5,6], y=2:[3,4,7,8]   

ALTERNATIVE TRANSCRIPTION: 

out2[2][aggregate_state][point_on_the_grid]
     [
     iter_no::Int64, dist_mismatch::Float64, future_dist::Vector{Float64},
     local_solution[
               iter_to_find_K::Int64, convergence_for_K::Float64, last_K::Float64,
               solution[
                   solution exists::Bool,
                   intercepts_for_future_cons_and_portfolios::Vector{Float64}
                       ]
                  ]
     ]

=#


#=
save the ouput for later use
=#
serialize("output-STHAM-sec4-mod.jls", out2);
