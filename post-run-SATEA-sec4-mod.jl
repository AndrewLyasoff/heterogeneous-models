######################################################################################################
#
# Julia code
#
# Illustrates the method described in Sec. 4 in the paper
#       "Self-Aware Transport of Economic Agents" [SATEA] by Andrew Lyasoff
#
# Provides an alternative solution to the example from the paper
#       Krusell, Per, and Anthony Smith. (1998). "Income and wealth heterogeneity in the macroeconomy."
#                         Journal of Political Economy 106 867-896.
#                
# Supplements the paper
#                "The Time-Interlaced Self-Consistent Master System of Heterogeneous-Agent Models" [SATEA] by Andrew Lyasoff
#
#
# Copyright © 2019-2025 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################


# not needed if run-SATEA-sec4-mod.jl is executed
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
    include("functions-SATEA-sec4.jl");
    include("ini-setup-SATEA-sec4-mod.jl"); # include("ini-setup-SATEA-sec4-mod.jl") to change β from original
    gqx, gqw = gausslegendre( 100 );
    GQX, GQW = gausslegendre( 10_000 );
    IS=[1,2]; # idiosyncratic states
    AS=[1,2]; # productivity states
    RRR=1.0;  # risk aversion
    ISpd=[[0.96, 0.04], [0.9, 0.1]] #same as PJL' (population distribution over employment in high and low)
    grid_ave=collect(0.4:0.01:1.0);
    len_grid=length(grid_ave)
end

# not needed if run-SATEA-sec4-mod.jl is executed
out2=deserialize("output-SATEA-sec4-mod.jls");


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

begin
    result_prev=out2[3];
    result = out2[2];
    interp = [[Spline1D(grid_ave, [result[x][iii][4][4][2][8+u] for iii=1:len_grid], k=3, bc = "extrapolate") for u in IS] for x in AS];
    K_prev_spl = [Spline1D(grid_ave, [result[x][iii][4][3] for  iii=1:len_grid], k=3, bc = "extrapolate") for x in AS];
    sol_prev_spl = [[Spline1D(grid_ave, [result[x][iii][4][4][2][inx] for iii=1:len_grid], k=3, bc = "extrapolate") for inx=1:length(result[x][1][4][4][2])] for x in AS];
end;


# convergence of portfolio intercepts
maximum([maximum([maximum([abs(result[x][iii][4][4][2][8+u]-result_prev[x][iii][4][4][2][8+u]) for iii=1:len_grid]) for u in IS]) for x in AS])
# 3.552713678800501e-15


#convergence of future consumption intercepts
maximum([maximum([maximum([abs(result[x][iii][4][4][2][ii]-result_prev[x][iii][4][4][2][ii]) for iii=1:len_grid]) for ii=1:8]) for x in AS])
# 6.938893903907228e-18


#convergence of installed capital
maximum([maximum([abs(result[x][iii][4][3]-result_prev[x][iii][4][3]) for iii=1:len_grid]) for x in AS])
# 3.552713678800501e-15

#convergence of future distributions
maximum([maximum([maximum(abs.(vcat((result[x][iii][3]-result_prev[x][iii][3])...))) for iii=1:len_grid]) for x in AS])
# 0.0

# maximum mismatch in future distribution and guessed ones in high and low states
[maximum(abs.([result[x][iii][2] for iii=1:len_grid])) for x in AS]
#=
2-element Vector{Float64}:
 0.0
 9.96594028390696e-6
=#

#largest consumption transfer line intercept
maximum([maximum([maximum([result[x][iii][4][4][2][j] for iii=1:len_grid]) for j=1:8]) for x in AS])
# 0.03242516394526393

#smallest consumption intercept
minimum([minimum([minimum([result[x][iii][4][4][2][j] for iii=1:len_grid]) for j=1:8]) for x in AS])
# -0.05426617609947015

# largest distance between portfolio lines intercepts in high and low state for employed
maximum([abs(result[1][iii][4][4][2][9]-result[2][iii][4][4][2][9]) for iii=1:len_grid])
# 0.012001269093659062

# largest distance between portfolio lines intercepts in high and low state for unemployed
maximum([abs(result[1][iii][4][4][2][10]-result[2][iii][4][4][2][10]) for iii=1:len_grid])
# 0.28686579418603486


# portfolio intercepts in high state
## provides the left plot in Figure 10 (not in SATEA)
begin
    #plot_grid=collect(grid_ave[1]:0.005:grid_ave[end])
    plot_grid=collect(0.5:0.005:1.0)
    @gp plot_grid [interp[1][1](arg) for arg in plot_grid] "w l t '' lw 2.5 lt rgb 'black'";
    @gp :- plot_grid [interp[1][2](arg) for arg in plot_grid] "w p lt rgb 'black' pt 7 ps 0.6 t ''";
end


# portfolio intercepts in how state
## provides the right plot in Figure 10 (not in SATEA)
begin
    #plot_grid=collect(grid_ave[1]:0.005:grid_ave[end])
    plot_grid=collect(0.5:0.005:1.0)
    @gp plot_grid [interp[2][1](arg) for arg in plot_grid] "w l t '' lw 2.5 lt rgb 'black'";
    @gp :- plot_grid [interp[2][2](arg) for arg in plot_grid] "w p lt rgb 'black' pt 7 ps 0.6 t ''";
end


# installed capital in high and low states
## provides the plot in Figure 11 (not in SATEA)
begin
    #plot_grid=collect(grid_ave[1]:0.005:grid_ave[end])
    plot_grid=collect(0.5:0.005:1.0)
    @gp plot_grid [K_prev_spl[1](arg) for arg in plot_grid] "w l t '' lw 2.5 lt rgb 'black'";
    @gp :- plot_grid [K_prev_spl[2](arg) for arg in plot_grid] "w p lt rgb 'black' pt 7 ps 0.6 t ''";
end

# uniform distance between installed capital in high and low states
maximum(abs.([K_prev_spl[1](arg)-K_prev_spl[2](arg) for arg in plot_grid]))
# 0.062207573066050514

transport_a(1, 2, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, 0.7)
# 0.6919438245183468

## transport of the average into state 1
## provides the left plot in Figure 12 (not in SATEA)
begin
    #plot_grid=collect(grid_ave[1]:0.005:grid_ave[end])
    plot_grid=collect(0.5:0.005:1.0)
    T11=[transport_a(1, 1, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg) for arg in plot_grid];
    T21=[transport_a(2, 1, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg) for arg in plot_grid];
    @gp plot_grid T11 "w l t '' lw 2.5 lt rgb 'black'";
    @gp :- plot_grid T21 "w p lt rgb 'black' pt 7 ps 0.6 t ''";
end


## transport of the average into state 2
## provides the right plot in Figure 12 (not in SATEA)
begin
    #plot_grid=collect(grid_ave[1]:0.005:grid_ave[end])
    plot_grid=collect(0.5:0.005:1.0)
    T12=[transport_a(1, 2, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg) for arg in plot_grid];
    T22=[transport_a(2, 2, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg) for arg in plot_grid];
    @gp plot_grid T12 "w l t '' lw 2.5 lt rgb 'black'";
    @gp :- plot_grid T22 "w p lt rgb 'black' pt 7 ps 0.6 t ''";
end

# largest uniform distance between transport of capital in the last two plots
(maximum(abs.((T11).-T21)), maximum(abs.((T12).-T22)))
# (0.005499037773428639, 0.005478140439179824)

#uniform distance between solid lines
maximum(abs.((T11).-T12))
# 0.010364586765959305

## transport in terms of capital into high state
## provides the left plot in Figure 13 (not in SATEA)
begin
    #plot_grid=collect(grid_ave[1]:0.005:grid_ave[end])
    plot_grid=collect(0.5:0.005:1.0);
    fut_a_1_1=[transport_a(1, 1, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg) for arg in plot_grid];
    fut_a_2_1=[transport_a(2, 1, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg) for arg in plot_grid];
    K_1_grid=[K_prev_spl[1](arg) for arg in plot_grid];
    K_2_grid=[K_prev_spl[2](arg) for arg in plot_grid];
    fut_K_1_1=[K_prev_spl[1](arg) for arg in fut_a_1_1];
    fut_K_2_1=[K_prev_spl[1](arg) for arg in fut_a_2_1];
    @gp K_1_grid fut_K_1_1 "w l t '' lw 2.5 lt rgb 'black'";
    @gp :- K_2_grid fut_K_2_1 "w p lt rgb 'black' pt 7 ps 0.9 t ''";
end


# largest distance between solid and dotted lines
begin
    fspl1=Spline1D(K_1_grid, fut_K_1_1, k=3, bc = "extrapolate");
    fspl2=Spline1D(K_2_grid, fut_K_2_1, k=3, bc = "extrapolate");
    lst1=[fspl1(arg) for arg in max(K_1_grid[1],K_2_grid[1]):0.001:min(K_1_grid[end],K_2_grid[end])];
    lst2=[fspl2(arg) for arg in max(K_1_grid[1],K_2_grid[1]):0.001:min(K_1_grid[end],K_2_grid[end])];
    maximum(abs.((lst1).-lst2))
end
# 8.625025005359888e-5


## transport in terms of capital into low state
## provides the right plot in Figure 13 (not in SATEA)
begin
    #plot_grid=collect(grid_ave[1]:0.005:grid_ave[end])
    plot_grid=collect(0.5:0.005:1.0);
    fut_a_1_2=[transport_a(1, 2, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg) for arg in plot_grid];
    fut_a_2_2=[transport_a(2, 2, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg) for arg in plot_grid];
    K_1_grid=[K_prev_spl[1](arg) for arg in plot_grid];
    K_2_grid=[K_prev_spl[2](arg) for arg in plot_grid];
    fut_K_1_2=[K_prev_spl[2](arg) for arg in fut_a_1_2];
    fut_K_2_2=[K_prev_spl[2](arg) for arg in fut_a_2_2];
    @gp K_1_grid fut_K_1_2 "w l t '' lw 2.5 lt rgb 'black'";
    @gp :- K_2_grid fut_K_2_2 "w p lt rgb 'black' pt 7 ps 0.9 t ''";
end

# largest distance between solid and dotted lines
begin
    fsspl1=Spline1D(K_1_grid, fut_K_1_2, k=3, bc = "extrapolate");
    fsspl2=Spline1D(K_2_grid, fut_K_2_2, k=3, bc = "extrapolate");
    lsst1=[fsspl1(arg) for arg in max(K_1_grid[1],K_2_grid[1]):0.001:min(K_1_grid[end],K_2_grid[end])];
    lsst2=[fsspl2(arg) for arg in max(K_1_grid[1],K_2_grid[1]):0.001:min(K_1_grid[end],K_2_grid[end])];
    maximum(abs.((lsst1).-lsst2))
end
# 7.699665445715098e-5

# largest distance between the two solid lines
maximum(abs.((lsst1).-lst1))
# 0.04847944993530273 

# examples of multivariate transport
last=transport_mult(1, 2, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, [0.8,0.6])
#=
2-element Vector{Float64}:
 0.787955084364198
0.6886814525447722
=#

# test agreement with transport of averages
transport_a(1, 2, α, β, δ, IS, ISpd, itpm, NN, XX, K_prev_spl, sol_prev_spl, ISpd[1]'*[0.8,0.6])
# 0.7780277211822555 

#which agrees with
ISpd[2]'*last
# 0.7780277211822555

#=
largest deviation in the kernel condition with comsumption =  the total average
=#
maximum([maximum([abs(kernel_err(x, u, α, β, δ, AS, IS, atpm, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg, arg)) for arg in grid_ave]) for u in IS, x in AS])
# 0.0014318389980028101

#=
largest deviation in the kernel condition with consumption = the borrowing limit
=#
maximum([maximum([abs(kernel_err(x, u, α,  β, δ, AS, IS, atpm, itpm, NN, XX, K_prev_spl, sol_prev_spl, arg, -((1-β)/β)*sol_prev_spl[x][8+u](arg))) for arg in grid_ave]) for u in IS, x in AS])
# 0.002447423058417275 


# simulation of group behavior 

begin
    using Random
    rng=MersenneTwister(42)
    #rng=RandomDevice()
end;


begin
    AR,A=gen_fut_dist(1, [0.8,0.7], 1100000, α, β, δ, IS, ISpd, itpm, atpm, NN, XX, K_prev_spl, sol_prev_spl);
    loc_var=filter(x->(x==2),AR);
    length(loc_var)/1100000 # should be close to 0.5 (steady state probability for low state)
end
# 0.5009827272727273


# plot the state of the population in terms of consumption in the last 10,000 periods
## provides the left plot in Figure 16
begin
    @gp [vec[1] for vec in A[end-10000:end]] [vec[2] for vec in A[end-10000:end]] "w p lt -1 pt 5 ps 0.1 t ''"
end


## transform consumption into financial wealth
## VW1 = investment of employed    VW2 = investment of unemployed
begin
    VW1=[sol_prev_spl[AR[ii]][9](ISpd[AR[ii]]'*A[ii])+(β/(1-β))*(A[ii][1]) for ii=(length(A)-10000):length(A)];
    VW2=[sol_prev_spl[AR[ii]][10](ISpd[AR[ii]]'*A[ii])+(β/(1-β))*(A[ii][2]) for ii=(length(A)-10000):length(A)];
end;


# plot the state of the population in terms of investment in the last 10,000 periods
## provides the right plot in Figure 16
begin
    @gp VW1 VW2 "w p lt -1 pt 5 ps 0.1 t ''"
end



##########################

#=
largest deviation in the kernel condition with comsumption =  employment class average
=#
maximum([maximum([abs(kernel_err(AR[ii], u, α, β, δ, AS, IS, atpm, itpm, NN, XX, K_prev_spl, sol_prev_spl, ISpd[AR[ii]]'*A[ii], A[ii][u])) for ii=(length(A)-10000):length(A)]) for u in IS])
# 0.0012384702068934939

#=
largest deviation in the kernel condition with consumption = the borrowing limit
=#
maximum([maximum([abs(kernel_err(x, u, α,  β, δ, AS, IS, atpm, itpm, NN, XX, K_prev_spl, sol_prev_spl, ISpd[AR[ii]]'*A[ii], -((1-β)/β)*sol_prev_spl[x][8+u](ISpd[AR[ii]]'*A[ii]))) for  ii=(length(A)-10000):length(A)]) for u in IS, x in AS])
# 0.0019844303932907703

