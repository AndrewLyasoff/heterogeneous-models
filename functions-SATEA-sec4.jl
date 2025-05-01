######################################################################################################
#
# Julia code
#
# Functions to implement the method described in Sec. 4 in the paper
#       "The Time-Interlaced Self-Consistent Master System of Heterogeneous-Agent Models" [SATEA] by Andrew Lyasoff
#
# The code provides an alternative solution to the example from the paper
#       Krusell, Per, and Anthony Smith. (1998). Income and wealth heterogeneity in the macroeconomy.
#                         Journal of Political Economy 106 867-896.
#                
# This code supplements the paper "The Time-Interlaced Self-Consistent Master System of Heterogeneous-Agent Models" [SATEA]
#                                                         by Andrew Lyasoff
#
# Copyright ©2019-2025 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################


#=
return on capital as a function of installed capital (arg_x)

XX is the list of aggregate shocks (in various aggregate states)
LL is the list of aggregate supplied labor
α is the risk aversion
x is the present state (superfluous in the next function)
y is the future state
=#

function ρ(y::Int64, α::Float64, LL::Array{Float64,1}, XX::Array{Float64,1})
    return (arg_K -> XX[y]*α*(arg_K/LL[y])^(α-1))
end

#
# wages in the future period as function of the capital installed in the present period
#

function w(y::Int64, α::Float64, LL::Array{Float64,1}, XX::Array{Float64,1})
    return ( arg_K -> XX[y]*(1-α)*(arg_K/LL[y])^α )
end

#=
NB: all slopes are in exact form; only the intercepts are unknown

UNKNOWNS (t and x are given):

g(u,y,v) corresponds to entry 4*(u-1)+2*(y-1)+v  in the list of 10 unknowns  u=1:2, y=1:2, v=1:2 | belongs to t+1

a(x,u) corresponds to entry 8+u  u=1:2 in the list of 10 unknowns | belongs to t

NB: the intercepts of the future consumption lines are the first 8 unknowns out of 10
NB: the intercepts of the present portfolio lines are the last 2 unknowns out of 10 

GIVENS:

a(y,v) , y=1:2, v=1:2 | belongs to t+1

NB: the intercepts of the portfolio lines for the next period are given
and depend on the future productivity state and the future employment state

cs::Float64 is c* -- see 4.1

=#

#=

system of first order Lagrange condition (FOLC) to solve -- see (4.14) and (4.16)
10 equations with 10 unknowns

written as a function from ℝ¹⁰ into ℝ¹⁰

period t capital and period t+1 portfolio intercepts are given parameters

=#
function sstm(n::Int64,  xxx::Int64, KK::Float64, fut_a_ini::Array{Array{Float64, 1}, 1}, α::Float64, β::Float64, δ::Float64,  NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1}, atpm::Matrix{Float64}, itpm::Array{Float64, 4})
    return (UUU->(vcat(vcat(vcat([[[(fut_a_ini[y][v] + UUU[4*(u-1)+2*(y-1)+v]*((1-β^n)/(1-β)) - UUU[8+u]*(ρ(y, α, NN, XX)(KK)+1-δ)-SS[v]*w(y, α, NN, XX)(KK)) for v in IS] for y in AS] for u in IS]...)...) , [(sum([-(1/(ρ(y, α, NN, XX)(KK)+1-δ))*(UUU[4*(u-1)+2*(y-1)+v])*atpm[xxx,y]*itpm[u,v,xxx,y] for v in IS, y in AS])) for u in IS ] )))
end


# the left side of the market clearing condition -- see (4.17)
function clear_mkt(n::Int64, xxx::Int64, sol_aaa::Array{Float64, 1}, c_ave::Float64, ISpd::Array{Array{Float64,1},1})
    return ((sum([β^i for i=1:n])*c_ave) + (ISpd[xxx]')*sol_aaa)
end


#=
solves the system composed of (4.14) and (4.16) for the 10 unknown intercepts (8 for the future consumption lines
                                              and 2 for the present portfolio lines)
=#
function solve_sstm(n::Int64, xxx::Int64, solve_ini_guess::Array{Float64, 1}, solve_fac::Float64, solve_ftol::Float64, KK::Float64, fut_a_ini::Array{Array{Float64, 1}, 1}, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1})
    local lcl_sln, test::Bool;
    test = false;
    lcl_sln=nlsolve(sstm(n, xxx, KK, fut_a_ini, α, β, δ, NN, XX, AS, IS, atpm, itpm), solve_ini_guess, factor=solve_fac, ftol=solve_ftol);
    test=converged(lcl_sln);
    if test
        return (css,1,lcl_sln.zero)
    else
        return false
    end
end


#=

Returns the solution to (4.14) & (4.16) for the present portfolio intercepts and
for the state transitions intercepts
(future portfolio intercepts are given)

=#
function get_cs(n::Int64, xxx::Int64, solve_ini_guess::Array{Float64, 1}, solve_fac::Float64, solve_ftol::Float64, KK::Float64, fut_a_ini::Array{Array{Float64, 1}, 1}, α::Float64,  β::Float64, δ::Float64, NN::Array{Float64,1}, XX::Array{Float64,1}, AS::Array{Int64,1}, IS::Array{Int64,1}, atpm::Matrix{Float64}, itpm::Array{Float64, 4}, kill_switch::Int64)
    local lcl_sln, test::Bool
    test = false;
    lcl_sln=nlsolve(sstm(n, xxx, KK, fut_a_ini, α, β, δ, NN, XX, AS, IS, atpm, itpm), solve_ini_guess, factor=solve_fac, ftol=solve_ftol);
    test=converged(lcl_sln);
    if test
        return (true,lcl_sln.zero)
    else
        return (false)
    end
end


#=

tâtonnement over capital:
performs iterations (2) through (4) in the generic backward step -- see [4.6]
(for fixed t ~ n, fixed productivity state xxx, and fixed population mean = c_ave)

returns:
the number of iterations performed, the market clearing achieved, the last K,
and the solution to (4.14) & (4.16)


=#
function find_K02(n::Int64, xxx::Int64, c_ave::Float64, K_ini::Float64, thresh_mc::Float64, solve_ini_guess1::Array{Float64, 1}, solve_fac::Float64, solve_ftol::Float64,  fut_a_ini::Array{Array{Float64, 1}, 1}, α::Float64, β::Float64, δ::Float64, NN::Array{Float64,1}, XX::Array{Float64,1}, IS::Array{Int64,1}, ISpd::Vector{Vector{Float64}}, atpm::Matrix{Float64}, itpm::Array{Float64, 4}, kill_switch::Int64, kill_switch2::Int64)
    local K_prev::Float64, K_last::Float64, mc_test::Float64, all_go::Bool, sol, iter::Int64;
    all_go=false;
    iter = 0;
    mc_test = Inf;
    K_last = K_ini;
    sol = get_cs(n, xxx, solve_ini_guess1, solve_fac, solve_ftol, K_last, fut_a_ini, α, β, δ, NN, XX, AS, IS, atpm, itpm, kill_switch);
    #println(sol[2])
    if (sol[1]) 
        K_last = clear_mkt(n, xxx, sol[2][9:10], c_ave, ISpd);
        if (K_last<1.0e-5) K_last=0.15 end;
        all_go = true;
        iter = 1;
    else
        println("solution falied at T = ", iter)
        all_go = false;
    end
    while (all_go & (iter<kill_switch2) & (mc_test>thresh_mc))
        sol=get_cs(n, xxx, solve_ini_guess1, solve_fac, solve_ftol, K_last, fut_a_ini, α, β, δ, NN, XX, AS, IS, atpm, itpm, kill_switch);
        if (sol[1])
            K_prev=K_last;
            K_last=clear_mkt(n, xxx, sol[2][9:10], c_ave, ISpd);
            all_go=true;
            iter+=1;
            if (K_last<1.0e-4)
                K_last=1.0e3;
                mc_test = Inf;
            else
                mc_test = abs(K_last-K_prev);
            end
        else
            all_go = false;
        end
    end
    if ((iter==kill_switch2)&(mc_test>thresh_mc))
        all_go=false;
        println("find_K02 reached the iterations maximum of ", kill_switch2)
    end;
    if all_go
        return (iter, mc_test, K_last, sol)
    else
        return false
    end
end



#=

tâtonnement over future distribution (an element of ℝ):

performs iterations (2) through (6) from the Metaprogram in [4.6]
meant to adjust the guess for the future distribution
(for fixed t ~ n, fixed productivity state xxx, and fixed population mean = c_ave)

fut_a_sp:: one spline for every employement state, for every aggreagare state (computed during the previous
iteration, i.e., in the future period) NB: represent portfolio slopes

K_prev_sp:: one spline for every aggregate state, represents total capital (given from the previous iter)

sol_prev_sp:: a list of 10 splines for every aggregate state * represents the slopes of the state transitions
and the present portfolios (given from the previous iter)

the present aggregate state = xxx and present total consumption mean = c_ave are fixed

the time period t ~ n is fixed

=#
function solve_loc(n::Int64, xxx::Int64, c_ave::Float64, ansatz_fut_ave::Array{Float64, 1}, fut_a_sp::Array{Array{Spline1D, 1}, 1}, K_prev_sp::Array{Spline1D, 1}, thresh_d::Float64, thresh_mc::Float64, sol_prev_sp::Array{Array{Spline1D, 1}, 1}, solve_fac::Float64, solve_ftol::Float64, α::Float64, β::Float64, δ::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1}, ISpd::Vector{Vector{Float64}}, atpm::Matrix{Float64}, itpm::Array{Float64, 4}, kill_switch::Int64, kill_switch2::Int64, kill_switch3::Int64)
    local  loc_sol, next_ave_iter::Vector{Float64}, distance_d::Float64, iter_no::Int64, K_ini::Float64, solve_ini_guess1::Array{Float64, 1}, no_alaram::Bool;
    fut_ave_guess = ansatz_fut_ave;
    fut_a = [[fut_a_sp[y][v](fut_ave_guess[y]) for v in IS] for y in AS];
    K_ini = K_prev_sp[xxx](c_ave); # step (2) 
    solve_ini_guess1 = [sol_prev_sp[xxx][inx](c_ave) for inx=1:length(sol_prev_sp[xxx])];
    loc_sol = find_K02(n, xxx, c_ave, K_ini, thresh_mc, solve_ini_guess1, solve_fac, solve_ftol, fut_a, α, β, δ, NN, XX, IS, ISpd, atpm, itpm, kill_switch, kill_switch2);
    if length(loc_sol) > 1
        no_alarm = true;
        next_ave_iter = [(β*(ρ(y,α,NN,XX)(loc_sol[3])+1-δ)*c_ave + sum([ISpd[xxx][u]*itpm[u,v,xxx,y]*loc_sol[4][2][4*(u-1)+2*(y-1)+v] for u in IS, v in IS])) for y in AS];
        dist_d = maximum(abs.(next_ave_iter.-fut_ave_guess));
        iter_no = 1;
    else
        no_alarm = false;
        dist_d = -Inf;
        iter_no = 2^50;
    end
    while no_alarm & (dist_d > thresh_d) & (iter_no < kill_switch3)
        fut_ave_guess = next_ave_iter;
        fut_a = [[fut_a_sp[y][v](fut_ave_guess[y]) for v in IS] for y in AS];
        loc_sol = find_K02(n, xxx, c_ave, K_ini, thresh_mc, solve_ini_guess1, solve_fac, solve_ftol, fut_a, α, β, δ, NN, XX, IS, ISpd, atpm, itpm, kill_switch, kill_switch2);
        if length(loc_sol) > 1
            no_alarm = true;
            next_ave_iter = [(β*(ρ(y,α,NN,XX)(loc_sol[3])+1-δ)*c_ave + sum([ISpd[xxx][u]*itpm[u,v,xxx,y]*loc_sol[4][2][4*(u-1)+2*(y-1)+v] for u in IS, v in IS])) for y in AS];
            iter_no+=1;
            dist_d = maximum(abs.(next_ave_iter.-fut_ave_guess));
        else
            no_alarm = false;
        end;
    end;
    if ((dist_d > thresh_d)&(iter_no == kill_switch3))
        no_alarm=false;
        println("solve_loc reached the iterations maximum of ", kill_switch3);
    end;
    if no_alarm
        return (iter_no, dist_d, next_ave_iter, loc_sol)
    else
        return false
    end
end

#=
completes the initial backward step for all points A* over a chosen grid:

OUTPUTS: the values for the 10 intercepts and the aggregate capital at ALL points on the grid
and for all productivity states; returns also the interpolated versions of those values (this redundacy is necssary
because splines cannot be dumped)

=#
function period_Tm1(grid_ave::Array{Float64, 1}, K_ini::Float64, solve_ini_guess1::Array{Float64, 1}, thresh_mc::Float64, solve_fac::Float64, solve_ftol::Float64, α::Float64, β::Float64, δ::Float64, NN::Array{Float64,1}, XX::Array{Float64,1}, AS::Array{Int64,1}, IS::Array{Int64,1}, ISpd::Vector{Vector{Float64}}, atpm::Matrix{Float64}, itpm::Array{Float64, 4}, kill_switch::Int64, kill_switch2::Int64)
    local fut_a::Vector{Vector{Float64}}, len_grid::Int64, alarm_off::Bool, loc_result, result1, result2, result_prev, result, fut_a_spl::Array{Array{Spline1D, 1}, 1}, K_prev_spl::Array{Spline1D, 1}, sol_prev_spl::Array{Array{Spline1D, 1}, 1}, test_K::Bool;
    test_K=false;
    len_grid=length(grid_ave);
    fut_a=[[0.0 for v in IS] for y in AS]; #no investment
    loc_result = [@spawn find_K02(1, 1, arg, K_ini, thresh_mc, solve_ini_guess1, solve_fac, solve_ftol, fut_a, α, β, δ, NN, XX, IS, ISpd, atpm, itpm, kill_switch, kill_switch2) for arg in grid_ave];
    #
    result1 = [fetch(loc_result[iii]) for iii=1:len_grid];
    #
    if (length(findall(x->(x==false) , result1)) > 0)
        prinln("solution failed at certain grid points")
    else
        for iii=1:len_grid
            if (abs(result1[iii][3]-clear_mkt(1,1,result1[iii][4][2][9:10],grid_ave[iii],ISpd))>0.001)
                test_K=true
            end
        end
    end
    if (test_K)
        println("state 1: wrong capital at n = 1")
    end
    test_K=false;
    #
    #
    if (length(findall(x->(x==false) , result1)) > 0)
        alarm_off = false;
        println("alarm 1 at period T = ", 1)
    else
        alarm_off = true;
        loc_result = [@spawn find_K02(1, 2, arg, K_ini, thresh_mc, solve_ini_guess1, solve_fac, solve_ftol, fut_a, α, β, δ, NN, XX, IS, ISpd, atpm, itpm, kill_switch, kill_switch2) for arg in grid_ave];
        result2 = [fetch(loc_result[iii]) for iii=1:len_grid];
        #
        for iii=1:len_grid
            if (abs(result2[iii][3]-clear_mkt(1,2,result2[iii][4][2][9:10],grid_ave[iii],ISpd))>0.001)
                test_K=true
            end
        end
        if (test_K)
            println("state 2: wrong capital at n = ", nn)
        end
        test_K=false;
        #
        if length(findall(x->(x==false) , result2)) > 0
            alarm_off = false;
            println("alarm 2 at period T = ", 1)
            return false
        else
            result=[result1,result2];
            fut_a_spl = [[Spline1D(grid_ave, [result[x][iii][4][2][8+u] for iii=1:len_grid], k=3, bc = "extrapolate") for u in IS] for x in AS];
            K_prev_spl = [Spline1D(grid_ave, [result[x][iii][3] for iii=1:len_grid], k=3, bc = "extrapolate") for x in AS];
            sol_prev_spl = [[Spline1D(grid_ave, [result[x][iii][4][2][inx] for iii=1:len_grid], k=3, bc = "extrapolate") for inx=1:length(result[x][1][4][2])] for x in AS];
            #cs_prev_spl = [Spline2D(xs, ys, [result[x][i,j][4][1] for i=1:lx, j=1:ly]) for x in AS];
            return (result, fut_a_spl, K_prev_spl, sol_prev_spl)
        end
    end
end

#=
repeats the generic backward step for T periods

EXPLANATION OF output:

output[1] = total number of completed iterations
output[2] = results from the last iteration
output[3] = results from next to last iterations
output[2][x][i] = solution from the last iteration at aggregate state x and grid point (average) [i]
       output from 'solve_loc'
output[2][x][i][1] = number of iterations between steps 3-7  ::Int64
output[2][x][i][2] = future distribution mismatch  ::Float64
output[2][x][i][3] = future distribution (as total average) in high and low states ::Vector{Float64}
              output from 'find_K02:
output[2][x][i][4=end] = local solution::Tuple{Int64, Float64, Float64, Tuple{Bool, Vector{Float64}}}
output[2][x][i][4][1] = number of iterations between steps 3-6
output[2][x][i][4][2] = convergence in the search for capital
output[2][x][i][4][3] = installed capital
                        output from 'get_cs'
output[2][x][i][4][4] = complete solution to the system
output[2][x][i][4][4][1] = solution exists::Bool
output[2][x][i][4][4][2][9:10] = portfolio lines intercepts in high and low states
output[2][x][i][4][4][2][1:8]  = furure consumption lines intercepts 4*(u-1)+2*(y-1)+v
                                   g(u,y,v) corresponds to entry 4*(u-1)+2*(y-1)+v
                        [g(1,1,1)=1, g(1,1,2)=2, g(1,2,1)=3, g(1,2,2)=4,
                                   g(2,1,1)=5, g(2,1,2)=6, g(2,2,1)=7, g(2,2,2)=8]
                                          y=1:[1,2,5,6], y=2:[3,4,7,8]   

ALTERNATIVE TRANSCRIPTION: 

output[2][aggregate_state][point_on_the_grid]
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
function period_T(T::Int64, start::Int64, grid_ave::Array{Float64, 1}, fut_a_spl::Array{Array{Spline1D, 1}, 1}, K_prev_spl::Array{Spline1D,1}, sol_prev_spl::Array{Array{Spline1D, 1}, 1}, thresh_mc::Float64, thresh_d::Float64, solve_fac::Float64, solve_ftol::Float64, α::Float64, β::Float64, δ::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1}, ISpd::Vector{Vector{Float64}}, atpm::Matrix{Float64}, itpm::Array{Float64, 4}, kill_switch::Int64, kill_switch2::Int64, kill_switch3::Int64)
    local nn::Int64, len_grid::Int64, first_success::Bool, alarm_off::Bool, loc_result, result, result1, result2, result_prev, result_prev_prev, test_K::Bool;
    test_K=false;
    len_grid=length(grid_ave);
    nn = start;
    println(nn);
    alarm_off = false;
    # first pass
    loc_result = [@spawn solve_loc(nn, 1, arg, [arg, arg], fut_a_spl, K_prev_spl, thresh_d, thresh_mc, sol_prev_spl, solve_fac, solve_ftol, α, β, δ, NN, XX, AS, IS, ISpd, atpm, itpm, kill_switch, kill_switch2, kill_switch3) for arg in grid_ave];
    #
    result1 = [fetch(loc_result[iii]) for iii=1:len_grid];
    #
    for iii=1:len_grid
        if (abs(result1[iii][4][3]-clear_mkt(nn,1,result1[iii][4][4][2][9:10],grid_ave[iii],ISpd))>0.001)
            test_K=true
        end
    end
    if (test_K)
        println("state 1: wrong capital at n = ", nn)
    end
    test_K=false;
    #
    if (length(findall(x->(x==false) , result1)) > 0)
        first_success = false;
        println("alarm 1 at period T = ", nn)
    else
        first_success = true;
        loc_result = [@spawn solve_loc(nn, 2, arg, [arg, arg], fut_a_spl, K_prev_spl, thresh_d, thresh_mc, sol_prev_spl, solve_fac, solve_ftol, α, β, δ, NN, XX,  AS, IS, ISpd, atpm, itpm, kill_switch, kill_switch2, kill_switch3) for arg in grid_ave];
        result2 = [fetch(loc_result[iii]) for iii=1:len_grid];
        #
        for iii=1:len_grid
            if (abs(result2[iii][4][3]-clear_mkt(nn,2,result2[iii][4][4][2][9:10],grid_ave[iii],ISpd))>0.001)
                test_K=true
            end
        end
        if (test_K)
            println("state 2: wrong capital at n = ", nn)
        end
        test_K=false;
        #
        if length(findall(x->(x==false) , result2)) > 0
            first_success = false;
            println("alarm 2 at period T = ", nn)
        else
            result=[result1,result2];
            fut_a_spl = [[Spline1D(grid_ave, [result[x][iii][4][4][2][8+u] for iii=1:len_grid], k=3, bc = "extrapolate") for u in IS] for x in AS];
            K_prev_spl = [Spline1D(grid_ave, [result[x][iii][4][3] for iii=1:len_grid], k=3, bc = "extrapolate") for x in AS];
            sol_prev_spl = [[Spline1D(grid_ave, [result[x][iii][4][4][2][inx] for iii=1:len_grid], k=3, bc = "extrapolate") for inx=1:length(result[x][1][4][4][2])] for x in AS];
            alarm_off = true;
        end
    end
    # start iterations
    while (alarm_off & (nn<T))
        nn+=1;
        println(nn);
        if (nn > start+1)
            result_prev_prev=result_prev;
        end
        result_prev = result; # available only if n ≥ start+1
        loc_result = [@spawn solve_loc(nn, 1, grid_ave[iii], result_prev[1][iii][3], fut_a_spl, K_prev_spl, thresh_d, thresh_mc, sol_prev_spl, solve_fac, solve_ftol, α, β, δ, NN, XX,  AS, IS, ISpd, atpm, itpm, kill_switch, kill_switch2, kill_switch3) for iii=1:len_grid];
        #
        result1 = [fetch(loc_result[iii]) for iii=1:len_grid];
        #
        for iii=1:len_grid
            if (abs(result1[iii][4][3]-clear_mkt(nn,1,result1[iii][4][4][2][9:10],grid_ave[iii],ISpd))>0.001)
                test_K=true
            end
        end
        if (test_K)
            println("state 1: wrong capital at n = ", nn)
        end
        test_K=false;
        #
        if (length(findall(x->(x==false) , result1)) > 0)
            alarm_off = false;
            println("alarm 1 at period T = ", nn)
        else
            alarm_off = true;
            loc_result = [@spawn solve_loc(nn, 2, grid_ave[iii], result_prev[1][iii][3], fut_a_spl, K_prev_spl, thresh_d, thresh_mc, sol_prev_spl, solve_fac, solve_ftol, α, β, δ, NN, XX,  AS, IS, ISpd, atpm, itpm, kill_switch, kill_switch2, kill_switch3) for iii=1:len_grid];
            result2 = [fetch(loc_result[iii]) for iii=1:len_grid];
            #
            for iii=1:len_grid
                if (abs(result2[iii][4][3]-clear_mkt(nn,2,result2[iii][4][4][2][9:10],grid_ave[iii],ISpd))>0.001)
                    test_K=true
                end
            end
            if (test_K)
                println("state 2: wrong capital at n = ", nn)
            end
            test_K=false;
            #
            if length(findall(x->(x==false) , result2)) > 0
                alarm_off = false;
                println("alarm 2 at period T = ", nn)
            else
                result=[result1,result2];
                fut_a_spl = [[Spline1D(grid_ave, [result[x][iii][4][4][2][8+u] for  iii=1:len_grid], k=3, bc = "extrapolate") for u in IS] for x in AS];
                K_prev_spl = [Spline1D(grid_ave, [result[x][iii][4][3] for iii=1:len_grid], k=3, bc = "extrapolate") for x in AS];
                sol_prev_spl = [[Spline1D(grid_ave, [result[x][iii][4][4][2][inx] for iii=1:len_grid], k=3, bc = "extrapolate") for inx=1:length(result[x][1][4][4][2])] for x in AS];
                #cs_prev_spl = [Spline2D(xs, ys, [result[x][i,j][4][4][1] for i=1:lx, j=1:ly]) for x in AS];
            end
        end
        #
    end
    if (nn > start)
        if alarm_off
            return (nn, result, result_prev)
        elseif (nn > start+1)
            return (nn-1, result_prev, result_prev_prev)
        else
            return (nn-1, result_prev)
        end
    else
        if first_success
            return (nn, result)
        else
            println("method failed")
            return false
        end
    end
end

#############

begin
    using Random
    rng=MersenneTwister(42)
    #rng=RandomDevice()
end;


function gen_state(sno::Int64, atpm::Array{Float64, 2})
    local uuu::Float64,rslt::Int64;
    uuu=rand(rng)
    if uuu<atpm[sno,1] rslt=1::Int64 else rslt=2::Int64 end
    return rslt
end

function gen_fut_dist(ini_state::Int64, ini_dist::Array{Float64, 1}, iter::Int64, α::Float64, β::Float64, δ::Float64, IS::Vector{Int64}, ISpd::Vector{Vector{Float64}}, itpm::Array{Float64, 4}, atpm::Array{Float64, 2}, LL::Vector{Float64}, XX::Vector{Float64}, K_spl::Vector{Spline1D}, sol_spl::Vector{Vector{Spline1D}})
    local now_state::Int64, now_dist::Array{Float64, 1}, i_no::Int64;
    AR=Array{Int64,1}(undef,iter);
    A=Array{Array{Float64,1},1}(undef,iter);
    now_state=ini_state;
    now_dist=ini_dist;
    AR[1]=now_state;
    A[1]=now_dist;
    i_no=1;
    while i_no<iter
        local fut_state::Int64, fut_dist::Array{Float64, 1};
        i_no+=1;
        fut_state=gen_state(now_state,atpm)
        fut_dist=transport_mult(now_state, fut_state, α, β, δ, IS, ISpd, itpm, LL, XX, K_spl, sol_spl, now_dist);
        AR[i_no]=fut_state;
        A[i_no]=fut_dist;
        now_state=fut_state;
        now_dist=fut_dist;
    end
    return AR, A
end


function transport_a(xxx::Int64, yyy::Int64, α::Float64, β::Float64, δ::Float64, IS::Vector{Int64}, ISpd::Vector{Vector{Float64}}, itpm::Array{Float64, 4}, LL::Vector{Float64}, XX::Vector{Float64}, K_spl::Vector{Spline1D}, sol_spl::Vector{Vector{Spline1D}}, current_ave::Float64)
    return β*(ρ(yyy,α,LL,XX)(K_spl[xxx](current_ave))+1-δ)*current_ave + sum([ISpd[xxx][u]*itpm[u,v,xxx,yyy]*sol_spl[xxx][4*(u-1)+2*(yyy-1)+v](current_ave) for u in IS, v in IS]);
end


function transport_mult(xxx::Int64, yyy::Int64, α::Float64, β::Float64, δ::Float64, IS::Vector{Int64}, ISpd::Vector{Vector{Float64}}, itpm::Array{Float64, 4}, LL::Vector{Float64}, XX::Vector{Float64}, K_spl::Vector{Spline1D}, sol_spl::Vector{Vector{Spline1D}}, current_avgs::Vector{Float64})
    local ave::Float64
    ave=ISpd[xxx]'*current_avgs
    return [sum([(ISpd[xxx][u]*itpm[u,v,xxx,yyy]/ISpd[yyy][v])*(β*(ρ(yyy,α,LL,XX)(K_spl[xxx](ave))+1-δ)*current_avgs[u] + sol_spl[xxx][4*(u-1)+2*(yyy-1)+v](ave)) for u in IS]) for v in AS];
end


function kernel_err(xxx::Int64, uuu::Int64, α::Float64, β::Float64, δ::Float64, AS::Vector{Int64}, IS::Vector{Int64}, atpm::Matrix{Float64}, itpm::Array{Float64, 4}, LL::Vector{Float64}, XX::Vector{Float64}, K_spl::Vector{Spline1D}, sol_spl::Vector{Vector{Spline1D}}, current_ave::Float64, ccc::Float64)
    return 1-sum([atpm[xxx,y]*itpm[uuu,v,xxx,y]*(β*(ρ(y,α,LL,XX)(K_spl[xxx](current_ave))+1-δ)/((1/ccc)*(sol_spl[xxx][4*(uuu-1)+2*(y-1)+v](current_ave))+β*(ρ(y,α,LL,XX)(K_spl[xxx](current_ave))+1-δ))) for y in AS, v in IS]);
end

