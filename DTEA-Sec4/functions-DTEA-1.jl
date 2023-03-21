######################################################################################################
#
# Julia code
#
# Functions to implement the method described in Sec. 4 in the paper
#       "Dynamic transportation of economic agents" [DTEA] by Andrew Lyasoff
#
# The code provides an alternative solution to the example from the paper
#       Krusell, Per, and Anthony Smith. (1998). Income and wealth heterogeneity in the macroeconomy.
#                         Journal of Political Economy 106 867-896.
#                
# This code supplements the paper "Dynamic transportation of economic agents" [DTEA]
#                                        by Andrew Lyasoff (www.andrewlyasoff.tech)
#
# Copyright © 2019-2023 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################

#=
NB: in [DTEA] present and future employments are u and v; in this code those are s and σ
    in [DTEA] present and future productivity states are x and y; in this code those are x and ξ
=#


function FGQT(f,a,b)
    return (arg->(1/2)*(b-a)*f( (1/2)*(b-a)*arg + (1/2)*(a+b) ))
end


## get the steady-state probabilities for a given transition matrix
function STSP(trans_prob::Array{Float64,2})
    local nos::Int64, mM::Array{Float64,2}, rhs::Array{Float64,1};
    nos=size(trans_prob)[1];
    if maximum(abs.([sum(trans_prob[i,:]) for i=1:nos]-ones(nos)))>10.0e-16
        println("not a TPM")
    end
    mM=[(trans_prob-I)';ones(nos)'];
    rhs=push!(zeros(nos),1.0);
    return (mM'*mM)\(mM'*rhs)
end

#=
return on capital as a function of installed capital (arg_x)

XX is the list of aggregate shocks (in various aggregate states)
NN is the list of aggregate supplied labor
α is the risk aversion
x is the present state (superfluous in the next function)
ξ is the future state
=#
function ρ(x::Int64, ξ::Int64, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1})
    return (arg_x -> XX[ξ]*α*(arg_x/NN[ξ])^(α-1))
end

#
# wages in the future period as function of the capital installed in the present period
#
function w(x::Int64, ξ::Int64, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1})
    return ( arg_y -> XX[ξ]*(1-α)*(arg_y/NN[ξ])^α )
end

#=
NB: all slopes are in exact form; only the intercepts are unknown

UNKNOWNS (t and x are given):

g(s,ξ,σ) corresponds to entry 4*(s-1)+2*(ξ-1)+σ    s=1:2, ξ=1:2, σ=1:2 | belong to t+1

a(x,s) corresponds to entry 8+s  s=1:2 | belongs to t

NB: the intercepts of the future consumption lines are the first 8 unknowns out of 10
NB: the intercepts of the present portfolio lines are the last 2 unknowns out of 10 

GIVENS:

a(ξ,σ) , y=1:2, v=1:2 | belong to t+1

NB: the intercepts of the portfolio lines for the next period are given
and depend on the future productivit state and the future employment state

=#

# system of forst order Lagrange condition (FOLC) to solve equations to solve -- see (4.11) and (4.10)
# 10 equations with 10 unknowns
#
function sstm(n::Int64,  xxx::Int64, cs::Vector{Float64}, KK::Float64, a_prev::Array{Array{Float64, 1}, 1}, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1})
    return (y->(vcat(vcat(vcat([[[(a_prev[ξ][σ] + y[4*(s-1)+2*(ξ-1)+σ]*((1-β^n)/(1-β)) - y[8+s]*(ρ(xxx, ξ, α, NN, XX)(KK)+1-δ)-SS[σ]*w(xxx,  ξ, α, NN, XX)(KK)) for σ in IS] for ξ in AS] for s in IS]...)...) , [(sum([β*(1/(y[4*(s-1)+2*(ξ-1)+σ]/cs[s] +β*(ρ(xxx, ξ, α, NN, XX)(KK)+1-δ)))*(ρ(xxx, ξ, α, NN, XX)(KK)+1-δ)*atpm[xxx,ξ]*itpm[s,σ,xxx,ξ] for σ in IS, ξ in AS])-1.0) for s in IS ] )))
end

# method for extrapolating a 2D-spline
function outside2D(ff::Spline2D,  mm::Int64, no_nodes::Int64, args::Array{Float64, 1})
    local b1::Array{Float64, 1}, b2::Array{Float64, 1}, abscissas::Array{Float64, 1}, vlus::Array{Float64, 1}, loc_spline::Spline1D, check_package::Bool;
    b1 = get_knots(ff)[1];
    b2 = get_knots(ff)[2];
    if (mm == 1) # the out-of-bounds argument
        if (args[2]>=b2[1])&(args[2]<=b2[end]) # arg2 is in bounds
            if (args[1]>=b1[1])&(args[1]<=b1[end]) #arg1 is inbounds
                return ff(Float64(args[1]),Float64(args[2]))
            else #arg1 is not in bounds 
                abscissas = [(b1[1]+i*(b1[end]-b1[1])/no_nodes) for i=0:no_nodes];
                vlus = [ff(Float64(abscissas[i]),Float64(args[2]))  for i=1:(1+no_nodes)]; # interpolate over arg1
                loc_spline = Spline1D(abscissas, vlus, k=3, bc = "extrapolate"); 
                return loc_spline(Float64(args[1]))
            end
        else #arg2 is not in boounds as stipulated by mm
            println("Arg 2 is not in bounds.")
            return NaN
        end
    elseif (mm == 2)
        if (args[1]>=b1[1])&(args[1]<=b1[end]) # arg1 is in bounds
            if (args[2]>=b2[1])&(args[2]<=b2[end]) # arg2 is in bounds
                return ff(Float64(args[1]),Float64(args[2]))
            else # arg2 is not in bounds
                abscissas = [(b2[1]+i*(b2[end]-b2[1])/no_nodes) for i=0:no_nodes];
                vlus = [ff(Float64(args[1]),Float64(abscissas[i]))  for i=1:(1+no_nodes)]; # interpolate over arg2
                loc_spline = Spline1D(abscissas, vlus, k=3, bc = "extrapolate");
                return loc_spline(Float64(args[2]))
            end
        else
            println("Arg 1 is not in bounds.")
            return NaN
        end
    else
        println("Argument label must be 1 or 2.")
        return NaN
    end
end

# second method for extrapolating a 2D-spline
function XSpline2D(ff::Spline2D, mm::Int64, no_nodes::Int64, args::Array{Float64, 1})
    local b1::Array{Float64, 1}, b2::Array{Float64, 1}, abscissas::Array{Float64, 1}, vlus::Array{Float64, 1}, loc_spline::Spline1D, check_package::Bool;
    b1 = get_knots(ff)[1];
    b2 = get_knots(ff)[2];
    if (mm == 1) # interpolation must be over arg1
        if (args[2]>=b2[1])&(args[2]<=b2[end]) # arg2 IS in bounds
            if (args[1]>=b1[1])&(args[1]<=b1[end]) #arg1 IS in bounds
                return ff(Float64(args[1]),Float64(args[2]))
            else # arg1 IS NOT in bounds
                return outside2D(ff, 1, no_nodes, args)
            end 
        else # arg2 IS NOT in bounds
            if (args[1]>=b1[1])&(args[1]<=b1[end]) #arg1 IS in bounds
                return outside2D(ff, 2, no_nodes, args)
            else #arg1 IS NOT in bounds
                abscissas = [(b1[1]+i*(b1[end]-b1[1])/no_nodes) for i=0:no_nodes];
                vlus = [outside2D(ff, 2, no_nodes, [abscissas[i],args[2]]) for i=1:(1+no_nodes)]; # interpolate over arg1
                loc_spline = Spline1D(abscissas, vlus, k=3, bc = "extrapolate"); 
                return loc_spline(Float64(args[1]))
            end
        end
    else (mm == 2) # interpolation must be over arg2
        if (args[1]>=b1[1])&(args[1]<=b1[end]) # arg1 IS in bounds
            if (args[2]>=b2[1])&(args[2]<=b2[end]) #arg2 IS in bounds
                return ff(Float64(args[1]),Float64(args[2]))
            else # arg2 IS NOT in bounds
                return outside2D(ff, 2, no_nodes, args)
            end 
        else # arg1 IS NOT in bounds
            if (args[2]>=b2[1])&(args[2]<=b2[end]) #arg2 IS in bounds
                return outside2D(ff, 1, no_nodes, args)
                else #arg2 IS NOT in bounds
                abscissas = [(b2[1]+i*(b2[end]-b2[1])/no_nodes) for i=0:no_nodes];
                vlus = [outside2D(ff, 1, no_nodes, [args[1],abscissas[i]]) for i=1:(1+no_nodes)]; # interpolate over arg2
                loc_spline = Spline1D(abscissas, vlus, k=3, bc = "extrapolate"); 
                return loc_spline(Float64(args[2]))
            end
        end
    end
end

# the left side of the market clearing condition -- see (4.6) and (4.8)
function clear_mkt(n::Int64, xxx::Int64, sol_sstm1a::Array{Float64, 1}, AA::Array{Float64, 1},  AS::Array{Int64,1}, IS::Array{Int64,1}, ISpd::Array{Array{Float64,1},1})
    return sum([ISpd[xxx][s]*(sol_sstm1a[s] + sum([β^i for i=1:n])*AA[s]) for s in IS])
end

#=
solves the system composed of (4.10) and (4.11) for the 10 unknown intercepts (8 for the future consumption lines
                                              and 2 for the present portfolio lines)
=#
function get_cs(n::Int64, xxx::Int64, cs0::Vector{Float64}, thresh::Float64, solve_ini_guess1::Array{Float64, 1}, solve_fac::Float64, solve_ftol::Float64, KK::Float64, a_prev::Array{Array{Float64, 1}, 1}, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1}, kill_switch::Int64)
    local css::Vector{Float64}, lcl_sln, test::Bool;
    css = cs0;
    test = false;
    lcl_sln=nlsolve(sstm(n, xxx, css, KK, a_prev, α, NN, XX, AS, IS), solve_ini_guess1, factor=solve_fac, ftol=solve_ftol);
    test=converged(lcl_sln);
    if test
        return (css,1,lcl_sln.zero)
    else
        return false
    end
end

#=
tâtonnement over capital:

at a given (fixed) point on the ℝ²-grid performs steps (2) through (6) in the generic backward step,
or steps (2) through (5) in the initial backward step  from 4.3
=# 
function find_K02(n::Int64, xxx::Int64, dist::Array{Float64, 1}, K_ini::Float64, cs0::Array{Float64, 1}, thresh_lin::Float64, thresh_mc::Float64, solve_ini_guess1::Array{Float64, 1}, solve_fac::Float64, solve_ftol::Float64,  a_prev::Array{Array{Float64, 1}, 1}, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1}, kill_switch::Int64, kill_switch2::Int64)
    local K_prev::Float64, K_last::Float64, mc_test::Float64, all_go::Bool, sol, iter::Int64;
    iter = 0;
    mc_test = Inf;
    K_last = K_ini;
    sol = get_cs(n, xxx, cs0, thresh_lin, solve_ini_guess1, solve_fac, solve_ftol, K_last, a_prev, α, NN, XX, AS, IS, kill_switch);
    #println(sol[2])
    if length(sol)>1 
        K_last = clear_mkt(n, xxx, sol[3][9:10], dist, AS, IS, ISpd);
        if (K_last<1.0e-4) K_last=0.2 end;
        all_go = true;
        iter = 1;
    else
        all_go = false;
    end
    while (all_go&(iter<kill_switch2)&(mc_test>thresh_mc))
        #println(iter)
        #sol=get_cs(n, xxx, 0.35*sol[1], thresh_lin, solve_ini_guess1, solve_fac, solve_ftol, K_last, a_prev, α, NN, XX, AS, IS, kill_switch);
        sol=get_cs(n, xxx, cs0, thresh_lin, solve_ini_guess1, solve_fac, solve_ftol, K_last, a_prev, α, NN, XX, AS, IS, kill_switch);
        #println(sol[2])
        if length(sol)>1
            K_prev=K_last;
            K_last=clear_mkt(n, xxx, sol[3][9:10], dist, AS, IS, ISpd);
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
tâtonnement over future distribution (an element of ℝ²):

at a given (fixed) point on the ℝ²-grid performs steps (1) through (8) in the generic backward step from 4.3
=#                
function solve_loc(n::Int64, xxx::Int64, loc_dist::Array{Float64, 1}, ansatz_fut_dist::Array{Array{Float64, 1}, 1}, interp::Array{Array{Spline2D, 1}, 1}, K_prev_sp::Array{Spline2D, 1}, cs0::Vector{Float64}, thresh_d::Float64, thresh_lin::Float64, thresh_mc::Float64, sol_prev_sp::Array{Array{Spline2D, 1}, 1}, solve_fac::Float64, solve_ftol::Float64, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1}, kill_switch::Int64, kill_switch2::Int64, kill_switch3::Int64)
    local  loc_sol, next_dist_i::Array{Array{Float64, 1}, 1}, distance_d::Float64, iter_no::Int64, K_ini::Float64, solve_ini_guess1::Array{Float64, 1}, no_alaram::Bool;
    fut_a = [[XSpline2D(interp[ξ][σ], 1, 50, [ansatz_fut_dist[ξ][1] , ansatz_fut_dist[ξ][2]/ansatz_fut_dist[ξ][1]]) for σ in IS] for ξ in AS];
    #K_ini = K_prev_sp[xxx](loc_dist[1],loc_dist[2]/loc_dist[1]);
    K_ini = XSpline2D(K_prev_sp[xxx], 1, 50, [loc_dist[1] , loc_dist[2]/loc_dist[1]]);
    solve_ini_guess1 = [XSpline2D(sol_prev_sp[xxx][inx], 1, 50, [loc_dist[1] , loc_dist[2]/loc_dist[1]]) for inx=1:length(sol_prev_sp[xxx])];
    loc_sol = find_K02(n, xxx, loc_dist, K_ini, cs0, thresh_lin, thresh_mc, solve_ini_guess1, solve_fac, solve_ftol, fut_a, α, NN, XX,  AS, IS, kill_switch, kill_switch2)
    if length(loc_sol) > 1
        no_alarm = true;
        next_dist_i = [[sum([(ISpd[xxx][s]*itpm[s,σ,xxx,ξ]/ISpd[ξ][σ])*(loc_sol[4][3][4*(s-1)+2*(ξ-1)+σ] + β*(ρ(xxx, ξ, α, NN, XX)(loc_sol[3])+1-δ)*loc_dist[s]) for s in IS]) for σ in IS] for ξ in AS];
        dist_d = maximum([maximum(abs.(next_dist_i[ξ].-ansatz_fut_dist[ξ])) for ξ in AS]);
        iter_no = 1;
    else
        no_alarm = false;
        dist_d = -Inf;
        iter_no = 2^50;
    end
    while no_alarm&(dist_d > thresh_d)&(iter_no < kill_switch3)
        ansatz_fut_dist = next_dist_i;
        fut_a = [[XSpline2D(interp[ξ][σ], 1, 50, [ansatz_fut_dist[ξ][1] , ansatz_fut_dist[ξ][2]/ansatz_fut_dist[ξ][1]]) for σ in IS] for ξ in AS];
        loc_sol = find_K02(n, xxx, loc_dist, K_ini, cs0, thresh_lin, thresh_mc, solve_ini_guess1, solve_fac, solve_ftol, fut_a, α, NN, XX,  AS, IS, kill_switch, kill_switch2);
        if length(loc_sol) > 1
            no_alarm = true;
            next_dist_i = [[sum([(ISpd[xxx][s]*itpm[s,σ,xxx,ξ]/ISpd[ξ][σ])*(loc_sol[4][3][4*(s-1)+2*(ξ-1)+σ] + β*(ρ(xxx, ξ, α, NN, XX)(loc_sol[3])+1-δ)*loc_dist[s]) for s in IS]) for σ in IS] for ξ in AS];
            iter_no+=1
            dist_d=maximum([maximum(abs.(next_dist_i[ξ].-ansatz_fut_dist[ξ])) for ξ in AS]);
        else
            no_alarm = false;
        end;
    end;
    if ((dist_d > thresh_d)&(iter_no == kill_switch3))
        no_alarm=false;
        println("solve_loc reached the iterations maximum of ", kill_switch3);
    end;
    if no_alarm
        return (iter_no, dist_d, next_dist_i, loc_sol)
    else
        return false
    end
end

#=
completes the initial backward step from 4.3 in [DTEA] for all points on the ℝ²-grid with parallelization:
OUTPUTS: the values for the 10 intercepts and the aggregate capital at ALL points on the grid
and for all productivity states; returns also the interpolated versions of those values (this redundacy is necssary
      because splines cannot be dumped)
=#
function period_Tm1(xs::Array{Float64, 1}, ys::Array{Float64, 1}, K_ini::Float64, cs0::Float64, solve_ini_guess1::Array{Float64, 1}, thresh_lin::Float64, thresh_mc::Float64, solve_fac::Float64, solve_ftol::Float64, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1}, kill_switch::Int64, kill_switch2::Int64)
    local fut_a::Array{Array{Float64, 1}, 1}, lx::Int64, ly::Int64, alarm_off::Bool, loc_result, result1, result2, result_prev, result, interp::Array{Array{Spline2D, 1}, 1}, K_prev_spl::Array{Spline2D, 1}, cs_prev_spl::Array{Spline2D, 1}, sol_prev_spl::Array{Array{Spline2D, 1}, 1};
    lx=length(xs);
    ly=length(ys);
    fut_a=[[0.0 for σ in IS] for ξ in AS];
    loc_result = [@spawn find_K02(1, 1, [x,x*y], K_ini, [x,x*y], thresh_lin, thresh_mc, solve_ini_guess1, solve_fac, solve_ftol, fut_a, α, NN, XX, AS, IS, kill_switch, kill_switch2) for x in xs, y in ys];
    #
    result1 = [fetch(loc_result[i,j]) for i=1:lx, j=1:ly];
    #
    if (length(findall(x->(x==false) , result1)) > 0)
        alarm_off = false;
        println("alarm 1 at period T = ", 1)
    else
        alarm_off = true;
        loc_result = [@spawn find_K02(1, 2, [x,x*y], K_ini, [x,x*y], thresh_lin, thresh_mc, solve_ini_guess1, solve_fac, solve_ftol, fut_a, α, NN, XX, AS, IS, kill_switch, kill_switch2) for x in xs, y in ys];
        result2 = [fetch(loc_result[i,j]) for i=1:lx, j=1:ly];
        if length(findall(x->(x==false) , result2)) > 0
            alarm_off = false;
            println("alarm 2 at period T = ", 1)
            return false
        else
            result=[result1,result2];
            interp = [[Spline2D(xs, ys, [result[x][i,j][4][3][8+s] for i=1:lx, j=1:ly]) for s in IS] for x in AS];
            K_prev_spl = [Spline2D(xs, ys, [result[x][i,j][3] for i=1:lx, j=1:ly]) for x in AS];
            sol_prev_spl = [[Spline2D(xs, ys, [result[x][i,j][4][3][inx] for i=1:lx, j=1:ly]) for inx=1:length(result[x][1,1][4][3])] for x in AS];
            #cs_prev_spl = [Spline2D(xs, ys, [result[x][i,j][4][1] for i=1:lx, j=1:ly]) for x in AS];
            return (result, interp, K_prev_spl, sol_prev_spl)
        end
    end
end

#=
repeats the generic backward step from 4.3 in [DTEA] for all points on the ℝ²-grid with parallelization,
                                                until the final backward step is reached:
OUTPUTS: the values for the 10 intercepts and the aggregate capital at ALL points on the grid
and for all productivity states; returns also the interpolated versions of those values (this redundacy is necssary
because splines cannot be dumped)

TO BE PRECISE

out[1] = total number of completed iterations
out[2] = results from the last iteration
out[3] = results from next to last iterations
out[2][x][i,j] = solution from the last iteration at aggregate state x and grid point [i,j]
out[2][x][i,j][1] = number of iterations between steps 3-7  ::Int64
out[2][x][i,j][2] = future distribution mismatch  ::Float64
out[2][x][i,j][3] = future distribution in high and in low states ::Vector{Vector{Float64}}
out[2][x][i,j][4=end] = ::Tuple{Int64, Float64, Float64, Tuple{Vector{Float64}, Int64, Vector{Float64}}}
out[2][x][i,j][4][1] = number of iterations between steps 3-6
out[2][x][i,j][4][2] = mismatch in the future distribution
out[2][x][i,j][4][3] = installed capital
out[2][x][i,j][4][4] = actual solution (list of 10 floats)
         out[2][x][i,j][4][4][9:10] = portfolio lines intercepts in high and low states
         out[2][x][i,j][4][4][1:8]  = furure consumption lines intercepts 4*(u-1)+2*(y-1)+v
                                   g(s,ξ,σ) corresponds to entry 4*(s-1)+2*(ξ-1)+σ
                        [g(1,1,1)=1, g(1,1,2)=2, g(1,2,1)=3, g(1,2,2)=4,
                                   g(2,1,1)=5, g(2,1,2)=6, g(2,2,1)=7, g(2,2,2)=8]
                                          y=1:[1,2,5,6], y=2:[3,4,7,8]   

=#
function period_T(T::Int64, start::Int64, xs::Array{Float64, 1}, ys::Array{Float64, 1}, interp::Array{Array{Spline2D, 1}, 1}, K_prev_spl::Array{Spline2D, 1}, cs0::Float64, sol_prev_spl::Array{Array{Spline2D, 1}, 1}, thresh_lin::Float64, thresh_mc::Float64, thresh_d::Float64, solve_fac::Float64, solve_ftol::Float64, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1}, kill_switch::Int64, kill_switch2::Int64, kill_switch3::Int64)
    local nn::Int64, lx::Int64, ly::Int64, first_success::Bool, alarm_off::Bool, loc_result, result, result1, result2, result_prev, result_prev_prev;
    lx=length(xs);
    ly=length(ys);
    nn = start;
    println(nn);
    alarm_off = false;
    # first pass
    loc_result = [@spawn solve_loc(nn, 1, [x,x*y], [[x,x*y],[x,x*y]], interp, K_prev_spl, [x,x*y], thresh_d, thresh_lin, thresh_mc, sol_prev_spl, solve_fac, solve_ftol, α, NN, XX,  AS, IS, kill_switch, kill_switch2, kill_switch3) for x in xs, y in ys];
    #
    result1 = [fetch(loc_result[i,j]) for i=1:lx, j=1:ly];
    #
    if (length(findall(x->(x==false) , result1)) > 0)
        first_success = false;
        println("alarm 1 at period T = ", nn)
    else
        first_success = true;
        loc_result = [@spawn solve_loc(nn, 2, [x,x*y], [[x,x*y],[x,x*y]], interp, K_prev_spl, [x,x*y], thresh_d, thresh_lin, thresh_mc, sol_prev_spl, solve_fac, solve_ftol, α, NN, XX,  AS, IS, kill_switch, kill_switch2, kill_switch3) for x in xs, y in ys];
        result2 = [fetch(loc_result[i,j]) for i=1:lx, j=1:ly];
        if length(findall(x->(x==false) , result2)) > 0
            first_success = false;
            println("alarm 2 at period T = ", nn)
        else
            result=[result1,result2];
            interp = [[Spline2D(xs, ys, [result[x][i,j][4][4][3][8+s] for i=1:lx, j=1:ly]) for s in IS] for x in AS];
            K_prev_spl = [Spline2D(xs, ys, [result[x][i,j][4][3] for i=1:lx, j=1:ly]) for x in AS];
            sol_prev_spl = [[Spline2D(xs, ys, [result[x][i,j][4][4][3][inx] for i=1:lx, j=1:ly]) for inx=1:length(result[x][1,1][4][4][3])] for x in AS];
            #cs_prev_spl = [Spline2D(xs, ys, [result[x][i,j][4][4][1] for i=1:lx, j=1:ly]) for x in AS];
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
        loc_result = [@spawn solve_loc(nn, 1, [xs[i],xs[i]*ys[j]], result_prev[1][i,j][3], interp, K_prev_spl, [xs[i],xs[i]*ys[j]], thresh_d, thresh_lin, thresh_mc, sol_prev_spl, solve_fac, solve_ftol, α, NN, XX,  AS, IS, kill_switch, kill_switch2, kill_switch3) for i=1:lx, j=1:ly];
        #
        result1 = [fetch(loc_result[i,j]) for i=1:lx, j=1:ly];
        #
        if (length(findall(x->(x==false) , result1)) > 0)
            alarm_off = false;
            println("alarm 1 at period T = ", nn)
        else
            alarm_off = true;
            loc_result = [@spawn solve_loc(nn, 2, [xs[i],xs[i]*ys[j]], result_prev[2][i,j][3], interp, K_prev_spl, [xs[i],xs[i]*ys[j]], thresh_d, thresh_lin, thresh_mc, sol_prev_spl, solve_fac, solve_ftol, α, NN, XX,  AS, IS, kill_switch, kill_switch2, kill_switch3) for i=1:lx, j=1:ly];
            result2 = [fetch(loc_result[i,j]) for i=1:lx, j=1:ly];
            if length(findall(x->(x==false) , result2)) > 0
                alarm_off = false;
                println("alarm 2 at period T = ", nn)
            else
                result=[result1,result2];
                interp = [[Spline2D(xs, ys, [result[x][i,j][4][4][3][8+s] for i=1:lx, j=1:ly]) for s in IS] for x in AS];
                K_prev_spl = [Spline2D(xs, ys, [result[x][i,j][4][3] for i=1:lx, j=1:ly]) for x in AS];
                sol_prev_spl = [[Spline2D(xs, ys, [result[x][i,j][4][4][3][inx] for i=1:lx, j=1:ly]) for inx=1:length(result[x][1,1][4][4][3])] for x in AS];
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
    if sno==1
        if uuu<atpm[sno,1] rslt=1::Int64 else rslt=2::Int64 end
    else
        if uuu<atpm[sno,1] rslt=1::Int64 else rslt=2::Int64 end
    end
    return rslt
end

function gen_fut_dist(ini_state::Int64, ini_dist::Array{Float64, 1}, iter::Int64, ii::Int64, inn::Int64, dist_spline_1_11::Spline2D, dist_spline_1_12::Spline2D, dist_spline_1_21::Spline2D, dist_spline_1_22::Spline2D, dist_spline_2_11::Spline2D, dist_spline_2_12::Spline2D, dist_spline_2_21::Spline2D, dist_spline_2_22::Spline2D, atpm::Array{Float64, 2})
    local now_state::Int64, now_dist::Array{Float64, 1}, i_no::Int64;
    AR=Array{Int64,1}(undef,iter);
    A=Array{Array{Float64,1},1}(undef,iter);
    now_state=ini_state;
    now_dist=ini_dist;
    now_dist_mod=[now_dist[1],now_dist[2]/now_dist[1]];
    AR[1]=now_state;
    A[1]=now_dist;
    i_no=1;
    while i_no<iter
        local fut_state::Int64, fut_dist::Array{Float64, 1};
        i_no+=1;
        fut_state=gen_state(now_state,atpm)
        if now_state==1
            if fut_state==1
                fut_dist=[XSpline2D(dist_spline_1_11, ii, inn, now_dist_mod),XSpline2D(dist_spline_1_12, ii, inn, now_dist_mod)];
            else
                fut_dist=[XSpline2D(dist_spline_1_21, ii, inn, now_dist_mod),XSpline2D(dist_spline_1_22, ii, inn, now_dist_mod)];
            end
        else
            if fut_state==1
                fut_dist=[XSpline2D(dist_spline_2_11, ii, inn, now_dist_mod),XSpline2D(dist_spline_2_12, ii, inn, now_dist_mod)];
            else
                fut_dist=[XSpline2D(dist_spline_2_21, ii, inn, now_dist_mod),XSpline2D(dist_spline_2_22, ii, inn, now_dist_mod)];
            end
        end
        AR[i_no]=fut_state;
        A[i_no]=fut_dist;
        now_state=fut_state;
        now_dist=fut_dist;
        now_dist_mod=[now_dist[1],now_dist[2]/now_dist[1]];
    end
    return AR, A
end


function fut_cdf(KK::Float64, cdff, nus::Array{Float64, 1}, xxx::Int64, ξ::Int64, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1})
    return [( u -> sum([(ISpd[xxx][s]*itpm[s,σ,xxx,ξ]/ISpd[ξ][σ])*cdff[s]( (1/(β*(ρ(xxx, ξ, α, NN, XX)(KK)+1-δ)))*( u - nus[4*(s-1)+2*(ξ-1)+σ]) ) for s in IS]) )    for σ in IS]
end


#=
function fut_cdf(KK::Float64, cdff, nus::Array{Float64, 3}, xxx::Int64, ξ::Int64, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1})
    return [( u -> sum([(ISpd[xxx][s]*itpm[s,σ,xxx,ξ]/ISpd[ξ][σ])*cdff[s]( (1/(β*(ρ(xxx, ξ, α, NN, XX)(KK)+1-δ)))*( u - nus[σ,ξ,s]) ) for s in IS]) )    for σ in IS]
end
=#

# expected value for a given CDF
# should be used for control purposes only
function Mcdf(F, steps::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}})
    return 0.5*dot( (F.(steps[2:end])-F.(steps[1:end-1])) , (steps[2:end] + steps[1:end-1]) )
end

        
function forward_dist(nnn::Int64, AR::Array{Int64, 1}, FU_ini::Array{Any, 1}, ub::Float64, ubs::Float64, numb_absc::Int64, bounds_tol1::Float64, bounds_tol2::Float64, iter_stop::Int64, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1})
    local xxx::Int64, f_a::Array{Float64,1}, K0::Float64, nunus::Array{Float64, 3}, ξξξ::Int64, FU::Array{Any, 1}, next_bounds::Array{Tuple{Float64, Float64}, 1}, spl1::Spline1D, spl2::Spline1D, absc1::Array{Any, 1}, absc2::Array{Any, 1}, vals1::Array{Any, 1}, vals2::Array{Any, 1};
    FU=FU_ini
    while nnn<iter_stop
        nnn+=1;
        xxx = AR[nnn];
        f_a=[Mcdf(FU[1],0.0:0.00001:ub),Mcdf(FU[2],0.0:0.00001:ub)];
        K0=K_prev_spl[xxx](f_a[1],f_a[2]/f_a[1]);
        nunus=[nu_spline[xxx][4*(s-1)+2*(ξ-1)+σ](f_a[1],f_a[2]/f_a[1]) for σ in IS, ξ in AS, s in IS];
        ξξξ=AR[nnn+1];
        FU=fut_cdf(K0, [FU[1],FU[2]], nunus, xxx, ξξξ, α, NN, XX);
        #println(nnn, " ", f_a)
        next_bounds = [(0.9*find_zero((arg->(FU[i](arg)-(bounds_tol1))), (0.0,f_a[i])), 1.1*find_zero((arg->(FU[i](arg)-(1-bounds_tol2))), (f_a[i],ubs))) for i=1:2];
        absc1=collect(next_bounds[1][1]:((next_bounds[1][2]-next_bounds[1][1])/numb_absc):next_bounds[1][2]);
        vals1=FU[1].(absc1);
        #
        pushfirst!(absc1,0.0);
        push!(absc1,next_bounds[1][2]+1.5);
        #
        pushfirst!(vals1,0.0);
        push!(vals1,1.0);
        ##
        absc2=collect(next_bounds[2][1]:((next_bounds[2][2]-next_bounds[2][1])/numb_absc):next_bounds[2][2]);
        vals2=FU[2].(absc2);
        #
        pushfirst!(absc2,0.0);
        push!(absc2,next_bounds[2][2]+1.5);
        #
        pushfirst!(vals2,0.0);
        push!(vals2,1.0);
        #
        spln1=Spline1D(absc1, vals1, k=1, bc = "nearest");
        spln2=Spline1D(absc2, vals2, k=1, bc = "nearest");
        FU=[spln1,spln2];
    end
    return nnn, FU
end
