######################################################################################################
#
# Julia code
#
# Illustrates the method described in Sec. 4 in the paper
#       "Another look at the distribution of income and wealth in the Macroeconomy" (DIW) by Andrew Lyasoff
#
# Provides an alternative solution to the example from the paper
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

#=
NB: in DIW present and future employments are u and v; in this code those are s and σ
    in DIW present and future productivity states are x and y; in this code those are x and ξ
=#

# not needed if run-DIW-1-mod.jl
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
    include("ini-setup-DIW-1-mod.jl"); # include("ini-setup-DIW-1-mod.jl") to change β from original
    gqx, gqw = gausslegendre( 100 );
    GQX, GQW = gausslegendre( 10_000 );
    IS=[1,2]; # idiosyncratic states
    AS=[1,2]; # productivity states
    RRR=1.0;  # risk aversion
    ISpd=[[0.96, 0.04], [0.9, 0.1]] #same as PJL' (population distribution over employment in high and low)
    xs=collect(0.5:0.05:1.0);
    ys=collect(0.8:0.02375:0.99);
    lx=length(xs);
    ly=length(ys);
    lx,ly
end

# not needed if run-DIW-1-mod.jl
out2=deserialize("output-DIW-1-mod.jls");

begin
    result_prev=out2[3];
    result = out2[2];
    interp = [[Spline2D(xs, ys, [result[n][i,j][4][4][3][8+s] for i=1:lx, j=1:ly]) for s=1:2] for n=1:2];
    K_prev_spl = [Spline2D(xs, ys, [result[x][i,j][4][3] for i=1:lx, j=1:ly]) for x in AS];
    sol_prev_spl = [[Spline2D(xs, ys, [result[x][i,j][4][4][3][inx] for i=1:lx, j=1:ly]) for inx=1:length(result[x][1,1][4][4][3])] for x in AS];
end;



# convergence of portfolio levels
maximum([maximum([maximum([abs(result[x][i,j][4][4][3][8+s]-result_prev[x][i,j][4][4][3][8+s]) for i=1:lx, j=1:ly]) for s in IS]) for x in AS])
# 8.881784197001252e-15


#convergence of future consumption levels
maximum([maximum([maximum([abs(result[x][i,j][4][4][3][ii]-result_prev[x][i,j][4][4][3][ii]) for i=1:lx, j=1:ly]) for ii=1:8]) for x in AS])
# 2.7755575615628914e-17


#convergence of capital
maximum([maximum([abs(result[x][i,j][4][3]-result_prev[x][i,j][4][3]) for i=1:lx, j=1:ly]) for x in AS])
# 8.881784197001252e-15


#convergence of future distributions
maximum([maximum([maximum(abs.(vcat((result[x][i,j][3]-result_prev[x][i,j][3])...))) for i=1:lx, j=1:ly]) for x in AS])
# 3.3306690738754696e-16


#largest consumption transfer line intercept
maximum([maximum([maximum([result[x][i,j][4][4][3][s] for i=1:lx, j=1:ly]) for s=1:8]) for x in AS])
# 0.0340283881451017

#smallest consumption intercept
maximum([minimum([minimum([result[x][i,j][4][4][3][s] for i=1:lx, j=1:ly]) for s=1:8]) for x in AS])
# -0.051881513696751536

# largest distance between portfolio lines intercepts in high and low state for employed
maximum([abs(result[1][i,j][4][4][3][9]-result[2][i,j][4][4][3][9]) for i=1:lx, j=1:ly])
# 0.17339691920379963

# plot the difference
# corresponds to the left plot in (DIW, Figure 14)
begin
    @gsp xs (100*ys) [(-result[1][i,j][4][4][3][9]+result[2][i,j][4][4][3][9]) for i=1:lx, j=1:ly] "w l t '' lw 2.5 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end

# largest distance between portfolio lines intercepts in high and low state for unemployed
maximum([abs(result[1][i,j][4][4][3][10]-result[2][i,j][4][4][3][10]) for i=1:lx, j=1:ly])
# 0.45637277264062526

#plot the difference
# corresponds to the right plot in (DIW, Figure 14)
begin
    @gsp xs (100*ys) [(-result[1][i,j][4][4][3][10]+result[2][i,j][4][4][3][10]) for i=1:lx, j=1:ly] "w l t '' lw 2.5 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end

## distance from the average to the borrowing limit in state 1
## corresponds to the left plot on (DIW, Figure 9)
begin
    @gsp xs (100*ys) [(xs[i]+result[1][i,j][4][4][3][9]*(1-β)/β) for i=1:lx, j=1:ly] "w l t '' lw 3 lt rgb 'black'";
    @gsp :- xs (100*ys) [(xs[i]*ys[j]+result[1][i,j][4][4][3][10]*(1-β)/β) for i=1:lx, j=1:ly] "w l t '' lw 2 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end


## borrowing bounds in state 2 relative to state 1
# corresponds to the right plot in (DIW, Figure 9)
begin
    @gsp xs (100*ys) [((result[2][i,j][4][4][3][9]*(1-β)/β)-result[1][i,j][4][4][3][9]*(1-β)/β) for i=1:lx, j=1:ly] "w l t '' lw 1.5 lt rgb 'black'";
    @gsp :- xs (100*ys) [((result[2][i,j][4][4][3][10]*(1-β)/β)-result[1][i,j][4][4][3][10]*(1-β)/β) for i=1:lx, j=1:ly] "w l t '' lw 2.5 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end

#plot K as a function of the distribution in aggregate state 1
# corresponds to the left plot in (DIW, Figure 10)
begin
    @gsp xs (100*ys) [result[1][i,j][4][3] for i=1:lx, j=1:ly] "w l t '' lw 2.5 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end

#plot K as a function of the distribution in aggregate state 2 relative to 1
#corresponds to the right plot in (DIW, Figure 10)
begin
    @gsp xs (100*ys) [result[2][i,j][4][3]-result[1][i,j][4][3] for i=1:lx, j=1:ly] "w l t '' lw 2.5 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end


function kernel_check(y::Array{Float64, 1}, xxx::Int64, ccs::Array{Float64,1}, KK::Float64, α::Float64, NN::Array{Float64,1}, XX::Array{Float64,1},  AS::Array{Int64,1}, IS::Array{Int64,1})
    return[(sum([β*(1/(y[4*(s-1)+2*(ξ-1)+σ]/ccs[s]+β*(ρ(xxx, ξ, α, NN, XX)(KK)+1-δ)))*(ρ(xxx, ξ, α, NN, XX)(KK)+1-δ)*atpm[xxx,ξ]*itpm[s,σ,xxx,ξ] for σ in IS, ξ in AS])-1.0) for s in IS ]
end

# largest aberration in the kernel condition at the population averages
[maximum([maximum(abs.(kernel_check(result[xxx][i,j][4][4][3], xxx, [Float64(xs[i]),Float64(xs[i]*ys[j])], result[xxx][i,j][4][3], α, NN, XX,  AS, IS))) for i=1:lx, j=1:ly]) for xxx in AS]
#=julia> 2-element Vector{Float64}:
 5.551115123125783e-16
 1.1102230246251565e-15
=#


# worst violation of the kernel condition at the no-borrwing threshold
[maximum([maximum(abs.(kernel_check(result[xxx][i,j][4][4][3], xxx, [Float64((-result[1][i,j][4][4][3][9]*(1-β)/β)),Float64((-result[1][i,j][4][4][3][10]*(1-β)/β))], result[xxx][i,j][4][3], α, NN, XX,  AS, IS))) for i=1:lx, j=1:ly]) for xxx in AS]
#=julia> 2-element Vector{Float64}:
 0.00040478540743360547
 0.0005353840398942822
=#


# distribution transport from 1 to 1
# corresponds to the left plot on (DIW, Figure 11)
begin
    @gsp xs (100*ys) [result[1][i,j][3][1][1] for i=1:lx, j=1:ly] "w l t '' lw 3.0 lt rgb 'black'";
    @gsp :- xs (100*ys) [result[1][i,j][3][1][2] for i=1:lx, j=1:ly] "w l t '' lw 1.5 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end

# distribution transport from 1 to 2 diff
# corresponds to the right plot on (DIW, Figure 11)
begin
    @gsp xs (100*ys) [result[1][i,j][3][2][1]-result[1][i,j][3][1][1] for i=1:lx, j=1:ly] "w l t '' lw 3.0 lt rgb 'black'";
    @gsp :- xs (100*ys) [result[1][i,j][3][2][2]-result[1][i,j][3][1][2] for i=1:lx, j=1:ly] "w l t '' lw 1.5 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end



# distribution transport from 2 to 1 as a difference to 1-to-1
# corresponds to the left plot on (DIW, Figure 12)
begin
    @gsp xs (100*ys) [result[2][i,j][3][1][1]-result[1][i,j][3][1][1] for i=1:lx, j=1:ly] "w l t '' lw 3.0 lt rgb 'black'";
    @gsp :- xs (100*ys) [result[2][i,j][3][1][2]-result[1][i,j][3][1][2] for i=1:lx, j=1:ly] "w l t '' lw 1.5 lt rgb 'black'";
    @gsp :- "set auto fix";
end




# distribution transport from 2 to 2 diff as a difference to 1-to-2
# corresponds to the right plot on (DIW, Figure 12)
begin
    @gsp xs (100*ys) [result[2][i,j][3][2][1]-result[1][i,j][3][2][1] for i=1:lx, j=1:ly] "w l t '' lw 3.0 lt rgb 'black'";
    @gsp :- xs (100*ys) [result[2][i,j][3][2][2]-result[1][i,j][3][2][2] for i=1:lx, j=1:ly] "w l t '' lw 1.5 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end



## SIMULATION

# define the transport mappings as 2D splines
begin
    dist_spline_1_11=Spline2D(xs, ys, [result[1][i,j][3][1][1] for i=1:lx, j=1:ly]);
    dist_spline_1_12=Spline2D(xs, ys, [result[1][i,j][3][1][2] for i=1:lx, j=1:ly]);
    dist_spline_1_21=Spline2D(xs, ys, [result[1][i,j][3][2][1] for i=1:lx, j=1:ly]);
    dist_spline_1_22=Spline2D(xs, ys, [result[1][i,j][3][2][2] for i=1:lx, j=1:ly]);
    dist_spline_2_11=Spline2D(xs, ys, [result[2][i,j][3][1][1] for i=1:lx, j=1:ly]);
    dist_spline_2_12=Spline2D(xs, ys, [result[2][i,j][3][1][2] for i=1:lx, j=1:ly]);
    dist_spline_2_21=Spline2D(xs, ys, [result[2][i,j][3][2][1] for i=1:lx, j=1:ly]);
    dist_spline_2_22=Spline2D(xs, ys, [result[2][i,j][3][2][2] for i=1:lx, j=1:ly]);
end;

begin
    rng=MersenneTwister(42)
    @time AR0, A0 = gen_fut_dist(1, [0.92,0.89], 1100000, 2, 100, dist_spline_1_11, dist_spline_1_12, dist_spline_1_21, dist_spline_1_22, dist_spline_2_11, dist_spline_2_12, dist_spline_2_21, dist_spline_2_22, atpm);
end;
#1.767172 seconds


#separate high and low states
begin
    A=A0[100001:end];
    AR=AR0[100001:end];
    length(A), length(AR)
    AR1=findall(x->(x==1),AR);
    AR2=findall(x->(x==2),AR);
    length(AR1),length(AR2),length(AR1)+length(AR2)
end

#define the portfolio lines intercepts as 2D splines
begin
    intr_spline_1_1=Spline2D(xs, ys, [result[1][i,j][4][4][3][9:10][1] for i=1:lx, j=1:ly]);
    intr_spline_1_2=Spline2D(xs, ys, [result[1][i,j][4][4][3][9:10][2] for i=1:lx, j=1:ly]);
    intr_spline_2_1=Spline2D(xs, ys, [result[2][i,j][4][4][3][9:10][1] for i=1:lx, j=1:ly]);
    intr_spline_2_2=Spline2D(xs, ys, [result[2][i,j][4][4][3][9:10][2] for i=1:lx, j=1:ly]);
end;


# plot the state of the population by consumption in the last 10,000 periods
# provides the left plot on (DIW, Figure 15)
begin
    @gp [vec[1] for vec in A[end-10000:end]] [vec[2] for vec in A[end-10000:end]] "w p lt -1 pt 5 ps 0.1 t ''"
end

## transform consumption into financial wealth
## VW1 = investment of employed    VW2 = investment of unemployed
begin
    VW1=[(if AR[ii]==1 (intr_spline_1_1(A[ii][1],A[ii][2]/A[ii][1])+(β/(1-β))*A[ii][1])  else (intr_spline_2_1(A[ii][1],A[ii][2]/A[ii][1])+(β/(1-β))*A[ii][1]) end) for ii=(length(A)-10000):length(A)];
    VW2=[(if AR[ii]==1 (intr_spline_1_2(A[ii][1],A[ii][2]/A[ii][1])+(β/(1-β))*A[ii][2])  else (intr_spline_2_2(A[ii][1],A[ii][2]/A[ii][1])+(β/(1-β))*A[ii][2]) end) for ii=(length(A)-10000):length(A)];
end;

# plot the state of the population by investment in the last 10,000 periods
# # provides the right plot on (DIW, Figure 15)
begin
    @gp VW1 VW2 "w p lt rgb 'black' pt 5 ps 0.1 t ''";
end

#portfolio intercept in states 1&2
# not included in DIW
begin
    @gsp xs (100*ys) [result[1][i,j][4][4][3][9] for i=1:lx, j=1:ly] "w l t '' lw 0.75 lt rgb 'black'";
    @gsp :- xs (100*ys) [result[1][i,j][4][4][3][10] for i=1:lx, j=1:ly] "w l t '' lw 0.75 lt rgb 'black'";
    @gsp :- xs (100*ys) [result[2][i,j][4][4][3][9] for i=1:lx, j=1:ly] "w l t '' lw 0.75 lt rgb 'black'";
    @gsp :- xs (100*ys) [result[2][i,j][4][4][3][10] for i=1:lx, j=1:ly] "w l t '' lw 0.75 lt rgb 'black'";
    @gsp :- "set auto fix";
    @gsp :- "set ytics format '%g%%'";
end


## SIMULATION OF CAPITAL

## RECREATING THE RESULTS FROM KRUSELL AND SMITH

## extract the aggregate installed capital from the simulated series

## constrict log-linear regression of installed capital in 1 & 2 agains preceiding capital

begin
    simK=[K_prev_spl[AR[i]](A[i][1],A[i][2]/A[i][1]) for i=1:length(A)];
    log_simK=log.(simK);
end;


begin
    ary1=AR1[2:end];
    ary2=AR2[2:end];
    length(ary1)+length(ary2)
end

@time begin
    data1 = DataFrame(X=log_simK[ary1.-1], Y=log_simK[ary1]);
    data2 = DataFrame(X=log_simK[ary2.-1], Y=log_simK[ary2]);
    ols1 = lm(@formula(Y ~ X), data1);
    ols2 = lm(@formula(Y ~ X), data2);
end;
#0.142809 seconds

## HIGH STATE

coef(ols1)
#=
julia> 2-element Vector{Float64}:
 0.09529355577139297
 0.9374438988630266
=#

stderror(ols1)
#=
julia> 2-element Vector{Float64}:
 1.5223240837069443e-6
 1.0415935050524544e-6
=#

r2(ols1)
#0.9999993843205861

1-r2(ols1)
#6.156794138956201e-7

deviance(ols1)
#0.00023838611399357314

sum((residuals(ols1)).^2) # same as deviance(ols1)
#0.00023838611399357298

sqrt(deviance(ols1)/(length(ary1)-1)) # regression error
#2.1863292642487016e-5


## LOW STATE

coef(ols2)
#=
julia> 2-element Vector{Float64}:
 0.07979444296224014
 0.9418290955002873
=#

stderror(ols2)
#=
julia> 2-element Vector{Float64}:
 1.4236470917501318e-6
 9.896915039403879e-7
=#

r2(ols2)
#0.9999994464737584

1-r2(ols2)
#5.535262416200837e-7

deviance(ols2)
#0.000221678635006257

sum((residuals(ols2)).^2) # same as deviance(ols1)
#0.00022167863500625732

sqrt(deviance(ols2)/(length(ary2)-1)) # regression error
#2.1029066665655758e-5
