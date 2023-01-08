######################################################################################################
#
# Julia code
#
# Illustrates the method described in Sec. 3 of the paper
#       "Another look at the distribution of income and wealth in the Macroeconomy" (DIW) by Andrew Lyasoff
#
# The code provides a verifiable solution to the example from
#                Sec. 18.7 in "Recursive Macroeconomic Theory" (RMT) by Lars Ljungqvist and Thomas Sargent 
#
# This code supplements the paper "Another look at the distribution of income and wealth in the Macroeconomy" (DIW)
#                     by Andrew Lyasoff  (www.andrewlyasoff.tech)
#
# Copyright © 2019-2022 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################

#=
If output-DIW-0.jls is available, there is no need for run-DIW-0.jl 
=#

# this block is superfluous if run-DIW-0.jl
begin
    using Serialization, FileIO
    using LinearAlgebra
    using Interpolations
    using Roots
    using Gnuplot
    using Optim
    include("functions-DIW-0.jl");
    include("ini-setup-RMT-ch18.jl");
    using FastGaussQuadrature
    gqx, gqw = gausslegendre( 100 );
    GQX, GQW = gausslegendre( 10_000 );
end;

# this block is superfluous if run-DIW-0.jl
begin
    saved_vals=deserialize("output-DIW-0.jls");
    SPOT=saved_vals[end]
    clearing=saved_vals[12] 
    g_step=saved_vals[13]
    gsz=saved_vals[14]
    upper_bb=saved_vals[15]
    Fval=saved_vals[16]
    sol=saved_vals[17];
end;

begin
    upper_b=upper_bb;
    cgrid=make_grid(g_step,gsz)
#
    θ=[CubicSplineInterpolation(cgrid,sol[k][:,nos+1],extrapolation_bc = Interpolations.Line()) for k=1:nos];
#
    GF=(0:199)*(upper_b/199)
    F_spline=[CubicSplineInterpolation(GF,Fval[k,:],extrapolation_bc = Interpolations.Line()) for k=1:nos];
    F=[x->(if x<GF[1] 0.0 elseif x>GF[end] 1.0 else F_spline[k](x) end) for k=1:nos];
#
    H=[CubicSplineInterpolation(cgrid,[x+θ[k](x)*SPOT for x in cgrid],extrapolation_bc = Interpolations.Line()) for k=1:nos];
    H_range_bounds=[(cgrid[1]+θ[k](cgrid[1])*SPOT,cgrid[end]+θ[k](cgrid[end])*SPOT) for k=1:nos];
    H_ran_step=[((H_range_bounds[k][2]-H_range_bounds[k][1])/549) for k=1:nos]
    H_range_grids=[H_range_bounds[k][1]:H_ran_step[k]:(H_range_bounds[k][1]+549*H_ran_step[k]) for k=1:nos];
    H_inv_vals=[[find_zero(x->H[k](x)-absc,absc,Order1()) for absc=H_range_grids[k]] for k=1:nos];
    invH=[CubicSplineInterpolation(H_range_grids[k],H_inv_vals[k],extrapolation_bc = Interpolations.Line()) for k=1:nos];
#
    cnext=[CubicSplineInterpolation(cgrid,sol[k0][:,k1],extrapolation_bc = Interpolations.Line()) for k0=1:nos, k1=1:nos];
    next_c_bounds=[(sol[k0][1,k1],sol[k0][end,k1]) for k0=1:nos, k1=1:nos];
    next_c_steps=[(next_c_bounds[k0,k1][2]-next_c_bounds[k0,k1][1])*(1/299) for k0=1:nos, k1=1:nos];
    next_c_range_grid=[next_c_bounds[k0,k1][1]:next_c_steps[k0,k1]:(next_c_bounds[k0,k1][1]+299*next_c_steps[k0,k1]) for k0=1:nos, k1=1:nos];
    inv_c_grid_vals=[[find_zero(x->cnext[k0,k1](x)-absc,absc,Order1()) for absc in next_c_range_grid[k0,k1]] for k0=1:nos, k1=1:nos];
    invcnext=[CubicSplineInterpolation(next_c_range_grid[k0,k1],inv_c_grid_vals[k0,k1],extrapolation_bc = Interpolations.Line()) for k0=1:nos, k1=1:nos];
end;

#market clearing at the last iteration
clearing
#returns -1.7387766322006504e-6

# plot the distribution of all households over the range of consumption
# provides the left plot on Figure 5
begin
    plotgrid=0.0:.001:0.5;
    pval=[F[1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[F[i](x) for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end


# same distribution as a density
# provides the right plot on Figure 5
begin
    plotgrid=cgrid[1]:.001:0.5;
    pval=[Interpolations.gradient(F_spline[1],x)[1] for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=pval=[Interpolations.gradient(F_spline[i],x)[1] for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end

## STATISTICS

begin
    cstep=(cgrid[end]-cgrid[1])*(1/20000);
    gridac=[cgrid[1]+i*cstep for i=0:20000];
    dF=[(F[i].(gridac[2:end])-F[i].(gridac[1:end-1])) for i=1:nos];
    gridacm=(gridac[1:end-1]+gridac[2:end])*(1/2);
    EE=[dot(gridacm,dF[i]) for i=1:nos];
    CM2=[dot((gridacm.-EE[i]).^2,dF[i]) for i=1:nos];
    CM3=[dot((gridacm.-EE[i]).^3,dF[i]) for i=1:nos];
    CM4=[dot((gridacm.-EE[i]).^4,dF[i]) for i=1:nos];
end;

# consumption avergaes across the employment categories
(x->round(x,digits=5)).(EE)
#=
 0.2044
 0.208
 0.21224
 0.21738
 0.2238
 0.23206
 0.24307
=#

# standard deviation
(x->round(x,digits=5)).(sqrt.(CM2))
#=
7-element Vector{Float64}:
 0.0458
 0.04422
 0.04271
 0.04132
 0.04011
 0.03909
 0.03827
=#

#skewness
(x->round(x,digits=5)).([CM3[i]/(CM2[i]^(3/2)) for i=1:nos])
#=
7-element Vector{Float64}:
 0.11558
 0.27624
 0.43391
 0.57253
 0.68702
 0.77908
 0.84976
=#

#kurtosis
(x->round(x,digits=5)).([CM4[i]/(CM2[i]^2) for i=1:nos])
#=
7-element Vector{Float64}:
 3.86671
 3.71633
 3.63201
 3.63853
 3.71073
 3.81284
 3.91739
=#

## DISTRIBUTION OVER ENTERING WEALTH
#
# invert c ↝ c+θᵤ(c)*SPOT-wage*hours[u] so that c ↝ Fᵤ(c) can be expressed θᵤ*SPOT ↝ Fᵤ(c+θᵤ(c)*SPOT-wage*hours[u])
#
# N.B.: 0+θᵤ(0)*SPOT-wage*hours[u] is a theoretical lower bound on entering wealth
# N.B.: the lower bound on entering wealth is θᵤ(0)*ι (ι == face value)
#
begin
    ran_w_bounds=[(cgrid[1]+θ[i](cgrid[1])*SPOT-wage*hours[i],cgrid[end]+θ[i](cgrid[end])*SPOT-wage*hours[i]) for i=1:nos];
    ran_w_steps=[(ran_w_bounds[i][2]-ran_w_bounds[i][1])*(1/2000) for i=1:nos];
    ran_w_range=[(ran_w_bounds[i][1]:ran_w_steps[i]:(ran_w_bounds[i][1]+2000*ran_w_steps[i])) for i=1:nos];
    gcmid=(cgrid[1]+cgrid[end])*(1/2);
    ran_w_vals=[[find_zero(x->x+θ[i](x)*SPOT-wage*hours[i]-absc,gcmid,Order2()) for absc in ran_w_range[i]] for i=1:nos];
    FW_spline=[CubicSplineInterpolation(ran_w_range[i],F[i].(ran_w_vals[i]),extrapolation_bc = Interpolations.Line()) for i=1:nos];
    FW=[x->(if x<ran_w_range[i][1] 0.0 elseif x>ran_w_range[i][end] 1.0 else FW_spline[i](x) end) for i=1:nos];
end;

# actual lower bound on entering wealth
lb_w=minimum([ι*θ[i](0.0) for i=1:nos])
# -1.6885398766990274



# plot the distribution of all households over entering wealth
# not included in [DIW]
begin
    plotgrid=lb_w:.005:5.0;
    pval=[FW[1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 1.5 lt rgb 'black'";
    for i=2:length(FW)
        local pplotgrid, ppval
        pplotgrid=lb_w:.005:5.0;
        ppval=[FW[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 1.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end


# detailed view close to the borrowing limit
# not included in [DIW]
begin
    plotgrid=lb_w:.005:-0.8;
    pval=[FW[1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 1.5 lt rgb 'black'";
    for i=2:length(FW)
        local pplotgrid, ppval
        pplotgrid=lb_w:.005:-0.8;
        ppval=[FW[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 1.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end


# same distribution as above in the form of density
# not included in [DIW]
begin
    plotgrid=(max(ran_w_range[1][1],lb_w)):.005:5.0;
    pval=[Interpolations.gradient(FW_spline[1],x)[1] for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 1.5 lt rgb 'black'";
    for i=2:length(FW)
        local pplotgrid, ppval
        pplotgrid=(max(ran_w_range[i][1],lb_w)):.005:5.0;
        ppval=[Interpolations.gradient(FW_spline[i],x)[1] for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 1.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end



## STATISTICS

begin
    w_min=minimum([max(ran_w_range[i][1],lb_w) for i=1:nos])-0.1
    wstep=(5.0-w_min)*(1/20000);
    gridaw=[w_min+i*wstep for i=0:20000];
    dFW=[(FW[i].(gridaw[2:end])-FW[i].(gridaw[1:end-1])) for i=1:nos];
    gridawm=(gridaw[1:end-1]+gridaw[2:end])*(1/2);
    WEE=[dot(gridawm,dFW[i]) for i=1:nos];
    WCM2=[dot((gridawm.-WEE[i]).^2,dFW[i]) for i=1:nos];
    WCM3=[dot((gridawm.-WEE[i]).^3,dFW[i]) for i=1:nos];
    WCM4=[dot((gridawm.-WEE[i]).^4,dFW[i]) for i=1:nos];
end;

# wealth avergaes across the employment categories
(x->round(x,digits=5)).(WEE) 
#=
7-element Vector{Float64}:
 -0.06015
 -0.04142
 -0.02227
 -0.00194
  0.01976
  0.04307
  0.06916
=#

# standard deviation
(x->round(x,digits=5)).(sqrt.(WCM2))
#=
7-element Vector{Float64}:
 0.97125
 0.9723
 0.97339
 0.97455
 0.97579
 0.97711
 0.97857
=#

#skewness
(x->round(x,digits=5)).([WCM3[i]/(WCM2[i]^(3/2)) for i=1:nos])
#=
7-element Vector{Float64}:
 0.96578
 0.96173
 0.9576
 0.9532
 0.94852
 0.94351
 0.93789
=#

#kurtosis
(x->round(x,digits=5)).([WCM4[i]/(WCM2[i]^2) for i=1:nos])
#=
7-element Vector{Float64}:
 4.03384
 4.02408
 4.01405
 4.00334
 3.99187
 3.97954
 3.96572
=#


## DISTRIBUTION OVER EXITING WEALTH
#
# invert c ↝ θᵤ(c)*SPOT so that c ↝ Fᵤ(c) can be expressed θᵤ*SPOT ↝ Fᵤ(θᵤ*SPOT)
#
begin
    ran_w_ebounds=[(θ[i](cgrid[1])*SPOT,θ[i](cgrid[end])*SPOT) for i=1:nos];
    ran_w_esteps=[(ran_w_ebounds[i][2]-ran_w_ebounds[i][1])*(1/2000) for i=1:nos];
    ran_w_erange=[(ran_w_ebounds[i][1]:ran_w_esteps[i]:(ran_w_ebounds[i][1]+2000*ran_w_esteps[i])) for i=1:nos];
    gcmid=(cgrid[1]+cgrid[end])*(1/2);
    ran_w_evals=[[find_zero(x->θ[i](x)*SPOT-absc,gcmid,Order2()) for absc in ran_w_erange[i]] for i=1:nos];
    FWX_spline=[CubicSplineInterpolation(ran_w_erange[i],F[i].(ran_w_evals[i]),extrapolation_bc = Interpolations.Line()) for i=1:nos];
    FWX=[x->(if x<ran_w_erange[i][1] 0.0 elseif x>ran_w_erange[i][end] 1.0 else FWX_spline[i](x) end) for i=1:nos];
end;

lb_xw=minimum([ran_w_erange[i][1] for i=1:nos])
# -1.6274869069693472

# plot the distribution of all households over the range of exiting wealth
# not included in [DIW]
begin
    plotgrid=lb_xw:.005:5.0;
    pval=[FWX[1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 1.5 lt rgb 'black'";
    for i=2:length(FWX)
        local pplotgrid,ppval
        pplotgrid=lb_xw:.005:5.0;
        ppval=[FWX[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 1.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end

# detailed view of the above near the borrowing limit
# mot included in [DIW]
begin
    plotgrid=(ran_w_erange[1][1]):.001:-0.8;
    pval=[FWX[1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 1.5 lt rgb 'black'";
    for i=2:length(FWX)
        local pplotgrid,ppval
        pplotgrid=(ran_w_erange[i][1]):.001:-0.8;
        ppval=[FWX[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 1.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end


# detailed view of exiting AND entering wealth near the borrowing limit
# provides the left plot on Figure 7
begin
    plotgrid=(ran_w_erange[1][1]):.005:-0.6;
    pval=[FW[1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(FW)
        local pplotgrid, ppval
        pplotgrid=(ran_w_erange[i][1]):.005:-0.6;
        ppval=[FW[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    plotgrid=(ran_w_erange[1][1]):.01:-0.6;
    pval=[FWX[1](x) for x in plotgrid];
    @gp :- plotgrid pval "w p lt rgb 'black' pt 7 ps 1.25 t ''";
    for i=2:length(FWX)
        local pplotgrid,ppval
        pplotgrid=(ran_w_erange[i][1]):.01:-0.6;
        ppval=[FWX[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w p lt rgb 'black' pt 7 ps 1.25 t ''";
    end
    @gp :- "set auto fix";
end



# distribution over entering AND exiting wealth 
# Provides the left plot on Figure 4.
begin
    plotgrid=ran_w_erange[1][1]:.05:3.5;
    pval=[FWX[1](x) for x in plotgrid];
    @gp plotgrid pval "w p lt rgb 'black' pt 7 ps 1.25 t ''";#"w l t '' lw 1.5 lt rgb 'black'";
    for i=2:length(FWX)
        local pplotgrid, ppval
        pplotgrid=ran_w_erange[i][1]:.05:3.5;
        ppval=[FWX[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w p lt rgb 'black' pt 7 ps 1.25 t ''";#"w l t '' lw 1.5 lt rgb 'black'";
    end
    for i=1:length(FW)
        local pplotgrid,ppval
        pplotgrid=ran_w_erange[i][1]:.004:3.5;
        ppval=[FW[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 1.85 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end



# distribution over exiting wealth as a density
# not included in [DIW]
begin
    plotgrid=(max(ran_w_erange[1][4],lb_w)):.005:5.0;
    pval=[Interpolations.gradient(FWX_spline[1],x)[1] for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 1.5 lt rgb 'black'";
    for i=2:length(FWX)
        local pplotgrid, ppval
        pplotgrid=(max(ran_w_erange[i][4],lb_w)):.005:5.0;
        ppval=[Interpolations.gradient(FWX_spline[i],x)[1] for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 1.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end



# distribution over entering AND exiting wealth as a density
# provides the right plot on Figure 7
begin
    plotgrid=(max(ran_w_erange[1][3],lb_w)):.005:5.0;
    pval=[Interpolations.gradient(FW_spline[1],x)[1] for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(FW)
        local pplotgrid, ppval
        pplotgrid=(max(ran_w_erange[i][3],lb_w)):.005:5.0;
        ppval=[Interpolations.gradient(FW_spline[i],x)[1] for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    plotgrid=(max(ran_w_erange[1][3],lb_w)):.025:5.0;
    pval=[Interpolations.gradient(FWX_spline[1],x)[1] for x in plotgrid];
    @gp :- plotgrid pval "w p lt rgb 'black' pt 7 ps 1.0 t ''";
    for i=2:length(FWX)
        local pplotgrid, ppval
        pplotgrid=(max(ran_w_erange[i][3],lb_w)):.025:5.0;
        ppval=[Interpolations.gradient(FWX_spline[i],x)[1] for x in pplotgrid];
        @gp :- pplotgrid ppval "w p lt rgb 'black' pt 7 ps 1.0 t ''";
    end
    @gp :- "set auto fix";
end



## STATISTICS

begin
    w_min=minimum([max(ran_w_erange[i][1],lb_w) for i=1:nos])-0.1
    wstep=(5.0-w_min)*(1/20000);
    gridaw=[w_min+i*wstep for i=0:20000];
    dFWX=[(FWX[i].(gridaw[2:end])-FWX[i].(gridaw[1:end-1])) for i=1:nos];
    gridawm=(gridaw[1:end-1]+gridaw[2:end])*(1/2);
    WEE=[dot(gridawm,dFWX[i]) for i=1:nos];
    WCM2=[dot((gridawm.-WEE[i]).^2,dFWX[i]) for i=1:nos];
    WCM3=[dot((gridawm.-WEE[i]).^3,dFWX[i]) for i=1:nos];
    WCM4=[dot((gridawm.-WEE[i]).^4,dFWX[i]) for i=1:nos];
end;

# wealth avergaes across the employment categories
(x->round(x,digits=5)).(WEE) 
#=
7-element Vector{Float64}:
 -0.20375
 -0.15901
 -0.09993
 -0.01888
  0.09461
  0.25604
  0.48922
=#

# standard deviation
(x->round(x,digits=5)).(sqrt.(WCM2))
#=
7-element Vector{Float64}:
 0.92818
 0.93036
 0.93253
 0.93465
 0.93657
 0.93805
 0.93861
=#

#skewness
(x->round(x,digits=5)).([WCM3[i]/(WCM2[i]^(3/2)) for i=1:nos])
#=
7-element Vector{Float64}:
 1.00653
 0.99688
 0.98671
 0.97588
 0.96403
 0.95015
 0.9316
=#

#kurtosis
(x->round(x,digits=5)).([WCM4[i]/(WCM2[i]^2) for i=1:nos])
#=
7-element Vector{Float64}:
 4.14594
 4.12455
 4.10058
 4.07227
 4.0367
 3.98815
 3.91411
=#





#################################################


# investment as a function of consumption in all 7 states
# not included in [DIW]
begin
    plotgrid=cgrid[1]:.001:cgrid[end];
    pval=[θ[1](x)*SPOT for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[θ[i](x)*SPOT for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
     end
    @gp :- "set auto fix";
end


# investment as a function of consumption in all 7 states
# provides the left plot on Figure 6
begin
    plotgrid=cgrid[1]:.001:0.4;
    pval=[θ[1](x)*SPOT for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 1.25 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[θ[i](x)*SPOT for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 1.25 lt rgb 'black'";
     end
    @gp :- "set auto fix";
end


# investment as a function of consumption in all 7 states near 0
# not included in [DIW]
begin
    #plotgrid=cgrid[1]:.0005:0.2;
    plotgrid=0.0:.0005:0.2;
    pval=[θ[1](x)*SPOT for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[θ[i](x)*SPOT for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
     end
    @gp :- "set auto fix";
end


# investment as a function of consumption in all 7 states extremely close to 0
# provides the right plot on Figure 6
begin
    #plotgrid=cgrid[1]:.0005:0.2;
    plotgrid=0.0:.0001:0.06;
    pval=[θ[1](x)*SPOT for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[θ[i](x)*SPOT for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
     end
    @gp :- "set auto fix";
end

# the slope of the investment near $0$ (cgrid[1]==0.006104783782620999)
(x->round(x,digits=5)).([Interpolations.gradient(θ[i],cgrid[1])[1] for i=1:nos]*SPOT)
#=
7-element Vector{Float64}:
 0.40615
 0.34429
 0.28827
 0.23839
 0.19469
 0.15701
 0.12502
=#

# marginal propenisty to consume in terms of exiting wealth
(x->round(x,digits=5)).([(1/(SPOT*Interpolations.gradient(θ[i],cgrid[1])[1])) for i=1:nos])
#=
7-element Vector{Float64}:
 2.46213
 2.90454
 3.46896
 4.19479
 5.13634
 6.36905
 7.99873
=#

# values at 0 (via extrapolation)
[θ[i](0.0) for i=1:7]
#=
7-element Vector{Float64}:
 -7.742462265862489
 -7.742440870497814
 -7.74242451293142
 -7.742412748758499
 -7.742404764730136
 -7.742399614117026
 -7.742396438820142
=#

# values at the first point on the grid
[θ[i](cgrid[1]) for i=1:7]
#=
7-element Vector{Float64}:
 -7.730672283746217
 -7.732446700451867
 -7.734056434204218
 -7.735492624443796
 -7.736753172403883
 -7.737841867581036
 -7.738767301581213
=#


###  future consumption as a function of entering wealth

begin
    ran_w_bounds_next=[(θ[i](cgrid[1])*ι,θ[i](cgrid[end])*ι) for i=1:nos,j=1:nos];
    ran_w_steps_next=[(ran_w_bounds_next[i,j][2]-ran_w_bounds_next[i,j][1])*(1/2000) for i=1:nos,j=1:nos];
    ran_w_range_next=[(ran_w_bounds_next[i,j][1]:ran_w_steps_next[i,j]:(ran_w_bounds_next[i,j][1]+2000*ran_w_steps_next[i,j])) for i=1:nos,j=1:nos];
    gcmid=(cgrid[1]+cgrid[end])*(1/2);
    ran_w_vals_next=[[find_zero(x->θ[i](x)*ι-absc,gcmid,Order2()) for absc in ran_w_range_next[i,j]] for i=1:nos,j=1:nos];
    cnw_spline=[CubicSplineInterpolation(ran_w_range_next[i,j],cnext[i,j].(ran_w_vals_next[i,j]),extrapolation_bc = Interpolations.Line()) for i=1:nos,j=1:nos];
end;


# as a function of entering wealth from state 2
# not included in [DIW]
begin
    plotgrid=ran_w_range_next[2,1][1]:.005:3.0;
    pval=[cnw_spline[2,1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 1.5 lt rgb 'black'";
    for j=2:nos
        local pplotgrid,ppval
        pplotgrid=ran_w_range_next[2,j][1]:.005:3.0;
        ppval=[cnw_spline[2,j](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 1.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end


# as a function of enetering wealth from state 1
# provides the first part of the right plot in Figure 8
begin
    plotgrid=ran_w_range_next[1,1][1]:.001:3;
    pval=[cnw_spline[1,1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for j=2:nos
        local pplotgrid,ppval
        pplotgrid=ran_w_range_next[1,j][1]:.001:3;
        ppval=[cnw_spline[1,j](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end

# in addition, as a function of entering wealth from state 7
# provides the second part of the right plot in Figure 8
begin
    plotgrid=ran_w_range_next[7,1][1]:.05:3;
    pval=[cnw_spline[7,1](x) for x in plotgrid];
    @gp :- plotgrid pval "w p lt rgb 'black' pt 7 ps 1.0 t ''";
    for j=2:nos
        local pplotgrid,ppval
        pplotgrid=ran_w_range_next[7,j][1]:.05:3;
        ppval=[cnw_spline[7,j](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w p lt rgb 'black' pt 7 ps 1.0 t ''";
    end
end


# in addition as a function of entering wealth in state 2
# not included in [DIW]
begin
    plotgrid=ran_w_range_next[2,1][1]:.05:3;
    pval=[cnw_spline[2,1](x) for x in plotgrid];
    @gp :- plotgrid pval "w p lt rgb 'black' pt 6 ps 1.0 t ''";#"w l t '' lw 2.5 lt rgb 'red'";
    for j=2:nos
        local pplotgrid,ppval
        pplotgrid=ran_w_range_next[2,j][1]:.05:3;
        ppval=[cnw_spline[2,j](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w p lt rgb 'black' pt 6 ps 1.0 t ''";#"w l t '' lw 2.5 lt rgb 'red'";
    end
    @gp :- "set auto fix";
end


# in addition as a function of entering wealth in state 4
# not included in [DIW]
begin
    plotgrid=ran_w_range_next[4,1][1]:.05:3.0;
    pval=[cnw_spline[4,1](x) for x in plotgrid];
    @gp :- plotgrid pval "w p lt rgb 'black' pt 5 ps 1.0 t ''";#"w l t '' lw 2.5 lt rgb 'red'";
    for j=2:nos
        local pplotgrid,ppval
        pplotgrid=ran_w_range_next[4,j][1]:.05:3.0;
        ppval=[cnw_spline[2,j](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w p lt rgb 'black' pt 5 ps 1.0 t ''";#"w l t '' lw 2.5 lt rgb 'red'";
    end
    @gp :- "set auto fix";
end


#
# future consumption in all 7 states as a function of present consumption in (present) state 1
# provides the first part of the left plot in Figure 8
begin
    plotgrid=cgrid[1]:.01:cgrid[end];
    pval=[cnext[1,1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[cnext[1,i](x) for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end



# in addition to the above
# future consumption in all 7 states as a function of present consumption in (present) state 7
# provides the second part fo the left plot in Figure 8
begin
    plotgrid=cgrid[1]:.01:cgrid[end];
    pval=[cnext[7,1](x) for x in plotgrid];
    @gp :- plotgrid pval "w p lt rgb 'black' pt 7 ps 1.0 t ''";
    for i=2:length(F)
        local ppval
        ppval=[cnext[7,i](x) for x in plotgrid];
        @gp :- plotgrid ppval "w p lt rgb 'black' pt 7 ps 1.0 t ''";
    end
end


###############################################################################



## CONSUMPTION AS A FUNCTION OF TOTAL WEALTH
#
# invert c ↝ c + θᵤ(c)*SPOT]
#
begin
    ran_w_cbounds=[(cgrid[1]+θ[i](cgrid[1])*SPOT,cgrid[end]+θ[i](cgrid[end])*SPOT) for i=1:nos];
    ran_w_csteps=[(ran_w_cbounds[i][2]-ran_w_cbounds[i][1])*(1/2000) for i=1:nos];
    ran_w_crange=[(ran_w_cbounds[i][1]:ran_w_csteps[i]:(ran_w_cbounds[i][1]+2000*ran_w_csteps[i])) for i=1:nos];
    gcmid=(cgrid[1]+cgrid[end])*(1/2);
    ran_w_cvals=[[find_zero(x->x+θ[i](x)*SPOT-absc,gcmid,Order2()) for absc in ran_w_crange[i]] for i=1:nos];
    CW_spline=[CubicSplineInterpolation(ran_w_crange[i],ran_w_cvals[i],extrapolation_bc = Interpolations.Line()) for i=1:nos];
    #CW=[x->(if x<ran_w_crange[i][1] 0.0 elseif x>ran_w_crange[i][end] 1.0 else CW_spline[i](x) end) for i=1:nos];
end;

# consumption as a function of total wealth
# generates the left plot in Figure 9
begin
    plotgrid=(ran_w_crange[1][1]):.005:0.5;
    pval=[CW_spline[1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(FW)
        local pplotgrid, ppval
        pplotgrid=(ran_w_crange[i][1]):.005:0.5;
        ppval=[CW_spline[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end


# marginal propensity to consume in terms of total wealth (consumption plus investment)
# generates the right plot in Figure 9
begin
    plotgrid=(ran_w_crange[1][1]):.005:0.5;
    pval=[Interpolations.gradient(CW_spline[1],x)[1] for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(FW)
        local pplotgrid, ppval
        pplotgrid=(ran_w_crange[i][1]):.005:0.5;
        ppval=[Interpolations.gradient(CW_spline[i],x)[1] for x in pplotgrid];
        @gp :- pplotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end




###############################################################################


## compute the endogenous bounds for the wealth at grid
(minimum([θ[k](cgrid[1]) for k=1:nos])*SPOT,maximum([θ[k](cgrid[end]) for k=1:nos])*SPOT)
# returns: (-1.6274869069693472, 17.937506971712203)

# endogenous borrowing limit for shares held at grid
minimum([θ[k](cgrid[1]) for k=1:nos])
# returns: -7.738767301581213

## check the distribution support at grid
(maximum([F[k](cgrid[1]) for k=1:nos]) ,minimum([F[k](cgrid[end]) for k=1:nos]))
# returns: (5.0147289107981385e-6, 0.9999999999999976)

## choosing a narrower interval where most of the distribution is amassed
(maximum([F[k](0.0065) for k=1:nos]),minimum([F[k](0.5) for k=1:nos]))
# returns: (5.990611174450551e-6, 0.9999961067541274)

## compute the corresponding endogenous bounds on the wealth
w_bounds=(minimum([θ[k](0.0065) for k=1:nos])*SPOT,maximum([θ[k](0.5) for k=1:nos])*SPOT)
# returns: (-1.6274374967927319, 7.182671906288613)

## the spread of the wealth 
w_bounds[end]-w_bounds[1]
#returns: 8.810109403081345

## the equilibrium interest rate is
true_r=ι/SPOT-1.0
#returns: 0.037018510722689246

################################################
#
# The rest of the code is illustrates the compatibility with the classical method:
# if the latter is initiated with the true r AND the true distribution  
#         it produces essentially the same result
#
###############################################

#the "ad hoc" borrowing limit from RMT
-min(3.0, wage*hours[1]/true_r)
# returns: -1.6272627182032762

# endogenous borrwoing limit for each employment class from DIW
([θ[k](cgrid[1]) for k=1:nos])*SPOT
#=
7-element Vector{Float64}:
 -1.6257844994637685
 -1.626157664832801
 -1.6264961968629765
 -1.6267982321509935
 -1.627063329319664
 -1.627292285311936
 -1.6274869069693472
=#

# endogenous borrwoing limit for each employment class from DIW
# uses extrapolation
([θ[k](0.0) for k=1:nos])*SPOT
#=
7-element Vector{Float64}:
 -1.6282639694852685
 -1.628259469973328
 -1.6282560299260829
 -1.6282535558838156
 -1.6282518768177423
 -1.6282507936277364
 -1.6282501258528954
=#

minimum(([θ[k](0.0) for k=1:nos])*SPOT)
# returns: -1.6282639694852685

# "ad hoc" from RMT minus truly endogenous from DIW
-min(3.0, wage*hours[1]/true_r)-minimum(([θ[k](0.0) for k=1:nos])*SPOT)
# returns: 0.0010012512819923547

# generate discrete initial distribution from the the true EXITING distribution 
begin
    grdcc=cgrid[1]:0.00025:(cgrid[1]+1999*0.00025);
    grdccx=cgrid[1]:0.00025:(cgrid[1]+2000*0.00025);
    grdaax=[θ[k](x)*SPOT for k=1:7, x in grdccx];
    grdaat=[θ[k](x)*SPOT for k=1:7, x in (grdcc+grdccx[2:end])*(1/2)];
    true_dist=[(FWX[k](grdaax[k,i])-FWX[k](grdaax[k,i-1]))*BSP[k] for k=1:7, i=2:2001];
end;

include("functions-RMT-ch18.jl");

# run the procedure described in RMT (may take some time -- runs on a single CPU)
@time Λ1t,tbl2t=LS_method(true_r,5,5,hours,wage,β,grdaat,true_dist,CPROB,R);
# 150.117135 seconds (23.68 M allocations: 294.403 GiB, 2.95% gc time);

begin
    aassetst=[grdaat[k,Int64(tbl2t[k,i])] for k=1:nos, i=1:2000];
    SLt=[sum(Λ1t[state,:]) for state=1:nos];
    LFt=[[Λ1t[state,1]/SLt[state]] for state=1:nos];
    for i=2:2000
        for state=1:nos
            LFt[state]=vcat(LFt[state],LFt[state][end]+Λ1t[state,i]/SLt[state]);
        end;
    end;
end;

# provides the first part of the right plot in Figure 4
begin
    @gp aassetst[1,1:10:1300] LFt[1][1:10:1300] "w l t '' lw 2.5 lt rgb 'black'";#"w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- "set auto fix";
    @gp :- aassetst[2,1:10:1300] LFt[2][1:10:1300] "w l t '' lw 2.5 lt rgb 'black'";#"w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[3,1:10:1300] LFt[3][1:10:1300] "w l t '' lw 2.5 lt rgb 'black'";#"w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[4,1:10:1300] LFt[4][1:10:1300] "w l t '' lw 2.5 lt rgb 'black'";#"w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[5,1:10:1300] LFt[5][1:10:1300] "w l t '' lw 2.5 lt rgb 'black'";#"w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[6,1:10:1300] LFt[6][1:10:1300] "w l t '' lw 2.5 lt rgb 'black'";#"w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[7,1:10:1300] LFt[7][1:10:1300] "w l t '' lw 2.5 lt rgb 'black'";#"w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- "set auto fix";
end

# provides the second part of the right plot in Figure 4
begin
    plotgrid=aassetst[1,1]:.05:aassetst[1,1300];
    pval=[FWX[1](x) for x in plotgrid];
    @gp :- plotgrid pval "w p lt rgb 'black' pt 7 ps 1.25 t ''";#"w l t '' lw 2.5 lt rgb 'black'";
     for i=2:length(FW)
        local pplotgrid,ppval
         pplotgrid=aassetst[i,1]:.05:aassetst[i,1300];
         ppval=[FWX[i](x) for x in pplotgrid];
        @gp :- pplotgrid ppval "w p lt rgb 'black' pt 7 ps 1.25 t ''";#"w l t '' lw 2.5 lt rgb 'black'";
    end
end

# market clearing produced with the classical method
begin
    demandt=[grdaat[kk,Int64(tbl2t[kk,ii])] for kk=1:nos, ii=1:2000];
    sum(Λ1t.*demandt)
end
#returns: -0.006734619163051014
