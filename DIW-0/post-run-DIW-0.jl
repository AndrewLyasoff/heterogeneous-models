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
#                                        by Andrew Lyasoff
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
end;

# this block is superfluous if run-DIW-0.jl
begin
    saved_vals=deserialize("output-DIW-0.jls");
    SPOT=saved_vals[end]
    clearing=saved_vals[12] 
    g_step=saved_vals[13]
    gsz=saved_vals[14]
    upper_b=saved_vals[15]
    Fval=saved_vals[16]
    sol=saved_vals[17];
#
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

#restore the grid
gridc=make_grid(g_step,gsz);

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
    plotgrid=gridc[1]:.001:0.5;
    pval=[Interpolations.gradient(F_spline[1],x)[1] for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=pval=[Interpolations.gradient(F_spline[i],x)[1] for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end


# plot the distribution of all households over the range of asset holdings (wealth)
# provides the left plot on Figure 6
begin
    plotgrid=gridc[1]:.001:0.43;
    pval=[F[1](x) for x in plotgrid];
    xval=[θ[1](x)*SPOT for x in plotgrid];
    @gp xval pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval, xxval
        ppval=[F[i](x) for x in plotgrid];
        xxval=[θ[i](x)*SPOT for x in plotgrid];
        @gp :- xxval ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end


# same distribution as above in the form of density
# provides the right plot on Figure 6
begin
    plotgrid=gridc[1]:.001:0.43;
    pval=[Interpolations.gradient(F_spline[1],x)[1] for x in plotgrid];
    xval=[θ[1](x)*SPOT for x in plotgrid];
    @gp xval pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval, xxval
        ppval=[Interpolations.gradient(F_spline[i],x)[1] for x in plotgrid];
        xxval=[θ[i](x)*SPOT for x in plotgrid];
        @gp :- xxval ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end


# investment as a function of consumption in all 7 states
# provides the left plot on Figure 7
begin
    plotgrid=gridc[1]:.001:gridc[end];
    pval=[θ[1](x)*SPOT for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[θ[i](x)*SPOT for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
     end
    @gp :- "set auto fix";
end


# investment as a function of consumption in all 7 states near 0
# provides the right plot on Figure 7
begin
    plotgrid=gridc[1]:.0005:0.2;
    pval=[θ[1](x)*SPOT for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[θ[i](x)*SPOT for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
     end
    @gp :- "set auto fix";
end

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

#
# future consumption in all 7 states as a function of present consumption in (present) state 2
# provides the left plot in Figure 8
begin
    plotgrid=gridc[1]:.001:gridc[end];
    pval=[cnext[2,1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[cnext[2,i](x) for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end



#
# future consumption in all 7 states as a function of present consumption in (present) state 7
# provides the right plot in Figure 8
begin
    plotgrid=gridc[1]:.001:gridc[end];
    pval=[cnext[7,1](x) for x in plotgrid];
    @gp plotgrid pval "w l t '' lw 2.5 lt rgb 'black'";
    for i=2:length(F)
        local ppval
        ppval=[cnext[7,i](x) for x in plotgrid];
        @gp :- plotgrid ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
    @gp :- "set auto fix";
end




## compute the endogenous bounds for the wealth at grid
(minimum([θ[k](gridc[1]) for k=1:nos])*SPOT,maximum([θ[k](gridc[end]) for k=1:nos])*SPOT)
# returns: (-1.6274869069693472, 17.937506971712203)

# endogenous borrowing limit for shares held at grid
minimum([θ[k](gridc[1]) for k=1:nos])
# returns: -7.738767301581213

## check the distribution support at grid
(maximum([F[k](gridc[1]) for k=1:nos]) ,minimum([F[k](gridc[end]) for k=1:nos]))
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
([θ[k](gridc[1]) for k=1:nos])*SPOT
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

# generate discrete initial distribution from the the true one
begin
    grdcc=cgrid[1]:0.00025:(cgrid[1]+1999*0.00025);
    grdccx=cgrid[1]:0.00025:(cgrid[1]+2000*0.00025);
    grdaat=[θ[k](x)*SPOT for k=1:7, x in grdcc];
    true_dist=[(F[k](grdccx[i])-F[k](grdccx[i-1]))*BSP[k] for k=1:7, i=2:2001];
end;

include("functions-RMT-ch18.jl");

# run the procedure described in RMT (may take some time -- runs on a single CPU)
@time Λ1t,tbl2t=LS_method(true_r,5,5,nos,hours,wage,β,grdaat,true_dist,2000,CPROB,R);
# 155.930313 seconds (23.68 M allocations: 294.404 GiB, 2.80% gc time)

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

# provides the first part of Figure 4
begin
    @gp aassetst[1,2:10:1300] LFt[1][2:10:1300] "w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- "set auto fix";
    @gp :- aassetst[2,2:10:1300] LFt[2][2:10:1300] "w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[3,2:10:1300] LFt[3][2:10:1300] "w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[4,2:10:1300] LFt[4][2:10:1300] "w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[5,2:10:1300] LFt[5][2:10:1300] "w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[6,2:10:1300] LFt[6][2:10:1300] "w p lt rgb 'black' pt 7 ps 0.75 t ''";
    @gp :- aassetst[7,2:10:1300] LFt[7][2:10:1300] "w p lt rgb 'black' pt 7 ps 0.75 t ''";
end

# provides the second part of Figure 4
begin
    plotgrid=cgrid[1]:.001:0.39;
    pval=[F[1](x) for x in plotgrid];
    xval=[θ[1](x)*SPOT for x in plotgrid];
    for i=2:length(F)
        local ppval, xxval
        ppval=[F[i](x) for x in plotgrid];
        xxval=[θ[i](x)*SPOT for x in plotgrid];
        @gp :- xxval ppval "w l t '' lw 2.5 lt rgb 'black'";
    end
end

# market clearing produced with the classical method
begin
    demandt=[grdaat[kk,Int64(tbl2t[kk,ii])] for kk=1:nos, ii=1:2000];
    sum(Λ1t.*demandt)
end
#returns: -0.009249065617648822
