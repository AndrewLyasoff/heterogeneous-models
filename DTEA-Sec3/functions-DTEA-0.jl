######################################################################################################
#
# Julia code
#
# Functions to implement the method described in Sec. 3 of the paper
#       "Dynamic transportation of economic agents" [DTEA] by Andrew Lyasoff
#
# The code provides a verifiable solution to the example from
#                Sec. 18.2 & 18.7 in "Recursive Macroeconomic Theory" (RMT) by Lars Ljungqvist and Thomas Sargent 
#
# This code supplements the paper "Dynamic transportation of economic agents" [DTEA]
#                                        by Andrew Lyasoff (www.andrewlyasoff.tech)
#
# Copyright © 2019-2023 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################

function FGQT(f,a,b)
    return (arg->(1/2)*(b-a)*f( (1/2)*(b-a)*arg + (1/2)*(a+b) ))
end


function dU(y::Float64,R::Float64)
    return 1/(y)^(R)
end;


function make_grid(step::Float64,size::Int64)
    return step:step:(step+(size-1)*step)
end

function sumU(nos::Int64,yyy::Float64,sno::Int64,ccc::Float64,spot::Float64,invh,CPROB::Array{Float64,2},inc::Array{Float64,1},ι::Float64,R::Float64)
    return sum([exp(-ρ)*dU(invh[k](inc[k]+yyy*ι)/ccc,R)*(ι)*CPROB[sno,k] for k=1:nos])-spot
end;

## for c ∈ grid, returns the solution *θᵤ(c) to equation (3.2); the argument invh is the inverse mapping in the denominator of (3.2)
function zroU(nos::Int64,sno::Int64,ccc::Float64,spot::Float64,invh,CPROB::Array{Float64,2},inc::Array{Float64,1},ι::Float64,R::Float64)
    local g, zro, zroPrev, ZRO
    g=(yyy->sum([0.96*dU(invh[k](inc[k]+yyy*ι)/ccc,R)*(ι)*CPROB[sno,k] for k=1:nos])-spot)
    if abs(g(0))<1e-14
        ZRO=0.0
    elseif g(0)>=1e-14
        zro=0.0;
        iii=0;
        while g(zro)>=0&&iii<100000
            iii+=1
            zroPrev=zro
            zro-=0.01/(g(zro+0.0001)-g(zro-0.0001))
        end
        if iii>=100000
            ZRO=NaN
        else
            ZRO=find_zero(g,(zroPrev,zro))
        end
    else
        zro=0.0;
        iii=0
        while g(zro)<=0&&iii<100000
            iii+=1
            zroPrev=zro
            zro+=max(min((1/10^6)*0.0002/(g(zro+0.0001)-g(zro-0.0001)),-0.001),-0.05)
        end
        if iii>=100000
            ZRO=NaN
        else
            ZRO=find_zero(g,(zro,zroPrev))
        end
    end
    return ZRO
end;

## implements Step 1 in 3.2; returns the upper bound on the effective grid over consumption
function FindUpperC(nos::Int64,spot::Float64,invh,CPROB::Array{Float64,2},inc::Array{Float64,1},ι::Float64,R::Float64)
    local tmp_bound, loc_inx, loc_iinx, loc_iiinx
    loc_inx=0;
    tmp_bound=loc_inx+0.0001;
    while tmp_bound>loc_inx && loc_inx<100100
        loc_inx+=10
        tmp_lst=[zroU(nos,k,Float64(loc_inx),spot,invh,CPROB,inc,ι,R) for k=1:nos];
        tmp_bound=maximum([invh[k](inc[k]+tmp_lst[k]*ι) for k=1:nos]);
    end
    if loc_inx<=10000
        loc_iinx=loc_inx-10;
        tmp_bound=loc_iinx+0.0001;
        while tmp_bound>Float64(loc_iinx) && loc_iinx<loc_inx
            loc_iinx+=1
            tmp_lst=[zroU(nos,k,Float64(loc_iinx),spot,invh,CPROB,inc,ι,R) for k=1:nos];
            tmp_bound=maximum([invh[k](inc[k]+tmp_lst[k]*ι) for k=1:nos]);
        end
        loc_iiinx=loc_iinx-1;
        tmp_bound=loc_iiinx+0.0001;
        while tmp_bound>Float64(loc_iiinx) && loc_iiinx<loc_iinx
            loc_iiinx+=1/10
            tmp_lst=[zroU(nos,k,Float64(loc_iiinx),spot,invh,CPROB,inc,ι,R) for k=1:nos];
            tmp_bound=maximum([invh[k](inc[k]+tmp_lst[k]*ι) for k=1:nos]);
        end
        return loc_iiinx
    else
        return NaN
    end
end;

#=
returns the maximal error in (3.1) & (3.2) for every c ∈ grid and for every employment state,
i.e., the maximal error in the supplied solution x, which is the list of future consumptions
(total of nos=7) followed by the present portfolio (a scalar)
=#
function gtest(nos::Int64,state::Int64,pflo,spot::Float64,ccc::Float64,x::Array{Float64,1},CPROB::Array{Float64,2},inc::Array{Float64,1},ι::Float64)
    return max(maximum([abs(x[i]+pflo[i](x[i])*spot-inc[i]-x[nos+1]*ι) for i=1:nos]), abs(sum([exp(-ρ)*dU(x[l]/ccc,R)*(ι)*CPROB[state,l] for l=1:nos])-spot))
end;

#=
solves equation (3.2) for *θᵤ(c);
returnes the denominators in (3.2) = future consumptions computed at c ∈ grid, followed by *θᵤ(c) itself 
sol = list of future consumptions Τᵤᵛ(c) followed by present portfolio *θᵤ(c) for c ∈ grid
(one list for every state u)
RETURNS: 'sol',  the largest error with which 'sol' satisfies the kernel condition (3.2) for all employment states and abscissas on the consumption grid; the largest error with which 'sol' satisfies all first order
conditions (3.1) and (3.2) for all states and abscissas
=#
function solve_sys(nos::Int64,g_step::Float64,gsz::Int64,spot::Float64,theta,invh,CPROB::Array{Float64,2},inc::Array{Float64,1},ι::Float64,R::Float64)
    local gg_step, cgrid, sol, accu
    gg_step=g_step
    cgrid=make_grid(gg_step,gsz)
    sol=[zeros(gsz,nos+1) for k=1:nos]
    accu=[zeros(gsz) for k=1:nos]
        for k=1:nos
            for i=1:gsz
                yloc=zroU(nos,k,Float64(cgrid[i]),spot,invh,CPROB,inc,ι,R)
                accu[k][i]=sumU(nos,yloc,k,Float64(cgrid[i]),spot,invh,CPROB,inc,ι,R)
                sol[k][i,:]=vcat([invh[kk](inc[kk]+yloc*ι) for kk=1:nos],yloc)   #* k ∼ s , kk ∼ σ
            end
        end
    return sol, maximum(maximum([abs.(accu[k]) for k=1:nos])),maximum([maximum([gtest(nos,k0,theta,spot,cgrid[i],sol[k0][i,:],CPROB,inc,ι) for i=1:gsz]) for k0=1:nos])
end;

#=
RETURNS: *θᵤ(·) as a spline over c ∈ grid for every s;
the future consumption mappings Τᵤᵛ(·) as splines over a grid of length NCGSZ;
the inverses of Τᵤᵛ(·) as splines
=#
function compute_next(nos::Int64,g_step::Float64,gsz::Int64,sol::Array{Array{Float64,2},1},NCGSZ::Int64)
    local gg_step, cgrid, θnew, cnext, invcnext, gcmid
    gg_step=g_step
    cgrid=make_grid(gg_step,gsz)
    gcmid=(cgrid[end]+cgrid[1])/2.0
    θnew=[CubicSplineInterpolation(cgrid,sol[k][:,nos+1],extrapolation_bc = Interpolations.Line()) for k=1:nos];

    # the rows in cnext are attached to the current states (k0), while the columns are attached to the next period states (k1)
    cnext=[CubicSplineInterpolation(cgrid,sol[k0][:,k1],extrapolation_bc = Interpolations.Line()) for k0=1:nos, k1=1:nos];
    next_c_bounds=[(sol[k0][1,k1],sol[k0][end,k1]) for k0=1:nos, k1=1:nos];
    next_c_steps=[(next_c_bounds[k0,k1][2]-next_c_bounds[k0,k1][1])*(1/NCGSZ) for k0=1:nos, k1=1:nos];
    next_c_range_grid=[next_c_bounds[k0,k1][1]:next_c_steps[k0,k1]:(next_c_bounds[k0,k1][1]+NCGSZ*next_c_steps[k0,k1]) for k0=1:nos, k1=1:nos];
    inv_c_grid_vals=[[find_zero(x->cnext[k0,k1](x)-absc,gcmid,Order2()) for absc in next_c_range_grid[k0,k1]] for k0=1:nos, k1=1:nos];
    invcnext=[CubicSplineInterpolation(next_c_range_grid[k0,k1],inv_c_grid_vals[k0,k1],extrapolation_bc = Interpolations.Line()) for k0=1:nos, k1=1:nos];
    return θnew, cnext, invcnext
end;

#=
iterates (3.3) starting with the uniform distribution;
FGSZ is the length of the interpolation grid for the population distribution 
=#
function generate_dist_0(nos::Int64,upper_b::Float64,invCfun,CPROB::Array{Float64,2},BSP::Array{Float64,1},FGSZ::Int64)
    local F,GF,Fval
    GF=(0:FGSZ)*(upper_b/FGSZ);lgf=length(GF);
    Fval=[Float64((j-1)*(1/(lgf-1))) for k=1:nos, j=1:lgf];

    F_spline=[CubicSplineInterpolation(GF,Fval[k,:],extrapolation_bc = Interpolations.Line()) for k=1:nos];
    F=[x->(if x<GF[1] 0.0 elseif x>GF[end] 1.0 else F_spline[k](x) end) for k=1:nos];
    nextFval=[sum([(BSP[jj]*CPROB[jj,k]/(sum([BSP[j]*CPROB[j,k] for j=1:nos])))*F[jj](invCfun[jj,k](GF[i])) for jj=1:nos]) for k=1:nos, i=1:lgf];

    nextF_spline=[CubicSplineInterpolation(GF,nextFval[k,:],extrapolation_bc = Interpolations.Line()) for k=1:nos]
    nextF=[x->(if x<GF[1] 0.0 elseif x>GF[end] 1.0 else nextF_spline[k](x) end) for k=1:nos];

    while (maximum([maximum(abs.(Fval[k,:]-nextFval[k,:])) for k=1:nos]) > 1/(10^8))
        Fval=nextFval;
        F=nextF;
        nextFval=[sum([(BSP[jj]*CPROB[jj,k]/(sum([BSP[j]*CPROB[j,k] for j=1:nos])))*F[jj](invCfun[jj,k](GF[i])) for jj=1:nos]) for k=1:nos, i=1:lgf];
        nextF_spline=[CubicSplineInterpolation(GF,nextFval[k,:],extrapolation_bc = Interpolations.Line()) for k=1:nos];
        nextF=[x->(if x<GF[1] 0.0 elseif x>GF[end] 1.0 else nextF_spline[k](x) end) for k=1:nos];
    end;
    return F,upper_b,Fval
end;

#=
iterates (3.3) starting with supplied distribution (performs Step 4 in 3.2)
returns the distribution calculated in Step 4 as a list of spline;
the largest abscissa in the interpolation grid for those splines
the list of tabulated values for the distribution over the list of abscissas in the grid
=#
function generate_dist(nos::Int64,upper_b::Float64,fval::Array{Float64,2},invCfun,CPROB::Array{Float64,2},BSP::Array{Float64,1},FGSZ::Int64)
    local F,GF,Fval
    Fval=fval
    GF=(0:FGSZ)*(upper_b/FGSZ);lgf=length(GF);

    F_spline=[CubicSplineInterpolation(GF,Fval[k,:],extrapolation_bc = Interpolations.Line()) for k=1:nos];
    F=[x->(if x<GF[1] 0.0 elseif x>GF[end] 1.0 else F_spline[k](x) end) for k=1:nos];
    nextFval=[sum([(BSP[jj]*CPROB[jj,k]/(sum([BSP[j]*CPROB[j,k] for j=1:nos])))*F[jj](invCfun[jj,k](GF[i])) for jj=1:nos]) for k=1:nos, i=1:lgf];

    nextF_spline=[CubicSplineInterpolation(GF,nextFval[k,:],extrapolation_bc = Interpolations.Line()) for k=1:nos]
    nextF=[x->(if x<GF[1] 0.0 elseif x>GF[end] 1.0 else nextF_spline[k](x) end) for k=1:nos];

    while (maximum([maximum(abs.(Fval[k,:]-nextFval[k,:])) for k=1:nos]) > 1/(10^8))
        Fval=nextFval;
        F=nextF;
        nextFval=[sum([(BSP[jj]*CPROB[jj,k]/(sum([BSP[j]*CPROB[j,k] for j=1:nos])))*F[jj](invCfun[jj,k](GF[i])) for jj=1:nos]) for k=1:nos, i=1:lgf];
        nextF_spline=[CubicSplineInterpolation(GF,nextFval[k,:],extrapolation_bc = Interpolations.Line()) for k=1:nos];
        nextF=[x->(if x<GF[1] 0.0 elseif x>GF[end] 1.0 else nextF_spline[k](x) end) for k=1:nos];
    end;
    return F,upper_b,Fval
end;

# computes the left side of 3.2(a)
function market_clr(nos::Int64,upper_b::Float64,pfloFun,distFun,BSP::Array{Float64,1})
nmbp=40000;
    istep=upper_b/nmbp;
    igrid=(0:(nmbp-1))*istep;

    locI=zeros(nos)
        for k=1:nos
            for u in igrid
                locI[k]+=pfloFun[k](u+istep/2)*(distFun[k](u+istep)-distFun[k](u))
            end
        end
    return locI'*BSP
end;

function dist_support_zzz(nos::Int64,FF,upper_b::Float64)
    x00=0.1
    while maximum([FF[k](x00) for k=1:nos])<0.75
        x00+=0.05
    end
    return maximum([find_zero(x->FF[k](x)-1+1/10^9,x00,Order2()) for k=1:nos])
end;

function dist_support(nos::Int64,FF,upper_b::Float64)
    return maximum([find_zero(x->FF[k](x)-1+1/10^9,(0.01,upper_b+0.2)) for k=1:nos])
end;

#=
repeats Steps 1-5 in 3.2 until the market clearing 3.2(b) holds
INPUTS: total number of idiosyncratic states; step-size of the grid over c;
number of grid points; previous portfolios †θ; inverse mappings in the denominator of (3.2) constructed from †θ;
the most recent approximation of the spot price of the bond; the transition probability matrix;
the vector of paychecks for the employment classes; the list of steady-state probabilities;
(constant) aggregate income in the economy; risk-aversion parameter;
previous distribution of households †f as a list of values over the interpolation grid over c;
the total number of grid points in the *ranges* of the right side of (3.1)
as a function of the abscissas Τ (used for inversion);
the total number of grid points in the *ranges* of the future consumption mappings Τᵤᵛ(·)
(used for inversion); the total number of abscissas in the grid used to interpolate the distribution of households

OUTPUTS: market clearing for the current step (the left side of 3.2(b));
the distribution of the population as a list of spline objects;
the largest abscissa in the common interpolation grid for those splines (upper bound on consumption)
obtained from the distribution support;
the list of tabualted values for the distribution over the grid;
the list of portfolios *θᵤ(·) as splines; the future consumptions Tᵤᵛ(·) as splines; the inverses of Tᵤᵛ(·) as splines; 
the step-size in the new interpolation grid over consumption;
the new solution (future consumptions and present portfolio) as lists of tabulated
values over the consumption grid;
the accuracy with which the kernel condition holds (uniformly over the grid and the employment states);
the accuracy with which all first order conditions are satisfied;
the list of spot prices that have been tried with the last one being the most recent one;
the list of market clearing values (left side of 3.2(b)) that have been achieved with the last
one being the most recent one (the closest to 0)
=#
function main_prog(nos::Int64,g_step::Float64,gsz::Int64,theta,invh,spot::Float64,CPROB::Array{Float64,2},inc::Array{Float64,1},BSP::Array{Float64,1},ι::Float64,R::Float64,fval::Array{Float64,2},HRGSZ::Int64,NCGSZ::Int64,FGSZ::Int64)
    local gg_step, cgrid, clearing, θ, invH, SPOT, F, upper_bb, Fval, θnew, cnext, invcnext, sol, accu0, accu, ANSATZ, CLR, c_high, gcmid

    Fval=fval
    SPOT=spot
    invH=invh
    θ=theta
    gg_step=g_step
    cgrid=make_grid(gg_step,gsz)
    gcmid=(cgrid[end]+cgrid[1])/2.0
    no_iiter=0
    clr_check=0
    ANSATZ=[]
    CLR=[]
    while clr_check==0&&no_iiter<25
        no_iiter+=1;
        sol,accu0,accu=solve_sys(nos,gg_step,gsz,SPOT,θ,invH,CPROB,inc,ι,R)
        θnew, cnext, invcnext = compute_next(nos,gg_step,gsz,sol,NCGSZ)
        F,upper_bb,Fval = generate_dist(nos,cgrid[end],Fval,invcnext,CPROB,BSP,FGSZ)
        clearing=market_clr(nos,cgrid[end],θnew,F,BSP)
        if no_iiter==1
            ANSATZ=[SPOT]
            CLR=[clearing]
        else
            ANSATZ=vcat(ANSATZ,SPOT)
            CLR=vcat(CLR,clearing)
        end;
        if abs(CLR[end])<1/10^5 # adjust the choice for spot
            clr_check=1 #The market is cleared; end of step 4.
            elseif CLR[end]>=1/10^5&&length(ANSATZ)<=1
                SPOT+=length(ANSATZ)*0.0003;
            elseif CLR[end]<=-1/10^5&&length(ANSATZ)<=1
                SPOT-=length(ANSATZ)*0.0003;
        elseif length(ANSATZ)>=2
            #sq2err(x)=sum([(x[1]*ANSATZ[end-i]^2+x[2]*ANSATZ[end-i]+x[3]-CLR[end-i])^2 for i=0:2])
            #cfn=optimize(sq2err, zeros(3), BFGS()).minimizer
            #sq2err2(x)=(cfn[1]*x[1]^2+cfn[2]*x[1]+cfn[3])^2
            ##SPOT=(optimize(sq2err2, [SPOT],BFGS()).minimizer)[1]
            #SPOT=find_zero(sq2err2,ANSATZ[end])
            SPOT=(ANSATZ[end]*CLR[end-1]-ANSATZ[end-1]*CLR[end])/(CLR[end-1]-CLR[end])
        end;

        #because SPOT changes invH has to be recomputed -- with the old theta

        H=[CubicSplineInterpolation(cgrid,[x+θ[k](x)*SPOT for x in cgrid],extrapolation_bc = Interpolations.Line()) for k=1:nos];
        H_range_bounds=[(cgrid[1]+θ[k](cgrid[1])*SPOT,cgrid[end]+θ[k](cgrid[end])*SPOT) for k=1:nos];
        H_ran_step=[((H_range_bounds[k][2]-H_range_bounds[k][1])/HRGSZ) for k=1:nos]
        H_range_grids=[H_range_bounds[k][1]:H_ran_step[k]:(H_range_bounds[k][1]+HRGSZ*H_ran_step[k]) for k=1:nos];
        H_inv_vals=[[find_zero(x->H[k](x)-absc,gcmid,Order2()) for absc=H_range_grids[k]] for k=1:nos];
        invH=[CubicSplineInterpolation(H_range_grids[k],H_inv_vals[k],extrapolation_bc = Interpolations.Line()) for k=1:nos];

        #LowBound=maximum([(0.0+θ[k](0.0)*SPOT-inc[k])/ι for k=1:nos])+0.0001

        c_high=dist_support(nos,F,upper_bb)+0.3
        gg_step=(c_high/gsz)
        cgrid=make_grid(gg_step,gsz) #the grid changes because SPOT is new
        gcmid=(cgrid[end]+cgrid[1])/2.0
    end
    return clearing, F, upper_bb, Fval, θnew, cnext, invcnext, gg_step, sol, accu0, accu, ANSATZ, CLR
end;

#=
same as main_prog() except that it does not take previous portfolios †θ on the input
(generates those internally as an initial ansatz choice)
=#
function main_prog_0(nos::Int64,g_step::Float64,gsz::Int64,theta,invh,spot::Float64,CPROB::Array{Float64,2},inc::Array{Float64,1},BSP::Array{Float64,1},ι::Float64,R::Float64,HRGSZ::Int64,NCGSZ::Int64,FGSZ::Int64)
    local gg_step, cgrid, clearing, θ, invH, SPOT, F, upper_bb, Fval, θnew, cnext, invcnext, sol, accu0, accu, ANSATZ, CLR, c_high, gcmid

    SPOT=spot
    invH=invh
    θ=theta
    gg_step=g_step
    cgrid=make_grid(gg_step,gsz)
    gcmid=(cgrid[end]+cgrid[1])/2.0
    no_iiter=0
    clr_check=0
    ANSATZ=[]
    CLR=[]
    while clr_check==0&&no_iiter<25
        no_iiter+=1;
        sol,accu0,accu=solve_sys(nos,gg_step,gsz,SPOT,θ,invH,CPROB,inc,ι,R)
        θnew, cnext, invcnext = compute_next(nos,gg_step,gsz,sol,NCGSZ)
        F,upper_bb,Fval = generate_dist_0(nos,cgrid[end],invcnext,CPROB,BSP,FGSZ)
        clearing=market_clr(nos,cgrid[end],θnew,F,BSP)
        if no_iiter==1
            ANSATZ=[SPOT]
            CLR=[clearing]
        else
            ANSATZ=vcat(ANSATZ,SPOT)
            CLR=vcat(CLR,clearing)
        end;
        if abs(CLR[end])<1/10^5 # adjust the choice for spot
            clr_check=1 #The market is cleared; end of step 4.
            elseif CLR[end]>=1/10^5&&length(ANSATZ)<=1
                SPOT+=length(ANSATZ)*0.0003;
            elseif CLR[end]<=-1/10^5&&length(ANSATZ)<=1
                SPOT-=length(ANSATZ)*0.0003;
        elseif length(ANSATZ)>=2
             SPOT=(ANSATZ[end]*CLR[end-1]-ANSATZ[end-1]*CLR[end])/(CLR[end-1]-CLR[end])
        end;

        #because SPOT changes invH has to be recomputed -- with the old theta

        H=[CubicSplineInterpolation(cgrid,[x+θ[k](x)*SPOT for x in cgrid],extrapolation_bc = Interpolations.Line()) for k=1:nos];
        H_range_bounds=[(cgrid[1]+θ[k](cgrid[1])*SPOT,cgrid[end]+θ[k](cgrid[end])*SPOT) for k=1:nos];
        H_ran_step=[((H_range_bounds[k][2]-H_range_bounds[k][1])/HRGSZ) for k=1:nos]
        H_range_grids=[H_range_bounds[k][1]:H_ran_step[k]:(H_range_bounds[k][1]+HRGSZ*H_ran_step[k]) for k=1:nos];
        H_inv_vals=[[find_zero(x->H[k](x)-absc,gcmid,Order2()) for absc=H_range_grids[k]] for k=1:nos];
        invH=[CubicSplineInterpolation(H_range_grids[k],H_inv_vals[k],extrapolation_bc = Interpolations.Line()) for k=1:nos];

        #LowBound=maximum([(0.0+θ[k](0.0)*SPOT-inc[k])/ι for k=1:nos])+0.0001

        c_high=dist_support(nos,F,upper_bb)+0.3
        gg_step=(c_high/gsz)
        cgrid=make_grid(gg_step,gsz) #the grid changes because SPOT is new
        gcmid=(cgrid[end]+cgrid[1])/2.0
    end
    return clearing, F, upper_bb, Fval, θnew, cnext, invcnext, gg_step, sol, accu0, accu, ANSATZ, CLR
end;

#=
runs Steps 1-6 in 3.2 for a specified number of iterations
after each iteration prints the iteration number followed by the pair of the second maximum in 3.2(c)
and the largest of the two maxima in 3.2(c)

INPUTS: number of employment states; starting iteration number;
last iteration number (after which the program stops); desired largest value of the two
maxima in 3.2(c) (the program stops is this threshold is acheived);
the step-size of the grid over the consumption range;
the number of grid points over the consumption range;
the most recently accepted portfolios (as splines over the consumption range) †θ;
the inverse mappings in the denominator of (3.2) as spline objects;
the most recently accepted spot price; the transition probability matrix;
the list of paychecks for the various employment classes;
list of steady-state probabilities; aggregate income in the economy;
risk aversion; the total number of grid points in the *ranges* of the right side of (3.1) as a function
of the abscissas ν (used for inversion);
the total number of grid points in the *ranges* of the future consumption mappings Τᵤᵛ(·) (used for inversion);
the total number of abscissas in the grid used to interpolate the distribution of households 

OUTPUTS: the number of the last iteration; the largest of the two maxima in 3.2(c);
the second maximum in 3.2(c); the last market clearing (left side of 3.2(b));
the upper bound on consumption (extracted from the support of the distribution);
the distribution of all households as a list of spline objects;
the distribution of households as a list of interpolated values over the abscissas;
the list of portfolios as spline objects;
the inverses of the right sides of (3.1) as functions of Τ (espressed as splines);
the most recently computed spot price;
the future consumptions Τᵤᵛ(·) as splines; the inverses of Τᵤᵛ(·) as splines;
the step-size of the last interpolation grid over consumption;
the number of grid points over the consumption range (same as input);
the actual solution (future consumptions and present portfolios) attached to states of employment
                                  and abscissas in the grid over the consumption range;
the largest error with which the last solution satisfies the kernel
                       condition (3.2) for all employment states and abscissas on the consumption grid;
the largest error with which the last solution satisfies all
                       first order conditions (3.1) and (3.2) for all states and abscissas;
the list of spot prices that have been tried during the last iteration so that SPOT=ANSATZ[end];
the list of market clearing values (left side of 3.2(a)) that have been achieved during the last iteration
with the last one being the most recent one (the closest to 0)
=#
function find_equil(nos::Int64,start_iter::Int64,nmb_of_iter::Int64,conv_tol::Float64,g_step::Float64,gsz::Int64,theta,invh,spot::Float64,CPROB::Array{Float64,2},inc::Array{Float64,1},BSP::Array{Float64,1},ι::Float64,R::Float64,HRGSZ::Int64,NCGSZ::Int64,FGSZ::Int64)
    local gg_step, cgrid, no_iter, solPrev, invH, SPOT, θ, SPOTprev, accuPrev, accu0Prev, θprev, conv_check, conv_check_c, clearing, F, upper_bb, Fval, θnew, cnext, invcnext, sol, accu0, accu, ANSATZ, CLR, c_high, gcmid

    SPOT=spot
    invH=invh
    θ=theta
    gg_step=g_step
    cgrid=make_grid(gg_step,gsz)
    gcmid=(cgrid[end]+cgrid[1])/2.0
    no_iter=start_iter;
    conv_check=100.0
    conv_check_c=100.0
    no_iter+=1;
    println(no_iter)
    clearing, F, upper_bb, Fval, θnew, cnext, invcnext, gg_step, sol, accu0, accu, ANSATZ, CLR = main_prog_0(nos,gg_step,gsz,θ,invH,SPOT,CPROB,inc,BSP,ι,R,HRGSZ,NCGSZ,FGSZ);
    SPOTprev=SPOT;
    SPOT=ANSATZ[end];
    accuPrev=accu;
    accu0Prev=accu0;
    solPrev=sol;
    θprev=θ;
    θ=θnew;
    H=[CubicSplineInterpolation(cgrid,[x+θ[k](x)*SPOT for x in cgrid],extrapolation_bc = Interpolations.Line()) for k=1:nos];
    H_range_bounds=[(cgrid[1]+θ[k](cgrid[1])*SPOT,cgrid[end]+θ[k](cgrid[end])*SPOT) for k=1:nos];
    H_ran_step=[((H_range_bounds[k][2]-H_range_bounds[k][1])/HRGSZ) for k=1:nos]
    H_range_grids=[H_range_bounds[k][1]:H_ran_step[k]:(H_range_bounds[k][1]+HRGSZ*H_ran_step[k]) for k=1:nos];
    H_inv_vals=[[find_zero(x->H[k](x)-absc,gcmid,Order2()) for absc=H_range_grids[k]] for k=1:nos];
    invH=[CubicSplineInterpolation(H_range_grids[k],H_inv_vals[k],extrapolation_bc = Interpolations.Line()) for k=1:nos];
    c_high=dist_support(nos,F,upper_bb)+0.3
    gg_step=(c_high/gsz)
    cgrid=make_grid(gg_step,gsz)
    gcmid=(cgrid[end]+cgrid[1])/2.0
    while no_iter<nmb_of_iter&&conv_check>conv_tol
        no_iter+=1;
        println(no_iter)
        clearing, F, upper_bb, Fval, θnew, cnext, invcnext, gg_step, sol, accu0, accu, ANSATZ, CLR = main_prog(nos,gg_step,gsz,θ,invH,SPOT,CPROB,inc,BSP,ι,R,Fval,HRGSZ,NCGSZ,FGSZ);
        cgrid=make_grid(gg_step,gsz)
        gcmid=(cgrid[end]+cgrid[1])/2.0
        conv_check=maximum([maximum([maximum(abs.(sol[k][:,i]-solPrev[k][:,i])) for i=1:(nos+1)]) for k=1:nos]);
        conv_check_c=maximum([maximum([maximum(abs.(sol[k][:,i]-solPrev[k][:,i])) for i=1:nos]) for k=1:nos]);
        println((conv_check_c,conv_check))
        SPOTprev=SPOT;
        SPOT=ANSATZ[end];
        accuPrev=accu;
        accu0Prev=accu0;
        solPrev=sol;
        θprev=θ;
        θ=θnew;
        H=[CubicSplineInterpolation(cgrid,[x+θ[k](x)*SPOT for x in cgrid],extrapolation_bc = Interpolations.Line()) for k=1:nos];
        H_range_bounds=[(cgrid[1]+θ[k](cgrid[1])*SPOT,cgrid[end]+θ[k](cgrid[end])*SPOT) for k=1:nos];
        H_ran_step=[((H_range_bounds[k][2]-H_range_bounds[k][1])/HRGSZ) for k=1:nos]
        H_range_grids=[H_range_bounds[k][1]:H_ran_step[k]:(H_range_bounds[k][1]+HRGSZ*H_ran_step[k]) for k=1:nos];
        H_inv_vals=[[find_zero(x->H[k](x)-absc,gcmid,Order2()) for absc=H_range_grids[k]] for k=1:nos];
        invH=[CubicSplineInterpolation(H_range_grids[k],H_inv_vals[k],extrapolation_bc = Interpolations.Line()) for k=1:nos];
        c_high=dist_support(nos,F,upper_bb)+0.3
        gg_step=(c_high/gsz)
        cgrid=make_grid(gg_step,gsz)
        gcmid=(cgrid[end]+cgrid[1])/2.0
    end
    return no_iter, conv_check, conv_check_c, clearing, upper_bb, F, Fval, θ, invH, SPOT, cnext, invcnext, gg_step, gsz, sol, accu0, accu, ANSATZ, CLR
end;
