######################################################################################################
#
# Julia code
#
# Functions to implement the method described in Sec. 18.7 in "Recursive Macroeconomic Theory" (RMT)
#                                     by Lars Ljungqvist and Thomas Sargent
#
# The code included here is an adaptation from the MATLAB code that accompanies RMT.
#
# This code supplements the paper "Another look at the distribution of income and wealth in the Macroeconomy" (DIW)
#                                        by Andrew Lyasoff (www.andrewlyasoff.tech)
# 
# Copyright © 2019-2023 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################

#=
The following are assigned from   ini-setup-RMT-ch18.jl

CPROB, nos, BSP, hours, wage, ι, ρ, β, R
=#



#=
Utility function
=#
function U(ξ::Float64,R::Float64)
    if ξ <= 0.
        -9e10
        else ((1/ξ^(R-1.0))/(1.0 - R))
    end
end;

#=
Indicator function
=#
function II(aa::Int64,k::Int64,a::Int64,tbl2::Array{Int64,2})
    if aa==tbl2[k,a]
        return 1
    else
        return 0
    end
end

#=
The next function calculates the average assets as a function of r

NB: In the RMT MATLAB program the iterations stop after the first repetition of the policy.
    Here the iterations stop after the first repetition of the policy AND the value -- see README.md .

INPUTS: interest rate to try; convergence tolerance for the value function as an integer power of 10;
convergence tolerance for the long-run distribution of the population; the number of states; the list of work hours;
the hourly wage; the discount factor; the grid on the range of wealth for each state;
initial probability distribution on the grid (the starting point for the iterations leading to
the long-run distribution of the population); the number of grid points on the range of wealth (assumed to be
the same for all idiosyncratic states; the transition probability matrix; the risk-aversion parameter

#OUTPUTS: the long run distribution of the population; the optimal policy rule
=#
function LS_method(intr::Float64,accu_pow1::Int64,accu_pow2::Int64,hours::Array{Float64,1},wage::Float64,β::Float64,grda::Array{Float64,2},grdpd::Array{Float64,2},CPROB::Array{Float64,2},R::Float64)
#    local TBL,TBL1,TBL2,tbl,tbl1,tbl2,Λ0,Λ1
    local TBL::Matrix{Tuple{Float64, Int64}},TBL1::Matrix{Float64},TBL2::Matrix{Int64},tbl::Matrix{Tuple{Float64, Int64}},tbl1::Matrix{Float64},tbl2::Matrix{Int64},Λ0::Matrix{Float64},Λ1::Matrix{Float64},nos::Int64,Lg::Int64
    #
    nos=(size(grda)[1]);
    Lg=(size(grda)[2]);
    if !(nos==size(grdpd)[1])
        println("incompatable arrays")
    end;
    if !(Lg==size(grdpd)[2])
        println("incompatable arrays")
    end;
    currentU=[U((1+intr)*grda[k,i]+wage*hours[k]-grda[k,j],R) for k=1:nos, i=1:Lg, j=1:Lg];
    nextV=[U((1+intr)*grda[k,j]+wage*hours[k],R) for k=1:nos, j=1:Lg];
#
    TBL=[findmax(currentU[k,i,:]+β*(CPROB[k,:]'*nextV)') for k=1:nos, i=1:Lg];
    TBL1=first.(TBL);TBL2=last.(TBL);
#
    tbl=[findmax(currentU[k,i,:]+β*(CPROB[k,:]'*TBL1)') for k=1:nos, i=1:Lg];
    tbl1=first.(tbl);tbl2=last.(tbl);
#
    while maximum(abs.(TBL2-tbl2))>0||maximum(abs.(TBL1-tbl1))>1/(10^accu_pow1)
        TBL1=tbl1;
        TBL2=tbl2;
        tbl=[findmax(currentU[k,i,:]+β*(CPROB[k,:]'*TBL1)') for k=1:nos, i=1:Lg];
        tbl1=first.(tbl);tbl2=last.(tbl);
    end
#
    Λ0=grdpd;
#
    Λ1=[sum([Λ0[k,i]*CPROB[k,kk]*II(ii,k,i,tbl2) for k=1:nos, i=1:Lg]) for kk=1:nos, ii=1:Lg];
#
    Λ0=Λ1;
    Λ1=[sum([Λ0[k,i]*CPROB[k,kk]*II(ii,k,i,tbl2) for k=1:nos, i=1:Lg]) for kk=1:nos, ii=1:Lg];
    while maximum(abs.(Λ1-Λ0))>1/10^accu_pow2
        Λ0=Λ1;
        Λ1=[sum([Λ0[k,i]*CPROB[k,kk]*II(ii,k,i,tbl2) for k=1:nos, i=1:Lg]) for kk=1:nos, ii=1:Lg];
    end
    return Λ1,tbl2
end


#=
The next function adjusts the interest rate to achieve market clearing. Every new rate that is tried is the arithmetic
average of the last rate that yields strictly positive demand and the last rate that yields strictly negative demand.
The input variable 'init' is the list of initial ansatz rates to try (must be of length ≥ 2). If the aggregate
demands associated with the rates supplied by 'init' are all positive (all negative) the programs prints the message
"bad initial choice" and returns the associated aggregate demands. 

INPUTS: 'no_iter' is the total number of interest rates to try; 'ι' is the aggregate income in the economy;
'b' and 'maxkap' are parameters that determine the range of private wealth to be covered by the grid (both parameters
have the same meaning as in RMT); BSP is the list of steady state probabilities associated
with CPROB -- All other inputs are as in the previous function

OUTPUTS: list if all interest rates that have been tried; list of all associated aggregate demands;
list of all checkpoints (whether or not the calculated limiting distribution is numerically acceptable) --
must return numbers that are numerical zeros; distribution with the most recent rate that yields a positive demand; optimal policy with the most recent rate that yields positive demand; distribution with the most recent rate that yields negative demand; optimal policy with the most recent rate that yields negative demand;
=#


function iterLS(no_iter::Int64,accu_pow1::Int64,accu_pow2::Int64,nos::Int64,hours::Array{Float64,1},wage::Float64,β::Float64,Lg::Int64,CPROB::Array{Float64,2},R::Float64,ι::Float64,init::Vector{Float64},b::Float64,maxkap::Float64,BSP::Vector{Float64})
    local test::Bool, ini_check::Bool, ini_dist::Matrix{Float64}, intr::Float64, minkap::Float64, grda::Vector{Float64}, grdaa::Matrix{Float64}, Λ1::Matrix{Float64}, tbl2::Matrix{Int64},  Λ1p::Matrix{Float64}, tbl2p::Matrix{Int64},Λ1m::Matrix{Float64}, tbl2m::Matrix{Int64}, bond_price::Float64, cons::Matrix{Float64}, aassets::Matrix{Float64}, shares::Matrix{Float64}, SL::Vector{Float64}, LF::Vector{Vector{Float64}}, demand::Matrix{Float64}, intr_trials::Vector{Float64}, Ed::Float64,last_pos::Float64,last_neg::Float64,d_pos::Float64,d_neg::Float64,w_pos::Float64,w_neg::Float64
    checkpoints=[]
    demands=Float64[];
    ini_check=false
    test=true
    for iii=1:min(length(init),no_iter)
        intr=init[iii]
        println(intr)
        minkap=-b;
        if intr<=0.0
            minkap = -b;
        else
            minkap=-min(b, wage*hours[1]/intr);
        end;
        grda=[minkap+i*(maxkap-minkap)/(Lg-1) for i=0:(Lg-1)];
        grdaa=[grda[i] for k=1:nos, i=1:length(grda)];
        ini_dist=[Float64(1/(nos*Lg)) for k=1:nos, i=1:Lg];
        Λ1,tbl2=LS_method(intr,accu_pow1,accu_pow2,hours,wage,β,grdaa,ini_dist,CPROB,R);
        checkpoints=vcat(checkpoints,[[sum(Λ1)-1.0,[sum(Λ1[i,:]) for i=1:nos]-BSP]]);
        bond_price=ι/(1+intr)
        cons=[(1+intr)*grda[i]+wage*hours[k]-grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
        aassets=[grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
        shares=[grda[Int64(tbl2[k,i])]/bond_price for k=1:nos, i=1:Lg];
        SL=[sum(Λ1[state,:]) for state=1:nos];
        LF=[[Λ1[state,1]/SL[state]] for state=1:nos];
        for i=2:Lg
            for state=1:nos
                LF[state]=vcat(LF[state],LF[state][end]+Λ1[state,i]/SL[state])
            end
        end
        demand=[grdaa[kk,Int64(tbl2[kk,ii])] for kk=1:nos, ii=1:Lg];
        Ed=sum(Λ1.*demand);
        println("     ",Ed)
        if Ed>=0.0
            Λ1p=Λ1;
            tbl2p=tbl2;
        else
            Λ1m=Λ1;
            tbl2m=tbl2;
        end
        demands=vcat(demands,Ed)
    end
    if length(findall(x->x<0,demands))*length(findall(x->x>0,demands))>0
        if length(init)<no_iter
            intr_trials=init;
            for iii=(length(init)+1):no_iter
                intr=(intr_trials[findlast(x->x>0,demands)]+intr_trials[findlast(x->x<0,demands)])/2;
                #=
                last_pos=intr_trials[findlast(x->x>0,demands)];
                last_neg=intr_trials[findlast(x->x<0,demands)];
                d_pos=abs(demands[findlast(x->x>0,demands)]);
                d_neg=abs(demands[findlast(x->x<0,demands)]);
                if d_pos>d_neg
                    w_pos=d_neg/(d_pos+d_neg);
                    w_neg=d_pos/(d_pos+d_neg);
                else
                    w_pos=d_pos/(d_pos+d_neg);
                    w_neg=d_neg/(d_pos+d_neg);
                end;
                intr=w_pos*last_pos+w_neg*last_neg;
                =#
                println(intr)
                intr_trials=vcat(intr_trials,intr);
                minkap=-b;
                if intr<=0.0
                    minkap = -b;
                else
                    minkap=-min(b, wage*hours[1]/intr);
                end;
                grda=[minkap+i*(maxkap-minkap)/(Lg-1) for i=0:(Lg-1)];
                grdaa=[grda[i] for k=1:nos, i=1:length(grda)];
                ini_dist=[Float64(1/(nos*Lg)) for k=1:nos, i=1:Lg];
                Λ1,tbl2=LS_method(intr,accu_pow1,accu_pow2,hours,wage,β,grdaa,ini_dist,CPROB,R);
                checkpoints=vcat(checkpoints,[[sum(Λ1)-1.0,[sum(Λ1[i,:]) for i=1:nos]-BSP]]);
                bond_price=ι/(1+intr)
                cons=[(1+intr)*grda[i]+wage*hours[k]-grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
                aassets=[grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
                shares=[grda[Int64(tbl2[k,i])]/bond_price for k=1:nos, i=1:Lg];
                SL=[sum(Λ1[state,:]) for state=1:nos];
                LF=[[Λ1[state,1]/SL[state]] for state=1:nos];
                for i=2:Lg
                    for state=1:nos
                        LF[state]=vcat(LF[state],LF[state][end]+Λ1[state,i]/SL[state])
                    end
                end
                demand=[grdaa[kk,Int64(tbl2[kk,ii])] for kk=1:nos, ii=1:Lg];
                Ed=sum(Λ1.*demand);
                println("     ",Ed)
                if Ed>=0.0
                    Λ1p=Λ1;
                    tbl2p=tbl2;
                else
                    Λ1m=Λ1;
                    tbl2m=tbl2;
                end
                demands=vcat(demands,Ed);
            end
        end
        Λ1p_ret=if (@isdefined Λ1p) Λ1p else Array{Float64}(undef, 0, 0) end
        tbl2p_ret=if (@isdefined tbl2p) tbl2p else Array{Int64}(undef, 0, 0) end
        Λ1m_ret=if (@isdefined Λ1m) Λ1m else Array{Float64}(undef, 0, 0) end
        tbl2m_ret=if (@isdefined tbl2m) tbl2m else Array{Int64}(undef, 0, 0) end
        return intr_trials, demands, checkpoints, Λ1p_ret, tbl2p_ret, Λ1m_ret, tbl2m_ret
    else
        println("bad initial choice")
        Λ1p_ret=if (@isdefined Λ1p) Λ1p else Array{Float64}(undef, 0, 0) end
        tbl2p_ret=if (@isdefined tbl2p) tbl2p else Array{Int64}(undef, 0, 0) end
        Λ1m_ret=if (@isdefined Λ1m) Λ1m else Array{Float64}(undef, 0, 0) end
        tbl2m_ret=if (@isdefined tbl2m) tbl2m else Array{Int64}(undef, 0, 0) end
        return init, demands, checkpoints, Λ1p_ret, tbl2p_ret, Λ1m_ret, tbl2m_ret
    end;
end;
