######################################################################################################
#
# Julia code
#
# Functions to implement the method described in Sec. 18.7 in "Recursive Macroeconomic Theory" (RMT)
#                                     by Lars Ljungqvist and Thomas Sargent
#
# The code included here is an adaptation from the MATLAB code that accompanies RMT.
#
# This code supplements the paper "Self-Consistent Transport in Heterogeneous Agent Models" [STHAM]
#                                                                     by Andrew Lyasoff
# 
# Copyright © 2019-2025 Andrew Lyasoff <mathema@laysoff.net>
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
the hourly wage; the discount factor; the grid over the range of wealth for each state;
initial probability distribution over the grid (the starting point for the iterations leading to
the long-run distribution of the population); the number of grid points over the range of wealth (assumed to be
the same for all idiosyncratic states; the transition probability matrix; the risk-aversion parameter

#OUTPUTS: the long run distribution of the population; the optimal policy rule
=#
function LS_method(intr::Float64,accu_pow1::Int64,accu_pow2::Int64,hours::Array{Float64,1},wage::Float64,β::Float64,grda::Array{Float64,2},grdpd::Array{Float64,2},CPROB::Array{Float64,2},R::Float64)
#    local TBL,TBL1,TBL2,tbl,tbl1,tbl2,Λ0,Λ1
    local TBL::Matrix{Tuple{Float64, Int64}},TBL1::Matrix{Float64},TBL2::Matrix{Int64},tbl::Matrix{Tuple{Float64, Int64}},tbl1::Matrix{Float64},tbl2::Matrix{Int64},Λ0::Matrix{Float64},Λ1::Matrix{Float64},nos::Int64,Lg::Int64
    #
    nos=(size(grda)[1]); # first dimension of grda is the number of states
    Lg=(size(grda)[2]); # the second dimension of grda is the number of grid points over wealth
    if !(nos==size(grdpd)[1]) # the grids used to track wealth and probability distribution must coincide
        println("incompatable arrays")
    end;
    if !(Lg==size(grdpd)[2]) # the grids used to track wealth and probability distribution must coincide
        println("incompatable arrays")
    end;
    #
    # table of local utilities from (consumption = income - investment)
    # k is the current employment state
    # i labels the previous investment ; j labels the current investment
    currentU=[U((1+intr)*grda[k,i]+wage*hours[k]-grda[k,j],R) for k=1:nos, i=1:Lg, j=1:Lg];
    #
    # the value function for the next period
    # initial value = utility from consumption in the last period
    nextV=[U((1+intr)*grda[k,j]+wage*hours[k],R) for k=1:nos, j=1:Lg];
    #
    # period T-1 present value in every employment state k
    #       for every choice of investment i
    # N.B. findmax() returns (value, position_in_the_array_of_choices)
    # N.B. postion_in_the_array_of_choices = optimal policy rule
    TBL=[findmax(currentU[k,i,:]+β*(CPROB[k,:]'*nextV)') for k=1:nos, i=1:Lg];
    TBL1=first.(TBL);TBL2=last.(TBL);
    #
    # period T-2 present value in every employment state k
    tbl=[findmax(currentU[k,i,:]+β*(CPROB[k,:]'*TBL1)') for k=1:nos, i=1:Lg];
    tbl1=first.(tbl);tbl2=last.(tbl);
    #
    # iterate the Bellman equation to convergence
    # N.B. tbl2 gives the optimal choice on the grid
    #         for every emplyment state k and every previus choice on the grid i
    # N.B. unlike the original code in RMT, the stopping rule
    #          involves both the value function and the policy rule
    while maximum(abs.(TBL2-tbl2))>0||maximum(abs.(TBL1-tbl1))>1/(10^accu_pow1)
        TBL1=tbl1;
        TBL2=tbl2;
        tbl=[findmax(currentU[k,i,:]+β*(CPROB[k,:]'*TBL1)') for k=1:nos, i=1:Lg];
        tbl1=first.(tbl);tbl2=last.(tbl);
    end
    #
    # initial probability distrubution (supplied as a function argument)
    Λ0=grdpd;
    #
    # iterates the probability distribution to convergence -- see equation (1.1) is STHAM
    # for a given interest and the policy function attached to it
    Λ1=[sum([Λ0[k,i]*CPROB[k,kk]*II(ii,k,i,tbl2) for k=1:nos, i=1:Lg]) for kk=1:nos, ii=1:Lg];
    #
    Λ0=Λ1;
    Λ1=[sum([Λ0[k,i]*CPROB[k,kk]*II(ii,k,i,tbl2) for k=1:nos, i=1:Lg]) for kk=1:nos, ii=1:Lg];
    while maximum(abs.(Λ1-Λ0))>1/10^accu_pow2
        Λ0=Λ1;
        Λ1=[sum([Λ0[k,i]*CPROB[k,kk]*II(ii,k,i,tbl2) for k=1:nos, i=1:Lg]) for kk=1:nos, ii=1:Lg];
    end
    #
    # return the probability distribution at T = ∞ , followed by the policy rule
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
    #
    # N.B. the lower bound on the wealth is determined as in RMT
    #
    for iii=1:min(length(init),no_iter) # pass through the list of ansatz choices for r
        intr=init[iii] # pick an interest rate from the list
        println(intr)
        minkap=-b;
        if intr<=0.0
            minkap = -b;
        else
            minkap=-min(b, wage*hours[1]/intr);
        end;
        #
        # construct uniform grid on the range of wealth
        # N.B. maxkap is a function argument
        grda=[minkap+i*(maxkap-minkap)/(Lg-1) for i=0:(Lg-1)];
        #
        # define the same grid for every employment category
        grdaa=[grda[i] for k=1:nos, i=1:length(grda)];
        #
        # start with the uniform distribution on the grid
        ini_dist=[Float64(1/(nos*Lg)) for k=1:nos, i=1:Lg];
        #
        # iterate the policy function and the distribution 
        Λ1,tbl2=LS_method(intr,accu_pow1,accu_pow2,hours,wage,β,grdaa,ini_dist,CPROB,R);
        #
        # test whether Λ1 has the properties of a "distribution"
        # e.g. the "weight" of every employment state must be the same as the
        # steady-state probability of that state
        checkpoints=vcat(checkpoints,[[sum(Λ1)-1.0,[sum(Λ1[i,:]) for i=1:nos]-BSP]]);
        bond_price=ι/(1+intr)
        #
        # optimal consumption levels and investment for every employment category k
        # and every previous choice of investment i
        cons=[(1+intr)*grda[i]+wage*hours[k]-grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
        aassets=[grda[Int64(tbl2[k,i])] for k=1:nos, i=1:Lg];
        shares=[grda[Int64(tbl2[k,i])]/bond_price for k=1:nos, i=1:Lg];
        #
        # For a given distribution of the population construct
        # the conditional discrete cummulative distribution function
        # for every category of employment
        SL=[sum(Λ1[state,:]) for state=1:nos];
        LF=[[Λ1[state,1]/SL[state]] for state=1:nos];
        for i=2:Lg
            for state=1:nos
                LF[state]=vcat(LF[state],LF[state][end]+Λ1[state,i]/SL[state])
            end
        end
        demand=[grdaa[kk,Int64(tbl2[kk,ii])] for kk=1:nos, ii=1:Lg];
        #
        # calculate the expected demand (for the chosen intr)
        Ed=sum(Λ1.*demand);
        println("     ",Ed)
        #
        # update the latest output from LS_method() the yields (+) expected demand
        # OR update the latest output from LS_method() the yields (--) expected demand
        if Ed>=0.0
            Λ1p=Λ1;
            tbl2p=tbl2;
        else
            Λ1m=Λ1;
            tbl2m=tbl2;
        end
        #
        # update the list of expected demands
        demands=vcat(demands,Ed)
    end
    #
    # test if there is at least one ++ demand AND at least one -- demand
    if length(findall(x->x<0,demands))*length(findall(x->x>0,demands))>0
        if length(init)<no_iter # is the list of ansatz choices smaller than the one requested
            intr_trials=init;
            for iii=(length(init)+1):no_iter # produce additional r's to be tried
                #
                # the next r is the average of the last one that gave positive expected demand
                # and the last one that gave negative expected demand
                intr=(intr_trials[findlast(x->x>0,demands)]+intr_trials[findlast(x->x<0,demands)])/2;
                println(intr)
                intr_trials=vcat(intr_trials,intr); # update the list of r's that are being tried
                #
                # the avove program is now repated with the NEW value for intr (next r to try)
                # N.B. new value for intr is generated (as an average) for as many times
                # as no_iter demands
                #
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
        println("bad initial choice") # either all demans are ++ OR all are --
        Λ1p_ret=if (@isdefined Λ1p) Λ1p else Array{Float64}(undef, 0, 0) end
        tbl2p_ret=if (@isdefined tbl2p) tbl2p else Array{Int64}(undef, 0, 0) end
        Λ1m_ret=if (@isdefined Λ1m) Λ1m else Array{Float64}(undef, 0, 0) end
        tbl2m_ret=if (@isdefined tbl2m) tbl2m else Array{Int64}(undef, 0, 0) end
        return init, demands, checkpoints, Λ1p_ret, tbl2p_ret, Λ1m_ret, tbl2m_ret
    end;
end;
