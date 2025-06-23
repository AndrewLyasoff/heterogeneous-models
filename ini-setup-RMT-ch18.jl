######################################################################################################
#
# Julia code
#
# Sets the parameter values from Sec. 18.7 in "Recursive Macroeconomic Theory" (RMT)
#                                     by Lars Ljungqvist and Thomas Sargent
#
# This code supplements the paper "The Time-Interlaced Self-Consistent Master System of Heterogeneous-Agent Models" [STHAM]
#                                        by Andrew Lyasoff 
# 
# Copyright © 2019-2025 Andrew Lyasoff <mathema@lyasoff.net>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################

#=
This code assigns values to:

CPROB, nos, BSP, hours, wage, ι, ρ, β, R
=#

#=
The MATLAB output from [prob,logs,invdist,alambda,asigma]=markovappr(rho,sigma,3,N);

with markovappr() loaded from markovapprox.m (see RMT) with parameters

N        = 7;             % number of discretized states
rho      = 0.2;           % first-order autoregressive coefficient
sigmaint = 0.4;           % intermediate value to calculate sigma
sigma    = sigmaint*sqrt(1-rho^2); % standard deviation of error_t

logs =

 -1.2000   -0.8000   -0.4000         0    0.4000    0.8000    1.2000

prob =

    0.0262    0.1529    0.3615    0.3286    0.1147    0.0153    0.0008
    0.0160    0.1147    0.3286    0.3615    0.1529    0.0247    0.0015
    0.0095    0.0828    0.2874    0.3828    0.1961    0.0384    0.0029
    0.0054    0.0575    0.2420    0.3902    0.2420    0.0575    0.0054
    0.0029    0.0384    0.1961    0.3828    0.2874    0.0828    0.0095
    0.0015    0.0247    0.1529    0.3615    0.3286    0.1147    0.0160
    0.0008    0.0153    0.1147    0.3286    0.3615    0.1529    0.0262

=#

#=
Thransition probability matrix for the idiosyncratic state
=#
CPROB=[0.0262397497796233 0.152923483594817 0.361483063911416 0.328567584707170 0.114741751017859 0.0152660804766739 0.000778286512441828; 0.0160443669891157 0.114741751017859 0.328567584707170 0.361483063911415 0.152923483594817 0.0247005562212874 0.00153919355833587; 0.00945177150556068 0.0828345049143809 0.287445156270025 0.382789297365280 0.196113758787202 0.0384369616149092 0.00292854954264210; 0.00536221880153007 0.0575309935179400 0.242023809517260 0.390165956326541 0.242023809517260 0.0575309935179400 0.00536221880153009; 0.00292854954264206 0.0384369616149092 0.196113758787202 0.382789297365280 0.287445156270025 0.0828345049143808 0.00945177150556065; 0.00153919355833585 0.0247005562212875 0.152923483594817 0.361483063911416 0.328567584707170 0.114741751017859 0.0160443669891157; 0.000778286512441806 0.0152660804766738 0.114741751017859 0.328567584707170 0.361483063911416 0.152923483594817 0.0262397497796233];

nos=size(CPROB)[1]; # total number of states
[sum(CPROB[i,:]) for i=1:nos];

mM=CPROB'-I;
mM[end,:]=[1 for i=1:nos];
rhs=[floor(i/nos) for i=1:nos];
BSP=mM\rhs; #the list of steady-state probabilities

hours = exp.([-1.20000000000000, -0.800000000000000, -0.400000000000000, 0, 0.400000000000000, 0.800000000000000, 1.20000000000000]);
wage=0.2;
inc=wage*hours;
ι=inc'*BSP; # payoff from the bond (= the aggregate output in the economy)

β=0.96; # impatience discount factor 
R=3.0; # risk aversion

