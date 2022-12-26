######################################################################################################
#
# Julia code
#
# Executes the method described in Sec. 18.7 in "Recursive Macroeconomic Theory" (RMT)
#                                     by Lars Ljungqvist and Thomas Sargent
#
# The code included here is an adaptation of the MATLAB code that accompanies RMT.
#
# It is discussed in the paper "Another look at the distribution of income and wealth in the Macroeconomy" (DIW)
#                                        by Andrew Lyasoff  (www.andrewlyasoff.tech)
# 
# Copyright © 2019-2022 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
####################################################################################################


begin
    using Serialization, FileIO
    using LinearAlgebra
    using Interpolations
    using Roots
    using Gnuplot
    using Optim
    include("functions-RMT-ch18.jl");
    include("ini-setup-RMT-ch18.jl");
end;

#=
The following are assigned from   ini-setup-RMT-ch18.jl

CPROB, nos, BSP, hours, wage, ι, ρ, β, R
=#


# 20 trials with 200 grid points
# initial trail rates are [0.03701851068729933,0.03,0.02]
@time A1, B1, C1, D1p, E1p, D1m, E1m = iterLS(20,5,5,nos,hours,wage,β,200,CPROB,R,ι,[0.03701851068729933,0.03,0.02],3.0,16.0,BSP);
# 73.198813 seconds

begin
    saved_vals=[A1,B1,C1,D1p,E1p,D1m,E1m];
    #dump to a file
    serialize("output-RMT-ch18-200.jls", saved_vals);
end;

# 20 trials with 2000 grid points
# initial trail rates are [0.03701851068729933,0.03,0.02]
@time A2, B2, C2, D2p, E2p, D2m, E2m = iterLS(20,5,5,nos,hours,wage,β,2000,CPROB,R,ι,[0.03701851068729933,0.03,0.02],3.0,16.0,BSP);
#4918.069815 seconds (1h 22min)

begin
    saved_vals=[A2,B2,C2,D2p,E2p,D2m,E2m];
    #dump to a file
    serialize("output-RMT-ch18-2K.jls", saved_vals);
end;
