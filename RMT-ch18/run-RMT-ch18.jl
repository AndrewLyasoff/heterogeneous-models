######################################################################################################
#
# Julia code
#
# Executes the method described in Sec. 18.7 in "Recursive Macroeconomic Theory" (RMT)
#                                     by Lars Ljungqvist and Thomas Sargent
#
# The code included here is an adaptation of the MATLAB code that accompanies RMT.
#
# It is discussed in the paper "Dynamic transportation of economic agents" [DTEA]
#                  by Andrew Lyasoff  (www.andrewlyasoff.tech)
# 
# Copyright © 2019-2023 Andrew Lyasoff <alyasoff@bu.edu>
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
#   69.069408 seconds (60.21 M allocations: 195.322 GiB, 6.44% gc time, 0.06% compilation time)

begin
    saved_vals=[A1,B1,C1,D1p,E1p,D1m,E1m];
    #dump to a file
    serialize("output-RMT-ch18-200.jls", saved_vals);
end;

# 20 trials with 2000 grid points
# initial trail rates are [0.03701851068729933,0.03,0.02]
@time A2, B2, C2, D2p, E2p, D2m, E2m = iterLS(20,5,5,nos,hours,wage,β,2000,CPROB,R,ι,[0.03701851068729933,0.03,0.02],3.0,16.0,BSP);
# 5117.190374 seconds (588.66 M allocations: 11.180 TiB, 2.08% gc time) 1h + 25min

begin
    saved_vals=[A2,B2,C2,D2p,E2p,D2m,E2m];
    #dump to a file
    serialize("output-RMT-ch18-2K.jls", saved_vals);
end;

# 3 trials with 3000 grid points
# initial trail rates are [0.024407690878163295,0.024407690847486237]
@time A3, B3, C3, D3p, E3p, D3m, E3m = iterLS(3,5,5,nos,hours,wage,β,3000,CPROB,R,ι,[0.024407690878163295,0.024407690847486237],3.0,16.0,BSP);
# 1590.487496 seconds (222.05 M allocations: 3.365 TiB, 1.70% gc time)

begin
    saved_vals=[A3,B3,C3,D3p,E3p,D3m,E3m];
    #dump to a file
    serialize("output-RMT-ch18-3K.jls", saved_vals);
end;


# 3 trials with 4000 grid points
# initial trail rates are [0.02151987413956899,0.021519782586834617]
@time A4, B4, C4, D4p, E4p, D4m, E4m = iterLS(3,5,5,nos,hours,wage,β,4000,CPROB,R,ι,[0.02151987413956899,0.021519782586834617],3.0,16.0,BSP);
# 3257.718288 seconds (313.38 M allocations: 6.801 TiB, 1.35% gc time)

begin
    saved_vals=[A4,B4,C4,D4p,E4p,D4m,E4m];
    #dump to a file
    serialize("output-RMT-ch18-4K.jls", saved_vals);
end;
