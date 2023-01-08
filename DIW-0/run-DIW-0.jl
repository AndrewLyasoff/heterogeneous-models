######################################################################################################
#
# Julia code
#
# Executes the method described in Sec. 3 of the paper
#       "Another look at the distribution of income and wealth in the Macroeconomy" (DIW) by Andrew Lyasoff
#
# The code provides a verifiable solution to the example from
#                Sec. 18.7 in "Recursive Macroeconomic Theory" (RMT) by Lars Ljungqvist and Thomas Sargent 
#
# This code supplements the paper "Another look at the distribution of income and wealth in the Macroeconomy" (DIW)
#                                        by Andrew Lyasoff (www.andrewlyasoff.com)
#
# Copyright © 2019-2023 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################


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

### INITIALIZATION

begin
    ρ=-log(β)
    # initial guess for the spot price
    SPOT=ι; # = face value = total income in the economy
    # initial guess for portfolio functions
    θ=[x->40.0*x-8.0 for k=1:nos];
    H=[x->x+(40.0*x-8.0)*SPOT for k=1:nos];
    # accordingly, initial choice for the inverse mapping in (3.3) 
    invH=[y->(y+8.0*SPOT)/(1+40.0*SPOT) for k=1:nos];
    # upper bound on the grid over consumption [Step 1 in (3.5)]
    c_high=FindUpperC(nos,SPOT,invH,CPROB,inc,ι,R);
    gsz=150 #The only place to set the number of grid points in the consumption range,
    g_step=(c_high/gsz);
end;



### THE MAIN ROUTINE IS NEXT

@time no_iter, conv_check, conv_check_c, clearing, upper_bb, F, Fval, θ, invH, SPOT, cnext, invcnext, g_step, gsz, sol, accu0, accu, ANSATZ, CLR=find_equil(nos,0,235,0.00001,g_step,gsz,θ,invH,SPOT,CPROB,inc,BSP,ι,R,549,299,199);
# 7941.893796 seconds (2h + 12 min)


#check the convergence
no_iter, conv_check, conv_check_c, clearing, accu0, accu
#=
returns: (235, 9.084474162968093e-5, 2.1632642125402057e-6, -1.7387766322006504e-6, 5.545564008002657e-14, 3.1645982615513546e-5)
=#

#### save the results produced by the main routine

#
# create the list of data to save
# N.B. spline objects cannot be saved and must be restored later from the respective abscissas and values
#

begin
    saved_vals=[CPROB, BSP, hours, wage, inc, ι, ρ, R, no_iter, conv_check, conv_check_c, clearing, g_step, gsz, upper_bb, Fval, sol, accu0, accu, ANSATZ, CLR, SPOT];
    serialize("output-DIW-0.jls", saved_vals); #dump to a file
end;
