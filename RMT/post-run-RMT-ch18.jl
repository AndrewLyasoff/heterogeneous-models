######################################################################################################
#
# Julia code
#
# Illustrates the method described in Sec. 18.7 in "Recursive Macroeconomic Theory" (RMT)
#                                     by Lars Ljungqvist and Thomas Sargent
#
# This code supplements the paper "Another look at the distribution of income and wealth in the Macroeconomy" (DIW)
#                                        by Andrew Lyasoff (www.andrewlyasoff.tech)
#
# Copyright © 2019-2022 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
####################################################################################################

#=
If output-RMT-ch18-*.jls are available, there is no need for run-RMT-ch18.jl 
=#

# # this block is superfluous if run-RMT-ch18.jl
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

# this block is superfluous if run-RMT-ch18.jl
begin
    saved_vals=deserialize("output-RMT-ch18-200.jls");
    A1=saved_vals[1];
    B1=saved_vals[2];
    C1=saved_vals[3];
    D1p=saved_vals[4];
    E1p=saved_vals[5];
    D1m=saved_vals[6];
    E1m=saved_vals[7];
end;

# the left plot in Figure 1 in DIW
begin
    lvl=collect(minimum(A1):0.001:maximum(A1));
    lvl0=[0.0 for x in lvl];
    @gp A1[1] B1[1] "w p lt rgb 'black' pt 7 ps 5.0 t ''";
    @gp :- A1[2:end] B1[2:end] "w p lt rgb 'black' pt 7 ps 2.5 t ''";
    @gp :- lvl lvl0 "w l t '' lw 1.25 lt rgb 'black'";
    @gp :- "set auto fix";
    @gp :- "set offsets graph .085, graph .085, graph .085, graph .1";
    @gp :- "set grid lw 0.25 lc '#010000000'";
end

# the right plot in Figure 1 in DIW
begin
    lvl=collect(minimum(A1[10:end]):0.0000001:maximum(A1[10:end]));
    lvl0=[0.0 for x in lvl];
    @gp A1[10:end] B1[10:end] "w p lt rgb 'black' pt 7 ps 2.5 t ''";
    @gp :- lvl lvl0 "w l t '' lw 1.25 lt rgb 'black'";
    @gp :- "set auto fix";
    @gp :- "set offsets graph .085, graph .085, graph .085, graph .1";
    @gp :- "set grid lw 0.25 lc '#010000000'";
end

# this block is not needed after run-RMT-ch18.jl
begin
    saved_vals=deserialize("output-RMT-ch18-2K.jls");
    A2=saved_vals[1];
    B2=saved_vals[2];
    C2=saved_vals[3];
    D2p=saved_vals[4];
    E2p=saved_vals[5];
    D2m=saved_vals[6];
    E2m=saved_vals[7];
end;

# the left plot in Figure 2 in DIW
begin
    lvl=collect(minimum(A2):0.0000001:maximum(A2));
    lvl0=[0.0 for x in lvl];
    @gp A2[1] B2[1] "w p lt rgb 'black' pt 7 ps 5.0 t ''";
    @gp :- A2[2:end] B2[2:end] "w p lt rgb 'black' pt 7 ps 2.5 t ''";
    @gp :- lvl lvl0 "w l t '' lw 1.25 lt rgb 'black'";
    @gp :- "set auto fix";
    @gp :- "set offsets graph .085, graph .085, graph .085, graph .1";
    @gp :- "set grid lw 0.25 lc '#010000000'";
end

# the right plot in Figure 2 in DIW
begin
    lvl=(minimum(A2[10:end]):0.0000001:maximum(A2[10:end]));
    lvl0=[0.0 for x in lvl];
    @gp A2[10:end] B2[10:end] "w p lt rgb 'black' pt 7 ps 2.5 t ''";
    @gp :- lvl lvl0 "w l t '' lw 1.25 lt rgb 'black'";
    @gp :- "set auto fix";
    @gp :- "set offsets graph .085, graph .085, graph .085, graph .1";
    @gp :- "set grid lw 0.25 lc '#010000000'";
end


# demonstrate the discontinuity in the distribution as a function of r

Lg=2000

begin
    intr=A2[end]
    b=3.0;
    minkap=-b;
    if intr<=0.0
        minkap = -b;
    else
        minkap=-min(b, wage*hours[1]/intr);
    end;
    maxkap=16.0;
    grda=[minkap+i*(maxkap-minkap)/(Lg-1) for i=0:(Lg-1)];
    grdaa=[grda[i] for k=1:nos, i=1:length(grda)];
end;


begin
    bond_price=ι/(1+intr)
    cons=[(1+intr)*grdaa[k,i]+wage*hours[k]-grda[Int64(E2m[k,i])] for k=1:nos, i=1:Lg];
    aassets=[grdaa[k,Int64(E2m[k,i])] for k=1:nos, i=1:Lg];
    shares=[grdaa[k,Int64(E2m[k,i])]/bond_price for k=1:nos, i=1:Lg];
    SL=[sum(D2m[state,:]) for state=1:nos];
    LF=[[D2m[state,1]/SL[state]] for state=1:nos];
    for i=2:Lg
        for state=1:nos
            LF[state]=vcat(LF[state],LF[state][end]+D2m[state,i]/SL[state])
        end
    end
end

LgX=1300


# first part of Figure 3 in DIW
begin
    @gp aassets[1,2:LgX] LF[1][2:LgX] "w l t '' dashtype '-  -  ' lw 1.25 lt rgb 'black'";
    @gp :- "set auto fix";
    @gp :- "set offsets graph .025, graph .025, graph .025, graph .025";
    @gp :- "set grid lw 0.25 lc '#010000000'";
    @gp :- aassets[2,2:LgX] LF[2][2:LgX] "w l t '0.027423171997070310' lw 1.25 dashtype '-  -  ' lt rgb 'black'";
    @gp :- aassets[3,2:LgX] LF[3][2:LgX] "w l t '' dashtype '-  -  ' lw 1.25 lt rgb 'black'";
    @gp :- aassets[4,2:LgX] LF[4][2:LgX] "w l t '' dashtype '-  -  ' lw 1.25 lt rgb 'black'";
    @gp :- aassets[5,2:LgX] LF[5][2:LgX] "w l t '' dashtype '-  -  ' lw 1.25 lt rgb 'black'";
    @gp :- aassets[6,2:LgX] LF[6][2:LgX] "w l t '' dashtype '-  -  ' lw 1.25 lt rgb 'black'";
    @gp :- aassets[7,2:LgX] LF[7][2:LgX] "w l t '' dashtype '-  -  ' lw 1.25 lt rgb 'black'";
    @gp :- "set key bottom right"
end


begin
    intr=A2[end-1]
    b=3.0;
    minkap=-b;
    if intr<=0.0
        minkap = -b;
    else
        minkap=-min(b, wage*hours[1]/intr);
    end;
    maxkap=16.0;
    grda=[minkap+i*(maxkap-minkap)/(Lg-1) for i=0:(Lg-1)];
    grdaa=[grda[i] for k=1:nos, i=1:length(grda)];
end;


begin
    bond_price=ι/(1+intr)
    cons=[(1+intr)*grdaa[k,i]+wage*hours[k]-grda[Int64(E2p[k,i])] for k=1:nos, i=1:Lg];
    aassets=[grdaa[k,Int64(E2p[k,i])] for k=1:nos, i=1:Lg];
    shares=[grdaa[k,Int64(E2p[k,i])]/bond_price for k=1:nos, i=1:Lg];
    SL=[sum(D2p[state,:]) for state=1:nos];
    LF=[[D2p[state,1]/SL[state]] for state=1:nos];
    for i=2:Lg
        for state=1:nos
            LF[state]=vcat(LF[state],LF[state][end]+D2p[state,i]/SL[state])
        end
    end
end

#second part of Figure 3 in DIW
begin
    @gp :- aassets[1,2:5:LgX] LF[1][2:5:LgX] "w l t '0.027423248291015622' lw 1.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[2,2:5:LgX] LF[2][2:5:LgX] "w l t '' lw 1.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[3,2:5:LgX] LF[3][2:5:LgX] "w l t '' lw 1.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[4,2:5:LgX] LF[4][2:5:LgX] "w l t '' lw 1.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[5,2:5:LgX] LF[5][2:5:LgX] "w l t '' lw 1.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[6,2:5:LgX] LF[6][2:5:LgX] "w l t '' lw 1.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[7,2:5:LgX] LF[7][2:5:LgX] "w l t '' lw 1.25 dashtype '.  .  ' lt rgb 'black'";
end

