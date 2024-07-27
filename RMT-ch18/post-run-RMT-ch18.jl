######################################################################################################
#
# Julia code
#
# Illustrates the method described in Sec. 18.7 in "Recursive Macroeconomic Theory" (RMT)
#                                     by Lars Ljungqvist and Thomas Sargent
#
# This code supplements the paper "Self-Aware Transport of Economic Agents" [SATEA]
#                                        by Andrew Lyasoff (www.andrewlyasoff.tech)
#
# Copyright © 2019-2024 Andrew Lyasoff <alyasoff@bu.edu>
# SPDX-License-Identifier: Apache-2.0
#
####################################################################################################

#=
If output-RMT-ch18-*.jls are available, there is no need for run-RMT-ch18.jl 
=#

# this block is not needed after run-RMT-ch18.jl
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

## experiment with 200 grid points

# this block is not needed after run-RMT-ch18.jl
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

# (rate,expected demand)
[(A1[i],B1[i]) for i=1:length(A1)]
#=
20-element Vector{Tuple{Float64, Float64}}:
 (0.03701851068729933, 7.3696566090458155)
 (0.03, 0.510066845946097)
 (0.02, -1.9729997109227826)
 (0.025, -1.0654201141611443)
 (0.0275, -0.505665361435287)
 (0.028749999999999998, -0.1113888746946741)
 (0.029375, 0.16311610386135325)
 (0.029062499999999998, 0.02379514486718422)
 (0.028906249999999998, -0.04304169618186463)
 (0.028984375, 0.01870022538380367)
 (0.0289453125, 0.015794907576101393)
 (0.028925781249999998, -0.04179338316455158)
 (0.028935546875, 0.01539762844520665)
 (0.0289306640625, -0.04146948169368391)
 (0.028933105468749998, -0.04131361382945803)
 (0.028934326171874997, -0.04123568976115094)
 (0.028934936523437496, -0.04119673019258489)
 (0.028935241699218746, -0.04117725102465344)
 (0.028935394287109374, 0.01538792345543593)
 (0.02893531799316406, 0.015383070922166776)
=#

A1[end]-A1[end-2],B1[end]-B1[end-2]
# (7.629394531416533e-8, 0.056560321946820216)


# the left plot in Figure 1 in [SATEA]
begin
    lvl=collect(minimum(A1):0.001:maximum(A1));
    lvl0=[0.0 for x in lvl];
    @gp [A1[1]] [B1[1]] "w p lt rgb 'black' pt 7 ps 5.0 t ''";
    @gp :- A1[2:end] B1[2:end] "w p lt rgb 'black' pt 7 ps 2.5 t ''";
    @gp :- lvl lvl0 "w l t '' lw 1.25 lt rgb 'black'";
end


# the right plot in Figure 1 in [SATEA]
begin
    lvl=collect(minimum(A1[10:end]):0.0000001:maximum(A1[10:end]));
    lvl0=[0.0 for x in lvl];
    @gp A1[10:end] B1[10:end] "w p lt rgb 'black' pt 7 ps 2.5 t ''";
    @gp :- lvl lvl0 "w l t '' lw 1.25 lt rgb 'black'";
end

## experiment with 2K grid points

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


# (rate,expected demand)
[(A2[i],B2[i]) for i=1:length(A1)]
#=
20-element Vector{Tuple{Float64, Float64}}:
 (0.03701851068729933, 6.901030700349697)
 (0.03, 4.752922406432852)
 (0.02, -2.3760213752509043)
 (0.025, -1.660077615449614)
 (0.0275, 1.9612350672273267)
 (0.026250000000000002, -1.4970516328079035)
 (0.026875000000000003, -1.4244100661711903)
 (0.027187500000000003, -1.384199047345656)
 (0.02734375, -1.3647950282856893)
 (0.027421875, -1.3540333498512596)
 (0.027460937499999998, 1.9123408467470753)
 (0.027441406249999998, 1.9052187085592882)
 (0.027431640624999996, 1.9014587979701378)
 (0.027426757812499997, 1.8995352213298011)
 (0.027424316406249996, 1.8985398756228542)
 (0.027423095703124997, -1.353910586114863)
 (0.027423706054687497, 1.898344736653897)
 (0.027423400878906247, 1.898108915981646)
 (0.027423248291015622, 1.8980739113636893)
 (0.02742317199707031, -1.3539040644306715)
=#

A2[end]-A2[end-2],B2[end]-B2[end-2]
# (-2.288818359355571e-7, -3.2520129804123172)


# the left plot in Figure 2 in [SATEA]
begin
    lvl=collect(minimum(A2):0.0000001:maximum(A2));
    lvl0=[0.0 for x in lvl];
    @gp [A2[1]] [B2[1]] "w p lt rgb 'black' pt 7 ps 5.0 t ''";
    @gp :- A2[2:end] B2[2:end] "w p lt rgb 'black' pt 7 ps 2.5 t ''";
    @gp :- lvl lvl0 "w l t '' lw 1.25 lt rgb 'black'";
end


# the right plot in Figure 2 in [SATEA]
begin
    lvl=(minimum(A2[10:end]):0.0000001:maximum(A2[10:end]));
    lvl0=[0.0 for x in lvl];
    @gp A2[10:end] B2[10:end] "w p lt rgb 'black' pt 7 ps 2.5 t ''";
    @gp :- lvl lvl0 "w l t '' lw 1.25 lt rgb 'black'";
end


# demonstrate the discontinuity in the distribution as a function of r with 2K grid points

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

intr
# 0.02742317199707031


# first part of Figure 3 Left in [SATEA]
begin
    @gp aassets[1,2:LgX] LF[1][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    #@gp :- "set auto fix";
    @gp :- aassets[2,2:LgX] LF[2][2:LgX] "w l t '0.027423171997070310' lw 2.25 dashtype '-  -  ' lt rgb 'black'";
    @gp :- aassets[3,2:LgX] LF[3][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[4,2:LgX] LF[4][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[5,2:LgX] LF[5][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[6,2:LgX] LF[6][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[7,2:LgX] LF[7][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- "set key bottom right font 'Linux Libertine O,40'"
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

intr
# 0.027423248291015622

#second part of Figure 3 Left in [SATEA]
begin
    @gp :- aassets[1,2:5:LgX] LF[1][2:5:LgX] "w l t '0.027423248291015622' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[2,2:5:LgX] LF[2][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[3,2:5:LgX] LF[3][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[4,2:5:LgX] LF[4][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[5,2:5:LgX] LF[5][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[6,2:5:LgX] LF[6][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[7,2:5:LgX] LF[7][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
end

## experiment with 3K grid points

# this block is not needed after run-RMT-ch18.jl
begin
    saved_vals=deserialize("output-RMT-ch18-3K.jls");
    A3=saved_vals[1];
    B3=saved_vals[2];
    C3=saved_vals[3];
    D3p=saved_vals[4];
    E3p=saved_vals[5];
    D3m=saved_vals[6];
    E3m=saved_vals[7];
end;

# (rate,expected demand)
[(A3[i],B3[i]) for i=1:length(A3)]
#=
3-element Vector{Tuple{Float64, Float64}}:
 (0.024407690878163295, 1.8459023007560789)
 (0.024407690847486237, -1.7162724004974983)
 (0.024407690862824766, 1.8459022995673882)
=#

A3[end]-A3[end-1],B3[end]-B3[end-1]
# (1.5338529057995487e-11, 3.5621747000648867)

# demonstrate the discontinuity in the distribution as a function of r 

Lg=3000

begin
    intr=A3[end]
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
    cons=[(1+intr)*grdaa[k,i]+wage*hours[k]-grda[Int64(E3m[k,i])] for k=1:nos, i=1:Lg];
    aassets=[grdaa[k,Int64(E3m[k,i])] for k=1:nos, i=1:Lg];
    shares=[grdaa[k,Int64(E3m[k,i])]/bond_price for k=1:nos, i=1:Lg];
    SL=[sum(D3m[state,:]) for state=1:nos];
    LF=[[D3m[state,1]/SL[state]] for state=1:nos];
    for i=2:Lg
        for state=1:nos
            LF[state]=vcat(LF[state],LF[state][end]+D3m[state,i]/SL[state])
        end
    end
end

LgX=2000

intr
# 0.024407690862824766

# first part of Figure 3 Right in [SATEA]
begin
    @gp aassets[1,2:LgX] LF[1][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    #@gp :- "set auto fix";
    @gp :- aassets[2,2:LgX] LF[2][2:LgX] "w l t '0.024407690862824766' lw 2.25 dashtype '-  -  ' lt rgb 'black'";
    @gp :- aassets[3,2:LgX] LF[3][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[4,2:LgX] LF[4][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[5,2:LgX] LF[5][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[6,2:LgX] LF[6][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[7,2:LgX] LF[7][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- "set key bottom right font 'Linux Libertine O,40'"
end

begin
    intr=A3[end-1]
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
    cons=[(1+intr)*grdaa[k,i]+wage*hours[k]-grda[Int64(E3p[k,i])] for k=1:nos, i=1:Lg];
    aassets=[grdaa[k,Int64(E3p[k,i])] for k=1:nos, i=1:Lg];
    shares=[grdaa[k,Int64(E3p[k,i])]/bond_price for k=1:nos, i=1:Lg];
    SL=[sum(D3p[state,:]) for state=1:nos];
    LF=[[D3p[state,1]/SL[state]] for state=1:nos];
    for i=2:Lg
        for state=1:nos
            LF[state]=vcat(LF[state],LF[state][end]+D3p[state,i]/SL[state])
        end
    end
end

intr
# 0.024407690847486237

#second part of Figure 3 Right in [SATEA]
begin
    @gp :- aassets[1,2:5:LgX] LF[1][2:5:LgX] "w l t '0.024407690847486237' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[2,2:5:LgX] LF[2][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[3,2:5:LgX] LF[3][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[4,2:5:LgX] LF[4][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[5,2:5:LgX] LF[5][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[6,2:5:LgX] LF[6][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[7,2:5:LgX] LF[7][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
end


## experiment with 4K grid points

# this block is not needed after run-RMT-ch18.jl
begin
    saved_vals=deserialize("output-RMT-ch18-4K.jls");
    A4=saved_vals[1];
    B4=saved_vals[2];
    C4=saved_vals[3];
    D4p=saved_vals[4];
    E4p=saved_vals[5];
    D4m=saved_vals[6];
    E4m=saved_vals[7];
end;

# (rate,expected demand)
[(A4[i],B4[i]) for i=1:length(A4)]
#=
julia> 3-element Vector{Tuple{Float64, Float64}}:
 (0.02151987413956899, 0.9481247517259364)
 (0.021519782586834617, -2.124200775454897)
 (0.021519828363201803, -2.124193122792572)
=#

A4[end]-A4[end-2],B4[end]-B4[end-2]
# (-4.577636718641753e-8, -3.0723178745185082)

# SOME ADDITIONAL RUNS

# 1 trial with 4000 grid points
@time A4x, B4x, C4x, D4px, E4px, D4mx, E4mx = iterLS(1,5,5,nos,hours,wage,β,4000,CPROB,R,ι,[0.021519851251385398],3.0,16.0,BSP);
#=
0.021519851251385398
     -2.12419013291833
bad initial choice
1141.987274 seconds (107.02 M allocations: 2.564 TiB, 1.10% gc time, 0.00% compilation time)
=#


# 1 trial with 4000 grid points
@time A4y, B4y, C4y, D4py, E4py, D4my, E4my = iterLS(1,5,5,nos,hours,wage,β,4000,CPROB,R,ι,[0.021519862695477194],3.0,16.0,BSP);
#=
0.021519862695477194
     0.9481235598541813
bad initial choice
869.767462 seconds (98.10 M allocations: 1.657 TiB, 1.45% gc time)
=#

0.021519862695477194-0.021519851251385398
# 1.1444091795737021e-8

0.9481235598541813-(-2.12419013291833)
# 3.0723136927725117

# demonstrate the discontinuity in the distribution as a function of r 

Lg=4000

begin
    intr=A4[end]
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
    cons=[(1+intr)*grdaa[k,i]+wage*hours[k]-grda[Int64(E4m[k,i])] for k=1:nos, i=1:Lg];
    aassets=[grdaa[k,Int64(E4m[k,i])] for k=1:nos, i=1:Lg];
    shares=[grdaa[k,Int64(E4m[k,i])]/bond_price for k=1:nos, i=1:Lg];
    SL=[sum(D4m[state,:]) for state=1:nos];
    LF=[[D4m[state,1]/SL[state]] for state=1:nos];
    for i=2:Lg
        for state=1:nos
            LF[state]=vcat(LF[state],LF[state][end]+D4m[state,i]/SL[state])
        end
    end
end

LgX=2667

intr
# 0.021519828363201803

# first part of plot not included in [SATEA]
begin
    @gp aassets[1,2:LgX] LF[1][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    #@gp :- "set auto fix";
    @gp :- aassets[2,2:LgX] LF[2][2:LgX] "w l t '0.021519828363201803' lw 2.25 dashtype '-  -  ' lt rgb 'black'";
    @gp :- aassets[3,2:LgX] LF[3][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[4,2:LgX] LF[4][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[5,2:LgX] LF[5][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[6,2:LgX] LF[6][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- aassets[7,2:LgX] LF[7][2:LgX] "w l t '' dashtype '-  -  ' lw 2.25 lt rgb 'black'";
    @gp :- "set key bottom right font 'Linux Libertine O,40'"
end


begin
    intr=A4[end-2]
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
    cons=[(1+intr)*grdaa[k,i]+wage*hours[k]-grda[Int64(E4p[k,i])] for k=1:nos, i=1:Lg];
    aassets=[grdaa[k,Int64(E4p[k,i])] for k=1:nos, i=1:Lg];
    shares=[grdaa[k,Int64(E4p[k,i])]/bond_price for k=1:nos, i=1:Lg];
    SL=[sum(D4p[state,:]) for state=1:nos];
    LF=[[D4p[state,1]/SL[state]] for state=1:nos];
    for i=2:Lg
        for state=1:nos
            LF[state]=vcat(LF[state],LF[state][end]+D4p[state,i]/SL[state])
        end
    end
end

intr
# 0.02151987413956899

#second part of the plot not included in [SATEA]
begin
    @gp :- aassets[1,2:5:LgX] LF[1][2:5:LgX] "w l t '0.021519874139568990' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[2,2:5:LgX] LF[2][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[3,2:5:LgX] LF[3][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[4,2:5:LgX] LF[4][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[5,2:5:LgX] LF[5][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[6,2:5:LgX] LF[6][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
    @gp :- aassets[7,2:5:LgX] LF[7][2:5:LgX] "w l t '' lw 2.25 dashtype '.  .  ' lt rgb 'black'";
end

