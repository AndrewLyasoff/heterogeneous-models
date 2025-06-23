######################################################################################################
#
# Julia code
#
# Sets the parameter values for the benchmark economy described in:
#       Krusell, Per, and Anthony Smith. (1998). Income and wealth heterogeneity in the macroeconomy.
#                         Journal of Political Economy 106 867-896.
#                
# This code supplements the paper "Self-Consistent Transport in Heterogeneous Agent Models" [STHAM]
#                                        by Andrew Lyasoff
#
# Copyright © 2019-2025 Andrew Lyasoff <mathema@lyasoff.net>
# SPDX-License-Identifier: Apache-2.0
#
###################################################################################################


function FGQT(f,a,b)
    return (arg->(1/2)*(b-a)*f( (1/2)*(b-a)*arg + (1/2)*(a+b) ))
end


## get the steady-state probabilities for a given transition matrix
function STSP(trans_prob::Array{Float64,2})
    local nos::Int64, mM::Array{Float64,2}, rhs::Array{Float64,1};
    nos=size(trans_prob)[1];
    if maximum(abs.([sum(trans_prob[i,:]) for i=1:nos]-ones(nos)))>10.0e-16
        println("not a TPM")
    end
    mM=[(trans_prob-I)';ones(nos)'];
    rhs=push!(zeros(nos),1.0);
    return (mM'*mM)\(mM'*rhs)
end



#transition probabilities borrowed from the FORTRAN code from KS

begin
    β=0.96;δ=0.025;α=0.36; #modified from the original β=0.99;δ=0.025;α=0.36; 
    durug=1.5;durub=2.5;unempg=0.04;unempb=0.1;
    durbd=8.0;durbd=8.0;durub=2.5;durgd=8.0;
    #
    pgg00 = (durug-1.0)/durug
    pbb00 = (durub-1.0)/durub
    pbg00 = 1.25*pbb00
    pgb00 = 0.75*pgg00
    pgg01 = (unempg - unempg*pgg00)/(1.0-unempg)
    pbb01 = (unempb - unempb*pbb00)/(1.0-unempb)
    pbg01 = (unempb - unempg*pbg00)/(1.0-unempg)
    pgb01 = (unempg - unempb*pgb00)/(1.0-unempb)
    pgg = (durgd-1.0)/durgd
    pgb = 1.0 - (durbd-1.0)/durbd
    #
    pgg10 = 1.0 - (durug-1.0)/durug
    pbb10 = 1.0 - (durub-1.0)/durub
    pbg10 = 1.0 - 1.25*pbb00
    pgb10 = 1.0 - 0.75*pgg00
    pgg11 = 1.0 - (unempg - unempg*pgg00)/(1.0-unempg)
    pbb11 = 1.0 - (unempb - unempb*pbb00)/(1.0-unempb)
    pbg11 = 1.0 - (unempb - unempg*pbg00)/(1.0-unempg)
    pgb11 = 1.0 - (unempg - unempb*pgb00)/(1.0-unempb)
    pbg = 1.0 - (durgd-1.0)/durgd
    pbb = (durbd-1.0)/durbd
    #
end;

pr=zeros(4,4)

#from the FORTRAN code the transpose of the transition
begin
    pr[1,1] = pgg*pgg11;
    pr[2,1] = pbg*pbg11;
    pr[3,1] = pgg*pgg01;
    pr[4,1] = pbg*pbg01;
    pr[1,2] = pgb*pgb11;
    pr[2,2] = pbb*pbb11;
    pr[3,2] = pgb*pgb01;
    pr[4,2] = pbb*pbb01;
    pr[1,3] = pgg*pgg10;
    pr[2,3] = pbg*pbg10;
    pr[3,3] = pgg*pgg00;
    pr[4,3] = pbg*pbg00;
    pr[1,4] = pgb*pgb10;
    pr[2,4] = pbb*pbb10;
    pr[3,4] = pgb*pgb00;
    pr[4,4] = pbb*pbb00;
end;



#the order of states is 1g 1b 0g 0b
#the order of states is 1g 1b 0g 0b
# 1g:1g 1g:1b 1g:0g 1g:0b
# 1b:1g 1b:1b 1b:0g 1b:0b
# 0g:1g 0g:1b 0g:0g 0g:0b
# 0b:1g 0b:1b 0b:0g 0b:0b

begin
    atpm=[pgg pbg;pgb pbb];
    tpmgg=[pgg11 pgg01;pgg10 pgg00];
    tpmbg=[pbg11 pbg01;pbg10 pbg00];
    tpmgb=[pgb11 pgb01;pgb10 pgb00];
    tpmbb=[pbb11 pbb01;pbb10 pbb00];
end;

begin
    itpm=zeros(2,2,2,2);
    itpm[:,:,2,2]=tpmbb;
    itpm[:,:,2,1]=tpmgb;
    itpm[:,:,1,2]=tpmbg;
    itpm[:,:,1,1]=tpmgg;
end;

#[sum(itpm[i,:,k,l]) for i=1:2,k=1:2,l=1:2]

astatesi=vec([(k,l) for l=1:2,k=1:2]);

tpm=[atpm[s1[2],s2[2]]*itpm[s1[1],s2[1],s1[2],s2[2]] for s1 in astatesi,s2 in astatesi];

sspa=STSP(atpm);

sspm=STSP(tpm);

pjl=[sspm[i]/sspa[astatesi[i][2]] for i=1:4];

function pos(j::Int64,l::Int64,asttsi::Vector{Tuple{Int64, Int64}})
    return findfirst(isequal((j,l)),asttsi)
end

PJL=[sspm[pos(j,l,astatesi)]/sspa[l]  for j=1:2,l=1:2]

hfix=0.3271;

#SS=[hfix*(1-0.04),hfix*(1-0.1)]; #old one

SS=[hfix,0.0]

XX=[1.01,0.99];

#SSXX=[(SS[s[1]],XX[s[2]]) for s in astatesi];

#NN=[SS[1]*pjl[1],SS[2]*pjl[2]]; #old one

NN=[hfix*pjl[1],hfix*pjl[2]]


#KK=[30.0,25.0];

