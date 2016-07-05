#Crude translation of
#http://www.mathworks.com/matlabcentral/fileexchange/22545-most-probable-number--mpn-
#itself a translation of
# http://www.i2workout.com/mcuriale/mpn/index.html


 MPNindex <-function( Dilution,Replicates,Positive)
 {  [MPN.R MPN.H MPN.L]
#  MPNindex <-function( Dilution,Replicates,Positive)
 {  [MPN.R MPN.H MPN.L]
# Temporary Bug Fix  - when division by 0 multiplying Dilution by 10
# 11.8.08
# -------------------------------------------------------------------------
# Purpose: getting the number of bactria in the original sample
# description: mpn algorithm
# Arguments: Dilution - an array contains the amount of original sample
#       Replicates - an array of the number of wells to concider (usually 8)
#       Positive - an array of the number of wells containing bact.
# Returns: MPN.R - how many bact were in the original sample
#       MPN.H, MPN.L - the error of the calculation - upper and lower limit
# -------------------------------------------------------------------------
#translated from Excel version MPN Calculator by by Mike Curiale
#(http://www.i2workout.com/mcuriale/mpn/index.html mcuriale@gmail.com)
#translated by Ofer Fridman (oferfrid@hotmail.com)
DilutionFactor<-1
numdil <- length(Dilution)

precision <- 0.0000000001
a <- precision


#loop to fix BUG
ConstDilution <- Dilution

flag <- 1
while flag==1
    Lower <- 0
Upper <- 1000000000
MPN <- Upper
#  'start values to    } # routine
vdp <- 2
vdn <- 1
    Dilution <- ConstDilution.*DilutionFactor
    flag<-0
    while vdn < vdp
        Upper <- Upper / 10
        MPN <- Upper
        vdp <- 0
        vdn <- 0
        for (  j in 1 : numdil
 ) {
            if exp(-Dilution( j) * MPN)==1
                DilutionFactor<-DilutionFactor*10
                flag<-1
                break
               } #
            vdn <- vdn + Dilution( j) * Replicates( j)
            vdp <- vdp + (Dilution( j) * Positive( j) / (1 - exp(-Dilution( j) * MPN)))
           } #
        if flag==1
            break
           } #
       } #
   } #



while vdn < vdp
    Upper <- Upper / 10
    MPN <- Upper
    vdp <- 0
    vdn <- 0
    for (  j in 1 : numdil
 ) {
        vdn <- vdn + Dilution( j) * Replicates( j)
        vdp <- vdp + (Dilution( j) * Positive( j) / (1 - exp(-Dilution( j) * MPN)))
       } #
   } #



Upper <- Upper * 10
while ((vdn - vdp) > a) || ((vdn - vdp) < -a)
    Midpoint <- (Lower + Upper) / 2
    MPN <- Midpoint
    #'reset vdn and vdp
    vdn <- 0
    vdp <- 0
    for (  j  in  1 : numdil
 ) {
        vdn <- vdn + Dilution( j) * Replicates( j)
        vdp <- vdp + (Dilution( j) * Positive( j) / (1 - exp(-Dilution( j) * MPN)))
        if MPN < 0.00001
            fprintf('Overflow calculation. Could not compute result. Reduce dilutions or number of tubes.')
           } #
       } #

    if vdp>vdn
        Lower <- Midpoint
       } #

    if vdp < vdn
        Upper <- Midpoint
       } #
   } #

#     'excel 97 does not have round function
MPN.R <- sigfig(MPN, 2)
g2n <- 0
#     Select Case MPNresult
#     Case 2  'return  lower limit of 95# ci
for (  j  in  1 : numdil
 ) {
    m <- (Dilution( j) ^ 2 * Replicates( j))
    n <- (exp(Dilution( j) * MPN) - 1)
    g2n <- g2n + m / n
   } #

g2n <- (1 / (MPN ^ 2 * g2n)) ^ 0.5 / log(10)
g2n <- 10 ^ ((log(MPN) / log(10) - 1.96 * g2n))
MPN.L <- sigfig(g2n, 2)

#     Case 3   'return upper limit of 95# ci
g2n <- 0
for (  j in 1 : numdil
 ) {
    m <- (Dilution( j) ^ 2 * Replicates( j))
    n <- (exp(Dilution( j) * MPN) - 1)
    g2n <- g2n + m / n
   } #

g2n <- (1 / (MPN ^ 2 * g2n)) ^ 0.5 / log(10)
g2n <- 10 ^ ((log(MPN) / log(10) + 1.96 * g2n))
MPN.H <- sigfig(g2n, 2)

#BUG FIX
MPN.R<-MPN.R*DilutionFactor
MPN.L<-MPN.L*DilutionFactor
MPN.H <- MPN.H*DilutionFactor


   } #

  sigfig <-function( GenNum, SigNum )
 {  sigfig
b<-  floor(log10(GenNum)) -(SigNum -1)
c<-GenNum/(10^b)
fc<-floor(c)
d<- c - fc
if d ><- 0.5,
    fc<-fc+1
   } #
#     'return value
sigfig <- fc*10^b
   } #
