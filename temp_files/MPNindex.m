function [MPN_R MPN_H MPN_L]= MPNindex(Dilution,Replicates,Positive)
% function [MPN_R MPN_H MPN_L]= MPNindex(Dilution,Replicates,Positive)
% Temporary Bug Fix  - when division by 0 multiplying Dilution by 10
% 11.8.08
% -------------------------------------------------------------------------
% Purpose: getting the number of bactria in the original sample
% description: mpn algorithm
% Arguments: Dilution - an array contains the amount of original sample
%       Replicates - an array of the number of wells to concider (usually 8)
%       Positive - an array of the number of wells containing bact.
% Returns: MPN_R - how many bact were in the original sample
%       MPN_H, MPN_L - the error of the calculation - upper and lower limit
% -------------------------------------------------------------------------
%translated from Excel version MPN Calculator by by Mike Curiale
%(http://www.i2workout.com/mcuriale/mpn/index.html mcuriale@gmail.com)
%translated by Ofer Fridman (oferfrid@hotmail.com)
DilutionFactor=1;
numdil = length(Dilution);

precision = 0.0000000001;
a = precision;


%loop to fix BUG
ConstDilution = Dilution;

flag = 1;
while flag==1
    Lower = 0;
Upper = 1000000000;
MPN = Upper;
%  'start values to end routine
vdp = 2;
vdn = 1;
    Dilution = ConstDilution.*DilutionFactor;
    flag=0;
    while vdn < vdp
        Upper = Upper / 10;
        MPN = Upper;
        vdp = 0;
        vdn = 0;
        for j=1:numdil
            if exp(-Dilution( j) * MPN)==1
                DilutionFactor=DilutionFactor*10;
                flag=1;
                break
            end
            vdn = vdn + Dilution( j) * Replicates( j);
            vdp = vdp + (Dilution( j) * Positive( j) / (1 - exp(-Dilution( j) * MPN)));
        end
        if flag==1
            break
        end
    end
end



while vdn < vdp
    Upper = Upper / 10;
    MPN = Upper;
    vdp = 0;
    vdn = 0;
    for j=1:numdil
        vdn = vdn + Dilution( j) * Replicates( j);
        vdp = vdp + (Dilution( j) * Positive( j) / (1 - exp(-Dilution( j) * MPN)));
    end
end



Upper = Upper * 10;
while ((vdn - vdp) > a) || ((vdn - vdp) < -a)
    Midpoint = (Lower + Upper) / 2;
    MPN = Midpoint;
    %'reset vdn and vdp
    vdn = 0;
    vdp = 0;
    for j = 1:numdil
        vdn = vdn + Dilution( j) * Replicates( j);
        vdp = vdp + (Dilution( j) * Positive( j) / (1 - exp(-Dilution( j) * MPN)));
        if MPN < 0.00001
            fprintf('Overflow calculation. Could not compute result. Reduce dilutions or number of tubes.');
        end
    end

    if vdp>vdn
        Lower = Midpoint;
    end

    if vdp < vdn
        Upper = Midpoint;
    end
end

%     'excel 97 does not have round function
MPN_R = sigfig(MPN, 2);
g2n = 0;
%     Select Case MPNresult
%     Case 2  'return  lower limit of 95% ci
for j = 1:numdil
    m = (Dilution( j) ^ 2 * Replicates( j));
    n = (exp(Dilution( j) * MPN) - 1);
    g2n = g2n + m / n;
end

g2n = (1 / (MPN ^ 2 * g2n)) ^ 0.5 / log(10);
g2n = 10 ^ ((log(MPN) / log(10) - 1.96 * g2n));
MPN_L = sigfig(g2n, 2);

%     Case 3   'return upper limit of 95% ci
g2n = 0;
for j=1:numdil
    m = (Dilution( j) ^ 2 * Replicates( j));
    n = (exp(Dilution( j) * MPN) - 1);
    g2n = g2n + m / n;
end

g2n = (1 / (MPN ^ 2 * g2n)) ^ 0.5 / log(10);
g2n = 10 ^ ((log(MPN) / log(10) + 1.96 * g2n));
MPN_H = sigfig(g2n, 2);

%BUG FIX
MPN_R=MPN_R*DilutionFactor;
MPN_L=MPN_L*DilutionFactor;
MPN_H = MPN_H*DilutionFactor;


end

function sigfig =  sigfig(GenNum, SigNum )
b=  floor(log10(GenNum)) -(SigNum -1);
c=GenNum/(10^b);
fc=floor(c);
d= c - fc;
if d >= 0.5,
    fc=fc+1;
end
%     'return value
sigfig = fc*10^b;
end
