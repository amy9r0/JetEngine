function [gam]=Gamma(T)
%T-     temperature (K)
%gam-   ratio of specific heats
if T<=600
    gam=    1.4;
else
    gam=    1.936-0.0838*log(T);
end