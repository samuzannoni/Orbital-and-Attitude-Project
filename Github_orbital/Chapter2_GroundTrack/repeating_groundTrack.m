function [a]=repeating_groundTrack(k, m, mu_E, w_E)
    T=(m/k)*(2*pi/w_E);
    a=(mu_E*T^2/4*pi^2)^(1/3);
end