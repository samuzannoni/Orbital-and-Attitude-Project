function [alpha, delta] = car2RADec(R)

Rnorm = vecnorm(R')'; % vecnorm() calculates column-wise norm, but we need it as row-wise

% Declination
delta = asin(R(:,3)./Rnorm);

% Right Ascension
alpha = acos(R(:,1)./(Rnorm.*cos(delta)));
if  (R(:,2)./Rnorm) <= 0
    alpha = 2*pi - alpha;
end

end
