function [km] = ly2km(ly)
    
    conversion_constant = 9.46e12; %[km]
    km = ly*conversion_constant;

end
