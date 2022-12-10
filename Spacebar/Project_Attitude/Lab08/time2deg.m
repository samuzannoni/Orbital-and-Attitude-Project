function [deg] = time2deg(hh,mm,ss)

    deg = (hh+mm/60+ss/3600)*360/24;

end
