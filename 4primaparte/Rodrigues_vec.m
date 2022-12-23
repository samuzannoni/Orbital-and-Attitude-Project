function [v_rot]= Rodrigues_vec(u,v,delta)
v_rot= v.*cos(delta)+cross(u,v).*sin(delta)+u.*dot(u,v).*(1-cos(delta));
end
