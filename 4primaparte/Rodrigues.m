function [v_rotated]= Rodrigues(u,v,delta)
v_rotated= [];
for i=1:length(u)
v_rot= v.*cos(delta)+cross(u(:,i),v).*sin(delta)+u(:,i).*dot(u(:,i),v).*(1-cos(delta));
v_rotated= [v_rotated,v_rot];
end
