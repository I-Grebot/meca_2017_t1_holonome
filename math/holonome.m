R=30
r=6
la=87
teta=45
j1a=[-R;0;0];
j1u=[r*sind(90+teta) la;-r*cosd(90+teta) 0;0 1];
delta_j1u=j1u*inv(transpose(j1u)*j1u)*transpose(j1u)-eye(3);
inv(transpose(j1a)*delta_j1u*j1a)*transpose(j1a)*delta_j1u;
q1a=ans*R;
j2a=[R/2;R*cosd(-30);0];
j2u=[r*sind(-30+teta) -la/2;-r*cosd(-30+teta) -la*cosd(-30);0 1];
delta_j2u=j2u*inv(transpose(j2u)*j2u)*transpose(j2u)-eye(3);
inv(transpose(j2a)*delta_j2u*j2a)*transpose(j2a)*delta_j2u;
q2a=ans*R;
j3a=[R/2;R*cosd(210);0];
j3u=[r*sind(210+teta) -la/2;-r*cosd(210+teta) -la*cosd(210);0 1];
delta_j3u=j3u*inv(transpose(j3u)*j3u)*transpose(j3u)-eye(3);
inv(transpose(j3a)*delta_j3u*j3a)*transpose(j3a)*delta_j3u;
q3a=ans*R;
Forward_solution=[q1a ;q2a ;q3a]
A=inv(delta_j1u+delta_j2u+delta_j3u);
p1s=A*(delta_j1u*j1a);
p2s=A*(delta_j2u*j2a);
p3s=A*(delta_j3u*j3a);
Revers_solution=[p1s/R p2s/R p3s/R]