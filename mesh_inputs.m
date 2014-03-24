function [h,l,m,n,option,E,nu,t,loadopt,V0,R,Mag,RadD,filen] = mesh_inputs 

h      = input('Height of rectangle ? ');
l      = input('Length of rectangle ? ');
m      = input('Increments along height ? ');
n      = input('Increments along length ? ');
option = input('0:Plane strain or 1:Plane stress ? ');
E      = input('Elasticity modulus ? ');
nu     = input('Poisson''s ratio ? ');
t      = input('Element thickness ? ');
disp('0: Parabolic shear stress or ');
disp('1: Linear bending stress at beam tip or');
disp('2: Beam forced into radius or');
disp('3: Beam transformed to 90 degree section loaded with pressure');
load_opt = input('4: Beam transformed to 90 degree section, inside surface displaced');
if load_opt == 0
    V0 = input('Magnitude of total shear force at beam tip ? ');
elseif load_opt == 1
    V0 = input('Magnitude of bending stress at top surface at beam tip ? ');
elseif load_opt==2
    R  = input('Enter radius that beam is deformed into. ');
elseif load_opt==3
    R = 2*l/pi;
    Mag = input('Internal pressure (in undeformed configuration) ? '); 
elseif load_opt==4
    R = 2*l/pi;
    RadD = input('Radial displacement ? ');
end
filen  = input('Write output to which file ? ','s');