clc
clear all

% This MATLAB code calculates the stiffness and compliance matrices in global coordinate system for three types of materials 
% (Isotropic, Transversely Isotropic & Orthotropic) if the material properties are given in the material coordinate system. 


% Choosing material type
prompt="Choose Material type: Type 1 for Isotropic / Type 2 for Transversely Isotropic / Type 3 for Orthotropic";
i=1:3;
i=input(prompt);


    if (i==1)
% 1. Isotropic Material (Only 2 independent constants are required)
E = input('Enter the value of E:');
NU = input('Enter the value of NU:');
Theta = input('Enter the value of angle in degree unit:');

ISOTROPIC_Compliance_Matrix = [1/E -NU/E -NU/E 0 0 0; -NU/E 1/E -NU/E 0 0 0; -NU/E -NU/E 1/E 0 0 0; 0 0 0 (2*(1+NU))/E 0 0; 
    0 0 0 0 (2*(1+NU))/E 0; 0 0 0 0 0 (2*(1+NU))/E];
ISOTROPIC_Stiffness_Matrix = inv(ISOTROPIC_Compliance_Matrix); % inverse of the compliance matrix

% Engineering constants in the global coordinate system:
c = cos(Theta*0.017453333333); % converting from degree to radian unit
s = sin(Theta*0.017453333333);
% Transformed reduced compliance matrix is [S_bar] 
S11_bar = (1/E)*(c^4) - (2*(-NU/E)+((2*(1+NU))/E))*(s^2)*(c^2) + (1/E)*(s^4);
S12_bar = (-NU/E)*((s^4)+(c^4)) + ((1/E)+(1/E)-((2*(1+NU))/E))*(s^2)*(c^2);
S22_bar = (1/E)*(s^4) - (2*(-NU/E)+((2*(1+NU))/E))*(s^2)*(c^2) + (1/E)*(c^4);
S16_bar = (2*(1/E)+2*(-NU/E)-((2*(1+NU))/E))*s*(c^3) - (2*(1/E)+2*(-NU/E)-((2*(1+NU))/E))*s*(c^3)*(s^3)*c;
S26_bar = (2*(1/E)+2*(-NU/E)-((2*(1+NU))/E))*(s^3)*c - (2*(1/E)+2*(-NU/E)-((2*(1+NU))/E))*s*(c^3)*s*(c^3);
S66_bar = 2*(2*(1/E)+2*(1/E)+4*(-NU/E)-(2*(1+NU))/E)*(s^2)*(c^2) + ((2*(1+NU))/E)*(s^4)*(c^4);

Ex = 1/S11_bar; % Elastic modulus in x-direction
NUxy = -(S12_bar/S11_bar); % Poisson's ratio
Mx = -(S16_bar*E); % Non-dimensional shear coupling parameter
Ey = 1/S22_bar; % Elastic modulus in y-direction
NUyx = -(S12_bar/S22_bar); % Poisson's ratio
My = -(S26_bar*E); % Non-dimensional shear coupling parameter
Gxy = 1/S66_bar; % Shear modulus

% Output
disp('ISOTROPIC material_Compliance Matrix:');
disp(ISOTROPIC_Compliance_Matrix); %6×6 matrix
disp('ISOTROPIC material_Stiffness Matrix:');
disp(ISOTROPIC_Stiffness_Matrix); %6×6 matrix
disp('Elastic modulus in x-direction, Ex = ');
disp(Ex);
disp('NUxy = ');
disp(NUxy);
disp('Elastic modulus in y-direction, Ey = ');
disp('Mx = ');
disp(Mx);
disp(Ey);
disp('NUyx = ');
disp(NUyx);
disp('My = ');
disp(My);
disp('Shear modulus, Gxy = ');
disp(Gxy);
    end

    
    if (i==2)
% 2. Transversely Isotropic Material (5 independent constants are required)
E1 = input('Enter the value of E1:');
E2 = input('Enter the value of E2:');
NU12 = input('Enter the value of NU12:');
NU23 = input('Enter the value of NU23:');
G12 = input('Enter the value of G12:');
Theta = input('Enter the value of angle in degree unit:');

TRANSVERSELY_ISOTROPIC_Compliance_Matrix = [1/E1 -NU12/E1 -NU23/E2 0 0 0; -NU12/E1 1/E1 -NU23/E2 0 0 0; 
    -NU23/E1 -NU23/E1 1/E2 0 0 0; 0 0 0 1/(2*G12) 0 0; 0 0 0 0 1/(2*G12) 0; 0 0 0 0 0 (1+NU12)/E1];
TRANSVERSELY_ISOTROPIC_Stiffness_Matrix = inv(TRANSVERSELY_ISOTROPIC_Compliance_Matrix); % inverse of the compliance matrix

% Engineering constants in the global coordinate system:
c = cos(Theta*0.017453333333); % converting from degree to radian unit
s = sin(Theta*0.017453333333);
% Transformed reduced compliance matrix is [S_bar] 
S11_bar = (1/E1)*(c^4) - ((-2)*(NU12/E1)+((1+NU12)/E1))*(s^2)*(c^2) + (1/E1)*(s^4);
S12_bar = -(NU12/E1)*((s^4)+(c^4)) + ((1/E1)+(1/E1)-((1+NU12)/E1))*(s^2)*(c^2);
S22_bar = (1/E1)*(s^4) - ((-2)*(NU12/E1)+((1+NU12)/E1))*(s^2)*(c^2) + (1/E1)*(c^4);
S16_bar = (2*(1/E1)+2*(NU12/E1)-((1+NU12)/E1))*s*(c^3) - (2*(1/E1)+2*(NU12/E1)-((1+NU12)/E1))*s*(c^3)*(s^3)*c;
S26_bar = (2*(1/E1)+2*(NU12/E1)-((1+NU12)/E1))*(s^3)*c - (2*(1/E1)+2*(NU12/E1)-((1+NU12)/E1))*s*(c^3)*s*(c^3);
S66_bar = 2*(2*(1/E1)+2*(1/E1)+4*(-NU12/E1)-((1+NU12)/E1))*(s^2)*(c^2) + ((1+NU12)/E1)*(s^4)*(c^4);

Ex = 1/S11_bar; % Elastic modulus in x-direction
NUxy = -(S12_bar/S11_bar); % Poisson's ratio
Mx = -(S16_bar*E1); % Non-dimensional shear coupling parameter
Ey = 1/S22_bar; % Elastic modulus in y-direction
NUyx = -(S12_bar/S22_bar); % Poisson's ratio
My = -(S26_bar*E1); % Non-dimensional shear coupling parameter
Gxy = 1/S66_bar; % Shear modulus

% Output
disp('TRANSVERSELY ISOTROPIC material_Compliance Matrix:');
disp(TRANSVERSELY_ISOTROPIC_Compliance_Matrix); %6×6 matrix
disp('TRANSVERSELY ISOTROPIC material_Stiffness Matrix:');
disp(TRANSVERSELY_ISOTROPIC_Stiffness_Matrix); %6×6 matrix
disp('Elastic modulus in x-direction, Ex = ');
disp(Ex);
disp('NUxy = ');
disp(NUxy);
disp('Elastic modulus in y-direction, Ey = ');
disp('Mx = ');
disp(Mx);
disp(Ey);
disp('NUyx = ');
disp(NUyx);
disp('My = ');
disp(My);
disp('Shear modulus, Gxy = ');
disp(Gxy);
    end

        
      if (i==3)
% 3. Orthotropic Material (9 independent constants are required)
E1 = input('Enter the value of E1:');
E2 = input('Enter the value of E2:');
E3 = input('Enter the value of E3:');
NU12 = input('Enter the value of NU12:');
NU23 = input('Enter the value of NU23:');
NU13 = input('Enter the value of NU13:');
G12 = input('Enter the value of G12:');
G23 = input('Enter the value of G23:');
G13 = input('Enter the value of G13:');
Theta = input('Enter the value of angle in degree unit:');

ORTHOTROPIC_Compliance_Matrix = [1/E1 -NU12/E2 -NU13/E3 0 0 0; -NU12/E1 1/E2 -NU23/E3 0 0 0; -NU13/E1 -NU23/E2 1/E3 0 0 0; 
    0 0 0 1/(2*G23) 0 0; 0 0 0 0 1/(2*G13) 0; 0 0 0 0 0 1/(2*G12)];
ORTHOTROPIC_Stiffness_Matrix = inv(ORTHOTROPIC_Compliance_Matrix); % inverse of the compliance matrix

% Engineering constants in the global coordinate system:
c = cos(Theta*0.017453333333); % converting from degree to radian unit
s = sin(Theta*0.017453333333);
% Transformed reduced compliance matrix is [S_bar] 
S11_bar = (1/E1)*(c^4) - ((-2)*(NU12/E1)+((1+NU12)/E1))*(s^2)*(c^2) + (1/E1)*(s^4);
S12_bar = -(NU12/E1)*((s^4)+(c^4)) + ((1/E1)+(1/E1)-((1+NU12)/E1))*(s^2)*(c^2);
S22_bar = (1/E1)*(s^4) - ((-2)*(NU12/E1)+((1+NU12)/E1))*(s^2)*(c^2) + (1/E1)*(c^4);
S16_bar = (2*(1/E1)+2*(NU12/E1)-((1+NU12)/E1))*s*(c^3) - (2*(1/E1)+2*(NU12/E1)-((1+NU12)/E1))*s*(c^3)*(s^3)*c;
S26_bar = (2*(1/E1)+2*(NU12/E1)-((1+NU12)/E1))*(s^3)*c - (2*(1/E1)+2*(NU12/E1)-((1+NU12)/E1))*s*(c^3)*s*(c^3);
S66_bar = 2*(2*(1/E1)+2*(1/E1)+4*(-NU12/E1)-((1+NU12)/E1))*(s^2)*(c^2) + ((1+NU12)/E1)*(s^4)*(c^4);

Ex = 1/S11_bar; % Elastic modulus in x-direction
NUxy = -(S12_bar/S11_bar); % Poisson's ratio
Mx = -(S16_bar*E1); % Non-dimensional shear coupling parameter
Ey = 1/S22_bar; % Elastic modulus in y-direction
NUyx = -(S12_bar/S22_bar); % Poisson's ratio
My = -(S26_bar*E1); % Non-dimensional shear coupling parameter
Gxy = 1/S66_bar; % Shear modulus

% Output
disp('ORTHOTROPIC material_Compliance Matrix:');
disp(ORTHOTROPIC_Compliance_Matrix); %6×6 matrix
disp('ORTHOTROPIC material_Stiffness Matrix:');
disp(ORTHOTROPIC_Stiffness_Matrix); %6×6 matrix
disp('Elastic modulus in x-direction, Ex = ');
disp(Ex);
disp('NUxy = ');
disp(NUxy);
disp('Elastic modulus in y-direction, Ey = ');
disp('Mx = ');
disp(Mx);
disp(Ey);
disp('NUyx = ');
disp(NUyx);
disp('My = ');
disp(My);
disp('Shear modulus, Gxy = ');
disp(Gxy);
        end
        