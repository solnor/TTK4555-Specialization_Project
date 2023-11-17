% Mechanical constants
J_m = 2.4059e-2;           % [kgm^2] - Moment of inertia of rotor
J_m_inv = 3000*1/J_m;           %
n_phases = 3;              % Motor winding phases
n_poles = 14;              % Number of poles in the stator
n_pole_pairs = n_poles/2;  % Number of pole pairs in the stator

% Electric constants
R_m   = 19.00e-3; % [Ohm] - Stator resistance
% R_m   = 0.04068; % [Ohm] - Stator resistance
L_m   = 8.97e-6;  % [H]   - Stator self inductance
% L_m   = 40.43e-6;  % [H]   - Stator self inductance
psi_m = 4.88e-3;  % [Wb]  - Magnetic flux linkage between rotor magnets/stator coils

% https://ieeexplore.ieee.org/mediastore_new/IEEE/content/media/63/7394209/6825903/liu.t4-2328655-large.gif
% R_m = 1.5e-2;
% L_m = 53e-6;
K_e = 0.052;
% n_pole_pairs = 7;


R_m_mat = R_m*eye(n_phases);
L_m_mat = L_m*eye(n_phases);
L_m_mat_inv = inv(L_m_mat);

fourier_coeff_a = [ 0.08333    0.00000     0.00000;
                    0.00000    1.05296    -1.05296;
                    0.07599    0.00000     0.00000;
                    0.00000    0.00000     0.00000;
                    0.05699    0.00000     0.00000;
                    0.00000   -0.04211     0.04211;
                    0.03377    0.00000     0.00000;
                    0.00000   -0.02148     0.02148;
                    0.01424    0.00000     0.00000;
                    0.00000    0.00000     0.00000;
                    0.00303    0.00000     0.00000;
                   -0.00067    0.00870    -0.00870;
                    0.00000    0.00000     0.00000;
                   -0.00048    0.00623    -0.00623;
                    0.00155    0.00000     0.00000;
                   -0.00270    0.00000     0.00000;
                    0.00356    0.00000     0.00000;
                   -0.00392   -0.00364     0.00364;
                    0.00375    0.00000     0.00000];

fourier_coeff_b = [-1.23020    0.60792     0.60792;
                    0.00000    0.00000     0.00000;
                   -0.30874   -0.27018    -0.27018;
                    0.04667    0.00000     0.00000;
                   -0.10013    0.02431     0.02431;
                    0.05305    0.00000     0.00000;
                   -0.02686   -0.01240    -0.01240;
                    0.04801    0.00000     0.00000;
                   -0.01285    0.03002     0.03002;
                    0.03709    0.00000     0.00000;
                   -0.02140   -0.00502    -0.00502;
                    0.02652    0.00000     0.00000;
                   -0.02988    0.00359     0.00359;
                    0.02005    0.00000     0.00000;
                   -0.02932   -0.01080    -0.01080;
                    0.01783    0.00000     0.00000;
                   -0.02187    0.00210     0.00210;
                    0.01768    0.00000     0.00000;
                   -0.01422   -0.00168    -0.00168];


% Clarke Transformation
Clarke = 2/3*[1 -1/2       -1/2;
              0 sqrt(3)/2  -sqrt(3)/2];