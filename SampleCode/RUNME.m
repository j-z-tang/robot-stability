% This script invokes the various steps involved in the computation of the
% stability funnel outlined in the paper:
% "Invariant Funnels for Underactuated Dynamic Walking Robots: New Phase Variable and Experimental Validation"
% by Justin Tang

%% Check for required softwares - Spotless and MOSEK
display('Checking for required software')
try
    mosekopt;
catch err
    error('MOSEK NOT PRESENT. Please visit Mosek.com for installation')
end
try
    msspoly;
catch err
    error('Spotless NOT installed. Please visit https://github.com/RobotLocomotion/spotless-pod for installation')
end

%% Run step 0 - Writing functions containing walker dynamics
display('===== Beginning step 0 =====')
prompt = 'Re-compute dynamics of the walker?  If no, we will use default files previously computed.  Y/N [N]';
userInput = input(prompt,'s');
if strcmp(userInput,'Y')
    display('Re-computing dynamics files')
    step_0_compute_model_abs_coord
else
    display('Using default pre-computed files for dynamics')
end

%% Run step 1 - Compute Bezier polynomials for virtual constraint defining nominal trajectory
display('===== Beginning step 1 =====')
% prompt = 'Re-compute Bezier polynomial for virtual constraint defining nominal trajectory?  If no, we will use default files previously computed.  (Y/N):';
% userInput = input(prompt,'s');
% if strcmp(userInput,'Y')
%     display('Re-computing Bezier polynomials')
%     
% else
%     display('Using default pre-computed files for Bezier polynomials')
% end

display('Computing Bezier polynomial for virtual constraint defining nominal trajectory...')
step_1_virtual_constraint_compute_abs_coord

%% Run step 2 - Compute nominal open loop trajectory and nominal input
display('===== Beginning step 2 =====')
% prompt = 'Re-compute nominal open loop trajectory and nominal input?  If no, we will use default files previously computed.  (Y/N):';
% userInput = input(prompt,'s');
% if strcmp(userInput,'Y')
%     display('Re-computing')
%     step_2_1_closed_loop_full_f_sim
% else
%     display('Using default pre-computed files for nominal trajectory and nominal input')
% end

display('Computing  nominal open loop trajectory and nominal input...')
step_2_1_closed_loop_full_f_sim

%% Run step 3 - Compute transverse coordinate
display('===== Beginning step 3 =====')
display('Computing Transverse Coordinates')
step_3_transverse_dir

display('Checking and plotting computed transverse coordinates...')
check_alpha_3_phase_surface_no_dependence
drawnow
%% Run step 4 - Compute Riccati equation solution
display('===== Beginning step 4 =====')
% prompt = 'Re-Compute Riccati equation solution?  If no, we will use default files previously computed.  (Y/N):';
% userInput = input(prompt,'s');
% if strcmp(userInput,'Y')
%     display('Re-computing')
%     
% else
%     display('Using default pre-computed files for Riccati equation solution')
% end

display('Computing Riccati equation solution...')
step_4_compute_lyapunov

%% Run step 5 - compute taylor series expansion at each sample point
display('===== Beginning step 5 =====')
display('Computing taylor series expansion...')
step_5_compute_taylor_dynamics

%% Run step 6 - SOS verification
display('===== Beginning step 6 =====')
display('Computing Lyapunov certificates via SOS verifications...')
step_6_verification

display('Visualising funnel')
check_plotROA_polynomial_tau_proper_plot

prompt = 'Sample inside the funnel to ensure Lyapunov condition is strictly satisfied with the true dynamics? (Y/N):';
userInput = input(prompt,'s');
if strcmp(userInput,'Y')
    display('Sampling...')
    check_alpha_6
else
end

display('DONE!')
