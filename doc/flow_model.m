##############################################
# INITIAL VALUES
##############################################
# Measurements in SI, but how knows, so let's see:
#   * length: in meters (m)
#   * density: in kg/m^3
#   * time: in seconds
#   * forces: in Newton (N)
#   * speed: in m/s

# CONSTANTS
##############################################
global window_height = 15; # in meters 
global surface_height= 5; # in meters

global N=20; # discretisation
global step_per_meter = window_height / N ;
global surface_level = floor(step_per_meter * surface_height);
# MAPS
##############################################

# A map discretisation:
##############################################
# (1, 1) ...... (1,N)
#  ...
# (surface_level, 1) ...., (surface_level, N) 
#
# (N, 1), ...., (N, N)

# (t) level
##############################################
# Grid
# The grid is a special grid, called staggered grid.
# u which is the x-axis direction of the velocity is given at every (i, j+1/2) points
# v which is the y-axis direction of the velocity is given at every (i +  1/2, j) points
# rho is given at (i + 1/2, j
global u_map = zeros(N, N);
global v_map = zeros(N,N);
global rho_map = zeros(N,N);
# (t-1) level
##############################################
global U_map = zeros(N, N);
global V_map = zeros(N,N);
global RHO_map = zeros(N,N);
global f_map = zeros(N,N);
##############################################
# FUNCTIONS 
##############################################

function init_globals()
	global window_height = 15 # in meters 
	global surface_height = 5 # in meters
	
	global N=20 # discretisation
	global step_per_meter = window_height / N 
	global surface_level = floor(step_per_meter * surface_height)
	
	global u_map = zeros(N, N);
	global v_map = zeros(N,N);
	global rho_map = zeros(N,N);
	global U_map = zeros(N, N);
	global V_map = zeros(N,N);
	global RHO_map = zeros(N,N);
	global f_map = zeros(N,N);
	u_map(:, end) = 2 * ones(N, 1);
	U_map = u_map;
	# v, and V can be zero everywhere
	v_map = zeros(N,N);
	V_map = v_map;
	# density of air, let approximate to 1 (1,2 at 20 degrees of Celsius)
	rho_map(1:surface_level - 1,1:N) = 1,2 * ones(surface_level - 1, N);
	# density of water, let assume 998 (20 degrees of Celsius)
	rho_map(surface_level:N, 1:N) = 998 * ones(N - surface_level + 1, N);
	RHO_map = rho_map;
	
	# Forces:
	# only gravity is counted
	f_map = 9.81 * ones(N,N);
endfunction
function next_iter()
	#globals
	global window_height; 
	global surface_height;
	global N;
	global step_per_meter;
	global surface_level;
	global u_map;
	global v_map;
	global rho_map;
	global U_map;
	global V_map;
	global RHO_map;
	global f_map;
	# locals
	u_next = u_map;
	v_next = v_next;
	rho_next = rho_next;
	
	
endfunction
##############################################
# MAIN
##############################################
init_globals()
rho_map
