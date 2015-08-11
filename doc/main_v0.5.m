# Script for testing Navier Stokes equation solver

1;

################################################################################
# Global variables:                                                            #
################################################################################

###############
# Parameters: #
###############

# Size of the time step
global T;

# Discretisation of the world
global H;

# Maximum number of iteration in the time line
global MAXITERATION;

# Width of the world
global WIDTH;

# Height of the world
global HEIGHT;

#########
# Maps: #
#########

# Map of the preassure 
# t -> [ x, y ]
global rhoMap;

# map of the velocity in the direction of X
# t -> [ x, y ]
global uMap;
# map of the velocity in the direction of Y
# t -> [ x, y ]
global vMap;

# Map of the forces
# [x, y] 
global fMap;
################################################################################
# Functions:                                                                   #
################################################################################

#####################
# Helper functions: #
#####################

# Traces the given map
# Params:
#   - map: map to trace
function traceMap(map, i)
  printf("time: %d\n", i)
  map{i}
endfunction

###############################################################################

# Get the value from a Map (Velocity U, V, or preassure Rho)
# Params:
#   - U: map for value calculation
#   - i: index for i :)
#   - i_offset: +-1 means j+-1/2
#   - j: index for j :)
#   - j_offset: +-1 means j+-1/2
function val = get_coord(U, i, i_offset, j, j_offset)
  val = U(2 * i + i_offset, 2 * j + j_offset);
endfunction

###############################################################################

# For D operators we need the velocity in the "preassure points". To calculate this
# we are using the average.
# This function has always a +1 offset on the "j" coordinate, therefore only the 
# i offset is relevant.
# Params:
#   - U: matrix of the velocity
#   - i: integer part of the index
#   - i_offset: +-1 means i+-1/2
#   - j: integer part of the index
#   - j_offset: +-1 means j+-1/2
function val = u_x(U, i, i_offset, j)
  val = 1/2 * (get_coord(U, i + i_offset, 0, j, 0) + get_coord(U, i, 0, j, 0));
endfunction

###############################################################################

# For D operators we need the velocity in integer points. To calculate this
# we are using the average.
# Params:
#   - U: matrix of the velocity
#   - i: integer part of the index
#   - j: integer part of the index
function val = u_y(U, i, j)
  val = 1/2 * (get_coord(U, i, 0, j, 1) + get_coord(U, i, 0, j - 1, 1));
endfunction

###############################################################################

# For D operators we need the velocity in the "preassure points". To calculate this
# we are using the average.
# Params:
#   - V: matrix of the velocity
#   - i: integer part of the index
#   - j: integer part of the index
function val = v_x(V, i, j)
  val = 1/2 * (get_coord(V, i, 1, j, 0) + get_coord(V, i - 1, 1, j, 0));
endfunction

###############################################################################

# For D operators we need the velocity in the "preassure points". To calculate this
# we are using the average.
# Params:
#   - V: matrix of the velocity
#   - i: integer part of the index
#   - j: integer part of the index
function val = v_y(V, i, j)
  val = 1/2 * (get_coord(V, i, 1, j + 1, 0) + get_coord(V, i, 1, j, 0));
endfunction

###############################
# Functions about derivation: #
###############################

# Calculates the gradient of the preasure.
# The preasure is given in the "half points" therefore in optimal case
# extra calculation does not necessary. 
# Params: 
#   - Rho: Map if the gradient.
#   - i, j: point to get the gradient
# Returns: The gradient of the vec.
function gradient = grad_rho(Rho, i, j)
  global H;
  gradient(1) = 1/H * (get_coord(Rho, i, 1, j, 1) - get_coord(Rho, i - 1, 1, j , 1));
  gradient(2) = 1/H * (get_coord(Rho, i, 1, j, 1) - get_coord(Rho, i, 1, j - 1, 1));
endfunction

#####################################
# Functions about stability (zeta): #
#####################################

# Helper function for stability calculation. On the border we apply the stability
# value with a factor of 2/3. This is necessery to keep the approximation in second
# order.
function val = gamma(i)
  global WIDTH;
  global H;
  if i == 1 || i == idivide(WIDTH, H, "fix")
    val = 2/3;
  else
    val = 1;
  endif
endfunction

###############################################################################

# Support value for the approximation operator. It garanties the self adjuction
# and positive semi definit properties
# Params:
#   - u: current value of velocity
#   - i: current place
function val = zeta(u, i)
  val = 1/2 * sqrt(gamma(i) * H^3 * abs(u));
endfunction

###############################################################################

# Zeta function derivation.
# Params:
#   - U: map for value calculation
#   - i: current position on the map
#   - j: current position on the map
function val = zeta_uu_x(U, i, j)
  global H;
  val =
    1/(H^2) * (zeta(u_x(U, i + 1, 1, j, 1)) * (get_coord(U, i + 2, 0, j, 1) - get_coord(U, i + 1, 0, j, 1)) / H 
               - 2 * zeta(u_x(U, i, 1, j, 1)) * (get_coord(U, i + 1, 0, j, 1) - get_coord(U, i, 0, j, 1)) / H 
               + zeta(u_x(U, i - 1, 1, j, 1)) * (get_coord(U, i, 0, j, 1) - get_coord(U, i - 1, 0, j, 1)) / H);
endfunction

###############################################################################

# Zeta function derivation.
# Params:
#   - U: map for value calculation
#   - i: current position on the map
#   - j: current position on the map
function val = zeta_vu_y(U, V, i, j)
  global H;
  1/(H^2) * (zeta(v_x(V, i, j + 1)) * (get_coord(U, i, 0, j + 1, 1) - get_coord(U, i, 0, j, 1)) / H
             - 2 * zeta(v_x(V, i, j)) * (get_coord(U, i, 0, j , 1) - get_coord(U, i, 0, j - 1, 1)) / H
             + zeta(v_x(V, i, j - 1)) * (get_coord(U, i, 0, j - 1, 1) - get_coord(U, i, 0, j - 2, 1)) / H);
endfunction

###############################################################################

# Zeta function derivation.
# Params:
#   - U: map for value calculation
#   - i: current position on the map
#   - j: current position on the map
function ret = zeta_uv_x(U, V, i, j)
  global H;
  ret = 1/(H^2) * (zeta(u_y(U, i + 1, j), i) * (get_coord(V, i + 1, 1, j, 0) - get_coord(V, i, 1, j, 0)) / H
                   - 2 * zeta(u_y(U, i, j), i) * (get_coord(V, i, 1, j, 0) - get_coord(V, i - 1, 1, j, 0)) / H
                   + zeta(u_y(U, i - 1, j), i) * (get_coord(V, i - 1, 1, j, 0) - get_coord(V, i - 2, j, 0)) / H);
endfunction

###############################################################################

# Zeta function derivation.
# Params:
#   - U: map for value calculation
#   - i: current position on the map
#   - j: current position on the map
function ret = zeta_vv_y(V, i, j)
  global H;
  ret = 1/(H^2) * (zeta(v_y(V, i, j + 1), i) * (get_coord(V, i, 1, j + 2, 0) - get_coord(V, i, 1, j + 1, 0)) / H
                   - 2 * zeta(v_y(V, i, j), i) * (get_coord(V, i, 1, j + 1, 0) - get_coord(V, i, 1, j, 0)) / H
                   + zeta(v_y(V, i, j - 1), i) * (get_coord(V, i, 1, j, 0) - get_coord(V, i, 1, j - 1, 0)) / H);
endfunction

###########################################
# Functions about Operator decomposition: #
###########################################
#   The L operator is consist of two part #
#   basicaly, a D operator, and a         #
#   Raynolds number part. These equations #
#   are different by direction.           #
###########################################

# Part of the opertator decomposition
# Params:
#   - U: Map of value calculation
#   - i: current position on the map
#   - j: current position on the map
function dx = Duu(U, i, j)
  global H;
  dx = 1 / H * (u_x(U, i, 1, j, 1)^2
                + zeta(u_x(U, i, 1, j, 1), i) 
                  * zeta_uu_x(U, i, j)) # i + 1/2
       - 1 / H * (u_x(U, i, -1, j, 1)^2
                  + zeta(u_x(U, i, -1, j, 1), i - 1) 
                    * zeta_uu_x(U, i - 1, j)); # i - 1/2
endfunction

###############################################################################

# Part of the opertator decomposition
# Params:
#   - U: Map of value calculation
#   - V: Map of value calculation
#   - i: current position on the map
#   - j: current position on the map
function dy = Dvu(U, V, i, j)
  global H;
  dy = 1 / H * (v_x(V, i, j + 1) * u_y(U, i, j + 1) 
                + zeta(v_x(U, i, j + 1), i) * zeta_vu_y(U, V, i, j + 1)) # j + 1
       - 1 / H * (v_x(V, i, j) * u_y(U, i, j) 
                  + zeta(v_x(U, i, j), i) * zeta_vu_y(U, V, i, j)) # j
endfunction

###############################################################################

# Part of the opertator decomposition
# Params:
#   - U: Map of value calculation
#   - i: current position on the map
#   - j: current position on the map
function ret = Duv(U, V, i, j)
  global H;
  ret = 1 / H * (u_y(U, i + 1, 0, j, 0) * v_x(V, i + 1, 0)
                 + zeta(u_y(U, i + 1, 0, j, 0), i + 1) 
                   * zeta_uv_x(U, V, i + 1, j)) # i + 1
        - 1 / H * (u_y(U, i, 0, j, 0) * v_x(V, i, 0)
                   + zeta(u_y(U, i, 0, j, 0), i) 
                     * zeta_uv_y(U, V, i, j)) # i
endfunction

###############################################################################

# Part of the opertator decomposition
# Params:
#   - U: Map of value calculation
#   - i: current position on the map
#   - j: current position on the map
function ret = Dvv(V, i, j)
  global H;
  ret = 1 / H * (v_y(V, i, 1, j, 1)^2
                 + zeta(v_y(U, i, 1, j, 1), i) 
                   * zeta_vv_y(V, i, j)) # j + 1/2
        - 1 / H * (v_y(V, i, 1, j - 1, 1)^2
                   + zeta(v_y(U, i, 1, j - 1, 1), i) 
                     * zeta_vv_y(V, i, j - 1)) # j + 1/2
endfunction

###############################################################################

# Convective part of the L operator. This function calculates it in the direction
# of x.
#   - U: Map of value calculation
#   - V: Map of value calculation
#   - i: current position on the map
#   - j: current position on the map
function val = D_x(U, V, i, j)
  val = Duu(U, i, j) + Dvu(U, V, i, j)
endfunction

###############################################################################

# Convective part of the L operator. This function calculates it in the direction
# of y.
#   - U: Map of value calculation
#   - V: Map of value calculation
#   - i: current position on the map
#   - j: current position on the map
function val = D_y(U, V, i, j)
  val = Duv(U, V, i, j) + Dvv(V, i, j)
endfunction

###############################################################################

# 1/RE * div(u) part of the operator.
#   - U: Map of value calculation
#   - i: current position on the map
#   - j: current position on the map
function res = div_u_part(U, i, j)
  global RE;
  global H;
  res = 
    1/RE * ((get_coord(U, i + 1, 0, j, 1) - 2 * get_coord(U, i, 0, j, 1) + get_coord(U, i - 1, 0, j, 1))/(H^2)
            + (get_coord(U, i, 0, j + 1, 1) - 2 * get_coord(U, i, 0, j, 1) + get_coord(U, i, 0, j - 1, 1))/(H^2));
endfunction

###############################################################################

# 1/RE * div(v) part of the operator.
#   - V: Map of value calculation
#   - i: current position on the map
#   - j: current position on the map
function res = div_v_part(V, i, j)
  global RE;
  global H;
  res = 
    1/RE * ((get_coord(V, i + 1, 1, j, 0) - 2 * get_coord(V, i, 1, j, 0) + get_coord(V, i - 1, 1, j, 0))/(H^2)
            + (get_coord(V, i, 1, j + 1, 0) - 2 * get_coord(V, i, 1, j, 0) + get_coord(V, i, 1, j - 1, 0))/(H^2));
endfunction

###############################################################################

# Calculates the core operator in the direction of x (it is the same direction
# as the v velocity)
# Params:
#   - U, V: vector-map to apply the operator
#   - i: x coordinate
#   - j: y coordinate
function res = L_y(U, V, i, j)
  res = D_y(U, V, i, j) - div_v_part(V, i, j);
endfunction

###############################################################################

# Calculates the core operator in the direction of x (it is the same direction
# as the u velocity)
# Params:
#   - U, V: vector-map to apply the operator
#   - i: x coordinate
#   - j: y coordinate
function res = L_x(U, V, i, j)
  res = D_x(U, V, i, j) - div_u_part(U, i, j);
endfunction

###################
# Main functions: #
###################

# Initialize global maps
function initMaps()
# TODO create fully calculated stredger grid. Efficiency is not part of the goal of this script.
#      therefore it is easier to store it once, and use it.
endfunction;

###############################################################################

# Do the predictor steps in the ith time frame.
# Params:
#   - t: current time frame
# Returns: w map [ x, y ] which is a solution wich not considers 
#          the Mess Equality.
function [w, ww] = predictorStep(t)
  global rhoMap;
  global uMap;
  global vMap;
  global fMap;
  global HEIGHT;
  global WIDTH;
  global H;
  w = zeros(HEIGHT / H, WIDTH / H);
  for i = 1 : idivide(HEIGHT, H, "fix")
    for j = 1 : idivide(WIDTH, H, "fix")
      if i != 1 && j != 1
        w(i, j) = fMap(i, j)(1) - grad_rho(rhoMap{t - 1}, i, j)(1) + (1 / T * uMap{t - 1}(i, j));
        w(i, j) = w(i, j) / (1 / T + L_x(uMap{t - 1}, vMap{t - 1}, i, j));
        ww(i, j) = fMap(i, j)(2) - grad_rho(rhoMap{t - 1}, i, j)(2) + (1 / T * uMap{t - 1}(i, j));
        ww(i, j) = ww(i, j) / (1 / T + L_y(uMap{t - 1}, vMap{t - 1}, i, j));
      else
        w(i, j) = fMap(i, j);
        ww(i, j) = fMap(i, j);
      endif
    endfor
  endfor     
endfunction

###############################################################################

# Solves the Poisson Equation to calculate the preassure changes between two
# time frame.
# Params: 
#   - i: current time frame
#   - w: result of the predictor step
# Returns: q map [ x, y ] which is the preassure delta.
function q = solvePoission(w, v, i)
#TODO
endfunction

###############################################################################

# Calculates the preassure in the ith time frame
# Params: 
#   - i: current time frame
#   - q: result of the predictor step
function calculatePreassure(q, i)
  global rhoMap;
#TODO
  rhoMap{i} = rhoMap{i - 1} + q;
endfunction

###############################################################################

# Calculates the velocity maps in the ith time frame
# Params: 
#   - i: current time frame
#   - w: result of the predictor step
#   - q: delta of the preassure
function calculateVelocity(w, v, q, i)
  global T;
#TODO
  uMap{i} = w - T * grad(q);
  vMap{i} = v - T * grad(q);
endfunction

###############################################################################

# Traces all the importent/global maps
# Params:
#   - i: current time frame
function traceMaps(i)
  global uMap;
  global vMap;
  global rhoMap;
  traceMap(uMap, i);
  traceMap(vMap, i);
  traceMap(rhoMap, i);
endfunction

################################################################################
# Script main                                                                  #
################################################################################

initMaps();

for i = 1:MAXITERATION
  [w, v] = predictorStep(i);
 # q = solvePoission(w, v, i);
 # calculatePreassure(q, i);
 # calculateVelocity(w, v, q, i);
 # traceMaps();
endfor







