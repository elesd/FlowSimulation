# FlowSimulation
2 dimensional flow simulation with Navier Stikes equations.

This is an ongoing project to create a simulator.

Plans for technology:
  * The implementation will be written in GLSL 4.0
  * For GUI imput it will use Qt4

Plans for the simulation:
  * A two dimensional sea map can be set in a Gui interface
  * Starting parameters will be given in the gui:
     * Reinolds number
     * Force of the flow
     * Pressaure

About the modell:
Assumptions:
 * The domain is related
 * At the beginning the pressaure, velocity is known.

Approximation:
Frjazinov approximation is used, which is based on the upwind (upstream shame), but it contains an correction to make it 2nd order.
