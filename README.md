# infection_network_reconstruction

You need the matlab software to solve optimization problems [CVX: Matlab Software for Disciplined Convex Programming](http://cvxr.com/cvx/)

The file figs_paper.m produces all the figures in the paper that are not schematics.

You should start with the file example_reconstruction.m, which of course have an example of reconstruction using the code.

The central functions are:
predator_prey_integrator.m, which integrates the dynamics, and
fun_net_recons.m, which takes measurements of the dynamics at discrete intervals and reconstructs the infection network M using cvx.
