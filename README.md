# SQHquantum

% Implementation of the SQH Method for solving a quantum optimal control 
problem for 2 uncoupled spin 1/2 

The governing model has the structure:  y' = (A + u B ) y, y(0)=y_0

The cost functional with L2, L1, and discontinuous L1 costs. 
The weights of these costs are nu, gamma, and beta (threshold s). 
The target at final time is y_d

 Main references:

T. Breitenbach, A. Borzì,
A sequential quadratic Hamiltonian scheme for solving non-smooth 
quantum control problems with sparsity,
Journal of Computational and Applied Mathematics, 369 (2019) 112583.

