  Matlab Toolbox for Submodular Function Optimization (v 2.0)
  Tutorial script and all implementations Andreas Krause (krausea@gmail.com)
  Please read license.txt for license information and disclaimers.
 
  Tested in MATLAB 7.0.1 (R14), 7.2.0 (R2006a), 7.4.0 (R2007a, MAC), 7.9.0 (MAC)
 
  First some information on conventions:
  --------------------------------------
  All algorithms will use function objects.
  For example, to measure variance reduction in a Gaussian model, call 
    F = sfo_fn_varred(sigma,V)
  where sigma is the covariance matrix and V is the ground set, e.g., 1:size(sigma,1)
 
  Implemented algorithms:
  ----------------------- 
 
  Minimization:
 
  sfo_min_norm_point: Fujishige's minimum norm point algorithm for minimizing
        general submodular functions 
  sfo_queyranne: Queyranne's algorithm for minimizing symmetric submodular
        functions
  sfo_ssp: Submodular-supermodular procedure of Narasimhan & Bilmes for 
        minimizing the difference of two submodular functions
  sfo_s_t_min_cut: For solving min F(A) s.t. s in A, t not in A
  sfo_minbound: Return an online bound on the minimum solution
  sfo_greedy_splitting: Greedy splitting algorithm for clustering of 
        Zhao et al
 
  Maximization:
 
  sfo_polyhedrongreedy: For solving an LP over the submodular polytope
  sfo_greedy_lazy: The greedy algorithm for constrained maximization / coverage
  sfo_greedy_welfare: The greedy algorithm for solving allocation problems
  sfo_cover: Greedy coverage algorithm
  sfo_celf: The CELF algorithm for budgeted maximization
  sfo_ls_lazy: Local search algorithm for maximizing nonnegative submodular functions
  sfo_saturate: The Saturate algorithm of Krause et al. for robust optimization of submodular
        functions
  sfo_max_dca_lazy: The Data Correcting algorithm of Goldengorin et al. for 
        maximizing general (not necessarily nondecreasing) submodular functions
  sfo_maxbound: Return an online bound on the maximum solution
  sfo_pspiel: pSPIEL algorithm for trading off information and
        communication cost
  sfo_pspiel_orienteering: pSPIEL algorithm for submodular orienteering
  sfo_balance: eSPASS algorithm for simultaneous placement and balanced scheduling
 
  Miscellaneous
 
  sfo_lovaszext: Computes the Lovasz extension for a submodular function
  sfo_mi_cluster: Clustering algorithm using both maximization and
        minimization
  sfo_pspiel_get_path: Convert a tree into a path using the MST heuristic
        algorithm
  sfo_pspiel_get_cost: Compute the Steiner cost of a tree / path
 
  Submodular functions (also try 'help sfo_fn')
  
  sfo_fn_cutfun: Cut function
  sfo_fn_detect: Outbreak detection / facility location
  sfo_fn_entropy: Entropy of Gaussian random variables
  sfo_fn_infogain: Information gain about gaussian random variables
  sfo_fn_mi: Gaussian mutual information I(A; V\A)
  sfo_fn_varred: Gaussian variance reduction (orthogonal matching pursuit)
  sfo_fn_example: 2 element example from tutorial
  sfo_fn_iwata: Iwata's test function for testing minimization
  sfo_fn_ising: Energy function for Ising model for image denoising
  sfo_fn_residual: For defining residual submodular functions
  sfo_fn_invert: For defining F(A) = F'(V\A)-F(V)
  sfo_fn_lincomb: For defining linear combinations of submodular functions
 
  Example: sfo_tutorial
 
  If you use the toolbox for your research, please cite
  A. Krause. "SFO: A Toolbox for Submodular Function Optimization". Journal
    of Machine Learning Research (2010). 
 
  Here is an overview reference for submodularity in AI
  A. Krause, C. Guestrin. "Near-optimal Observation Selection Using Submodular Functions". 
    Survey paper, Proc. of 22nd Conference on Artificial Intelligence (AAAI) 2007 -- Nectar Track

Change log:
===============
Version 2.0
* Modified specification of optional parameters (using sfo_opt)
* Added sfo_ls_lazy for maximizing nonnegative submodular functions
* Added sfo_fn_infogain, sfo_fn_lincomb, sfo_fn_invert, ...
* Added additional documentation and more examples

Version 1.0
* added pSPIEL for informative path planning
* added eSPASS for simultaneous placement and scheduling
* new convention for submodular functions (incremental computations, etc.) Much faster!

Version 0.99:
* Initial release

