=====================================================
              MATLAB 1D Time-Dependent PDE Solver
=====================================================

This project is written in MATLAB and provides a solver for 1D time-dependent linear equations. It supports two types of shape functions: P2 and P3. The solver includes validation scripts and the ability to perform both time and space convergence analysis.

-----------------------------------------------------
FILES
-----------------------------------------------------

- main.m         : Contains the main validation script.
- main.txt       : Contains a validation test using both time and space refinement. Errors are saved into Excel files.
- rapport.pdf    : Describes the equation being solved and provides details about the numerical method.
- commands.txt   : Lists the commands and parameters used to run the solver.
- [P2/P3 Scripts]: Located in the project folder, each script corresponds to a specific shape function.

-----------------------------------------------------
HOW TO USE THE SOLVER
-----------------------------------------------------

Boundary Conditions:

- cl1: Left-side boundary condition (vector of two values)
  - First value: 0 for Dirichlet, 1 for Neumann
  - Second value: Value of u(0) or u'(0)

- cl2: Right-side boundary condition (vector of two values)
  - First value: 0 for Dirichlet, 1 for Neumann
  - Second value: Value of u(L) or u'(L)

Solver Function:

    Solver1D_P2_Ref_T_EDP_a(a, b, h, t0, tf, ht, u, mrf, gordre, alpha)

Parameters:
    - a      : Left boundary node
    - b      : Right boundary node
    - h      : Space step size
    - t0     : Initial time
    - tf     : Final time
    - ht     : Time step size
    - u      : Exact solution (function handle)
    - mrf    : Use reference element technique (1) or not (0)
    - gordre : Gauss quadrature order
    - alpha  : Coefficient of the PDE (α * du/dt)

Example:

    sll = Solver1D_P2_Ref_T_EDP_a(0, 4, 0.01, 0, 1, 0.01, u, 1, 5, 1);

Object Methods:

    sll = sll.CLS(cl1, cl2)
        → Store boundary conditions for later use.

    sll = sll.initialize(u0)
        → Initialize the first time step with the initial condition.

    sll = sll.LoopTime()
        → Loop over time steps:
            - Assemble local and global matrices
            - Apply boundary conditions
            - Solve the system

    sll = sll.erreur()
        → Compute the error in space and time.
        → Result stored in: sll.er

-----------------------------------------------------
NOTES
-----------------------------------------------------

- The project includes comments in the MATLAB code to explain each section.
- Error data is saved in Excel format for both time and space validations.
- Make sure to define the exact solution function `u(t, x)` and initial condition `u0` before running the solver.

-----------------------------------------------------
