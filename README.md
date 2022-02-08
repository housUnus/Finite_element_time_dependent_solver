# Finite_element_time_dependent_solver
The "main.m" file contains the code of validation 

This code is written in MATLAB, and it has the ability to solve time-dependent linear diff equations.

The form of the equation that can be solved, the details of the solving can be found in the "rapport.pdf".

The commands to use will be found in the "commands.txt".

The project can use two type of shape functions P2 and P3, each one of these shape functions has a script that you can find in the project.

The "main.txt" contains the validation script, it uses time, space  validation and save the error for each one in two excel files.

The code contains comments to explain each part.	

** HOW TO USE THE SOLVER **

-cl1: the boundary condition row at the left side which contains row values, the first one is for the type of condition (0 for Dirichlet, 1 for Neumann) and the second value is for the value of the function at this point (u(0) for Dirichlet or u'(0) for Neumann)

-cl2: the boundary condition row at the right side which contains row values, the first one is for the type of condition (0 for Dirichlet, 1 for Neumann) and the second value is for the value of the function at this point (u(L) for Dirichlet or u'(L) for Neumann)

- Solver1D_P2_Ref_T_EDP_a(a,b,h,t0,tf,ht,u,mrf,gordre,alpha);

-- a  : the left side noude

-- b  : the right side noude

-- h  : the space step

-- t0 : the first moment

-- tf : the last moment 

-- ht : the time step

-- u  : the exacte solution

-- mrf : a parameter to choose between using a solver with a reference element technique (1) or just a normal solving (0)

-- gordre : a parameter to choose the order of gauss quadrature approximation 

-- alpha : a parameter of the EDP equation (alpha*du/dt) (1 or 0)

- sll = Solver1D_P2_Ref_T_EDP_a(0,4,0.01,0,1,.01,u,1,5,1): creation of an object with initialization of the necessary parameters.

- sll= sll.CLS(cl1,cl2) : store the boundry conditions on the object to use them during the post treatement.

- sll= sll.initialize(u0) : initialize the first columns of the approximated solution which contains the initial condition

- sll = sll.LoopTime() : start the post-treatment by starting a loop  through time and calculate the local matrices then global matrices, inserting the boundary conditions, then solving the system.

- sll = sll.erreur() : calculate the error for each point in space and time, they are stored in a variable called "er" of the "sll" object (sll.er)






