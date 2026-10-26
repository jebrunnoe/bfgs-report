# bfgs-report
Numerical results source code for "On the BFGS Method for Unconstrained Optimization"

## Running an example
To run a BFGS method on Rosenbrock's function:
- Set the method in `bfgs_parms.m`. For example, setting
  ```MATLAB
  parms.solver = bfgsR;
  parms.self_scaled = true;
  ```
  corresponds with the self-scaled factored Hessian method bfgsRS.
- Run the `runProb` script with no arguments:
  ```MATLAB
  [x, status, Stats] = runProb();
  ```
  This will run the configured solver on the function defined in `getCUSTomProblem.m`.
- The output from running `bfgsRS` on the Rosenbrock example is provided in `rosenbrock.txt` for reference.
