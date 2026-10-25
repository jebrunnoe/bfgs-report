function [parms] = bfgs_parms()
% ------------------------------------------------------------------------
% Authors : Philip E. Gill & Jeb Runnoe
% Date    : October 23, 2026
% Purpose : Sets the control parameters for variants of the bfgs method.
% ------------------------------------------------------------------------

% parms.solver        = 'bfgsM'   ;  % Solver
  parms.solver        = 'bfgsR'   ;  % Solver

% parms.self_scaled   = false     ;  % controls use of self-scaling.
  parms.self_scaled   = true      ;  % controls use of self-scaling.

  parms.maxIterations = 6000      ;  % maximum iterates allowed (Paper)

  parms.tiny          = 1.0e-12   ;  % tiny number
  parms.tolStat       = 1.0e-5    ;  % stationary tolerance     (Paper)

  parms.caution_level = 0         ;  % 0/1/2/3/4 caution level when accepting updates
  parms.updateTol     = 1.0e-6    ;  % Li-Fukushima Tolerance for cautious, 1.0e-16 for Powell
  parms.minCurv       = 0.8       ;  % Least positive curvature (mu)

  parms.jfMax         = 100       ;  % max functions per line search
  parms.gammaC        = 0.5       ;  % contraction factor in backtracking

  parms.searchType    = 'Wolfe'   ;  % Armijo/Wolfe     line search
  parms.wolfeTols     = zeros(1,2);
  parms.wolfeTols(1)  = 0.9       ;  % Wolfe gradient tolerance  (Paper)
  parms.wolfeTols(2)  = 1.0e-4    ;  % Wolfe function tolerance  (Paper)

  parms.condMax       = 1.0e+8    ;  % maximum condition estimate of R old version
  parms.condRMax      = 1.0e+8    ;  % maximum condition estimate of R

  parms.dxMax         = 10        ;  % limit on the line-search initial step
  parms.dxInf         = 1.0e+15   ;  % infinite change in x

  parms.ResetFreq     = NaN       ;  % Reset frequency
  parms.lines         = 20        ;  % lines between printing header
  parms.printSol      = 1         ;  % to print solution to screen
  parms.infinity      = 1.0e+15   ;  % user-defined +infinity
  parms.Assert        = 0         ;  % (0/1) check consistency
  parms.trials        = 1         ;  % Repeat each problem # times; metrics averaged
  parms.show_warnings = false     ;  % (0/1) print warnings
  parms.count_near    = true      ;  % count near optimal solutions
  parms.print_run_det = true      ;  % option to print run details
  parms.write_file    = true      ;  % write final run statistics to file
  parms.out_file      = false     ;  % write one line per iteration to file

  parms.unBoundedf    =-1.0e+10   ;  % definition of an unbounded objective value
  parms.debugItn      = 1         ;  % iteration for setting dbstop

end
