function [x,status,Stats] = bfgsTest( prob )
%        [x,status,Stats] = bfgsTest( prob ) finds a local minimizer
%        of a scalar-valued multivariate function f(x) using a
%        BFGS quasi-Newton method with a Wolfe linesearch.
%
%        This version dated 7-Dec-2019.
%        Philip E. Gill, University of California, San Diego.

  format compact

% -------------------------------
% Assign local control parameters
% -------------------------------
  parms       = bfgs_parms();

  solver      = parms.solver;      % Solver
  self_scaled = parms.self_scaled; % Self-scaled variant (true/false)
  trials      = parms.trials;      % recorded metric is averaged over #trials.
  write_file  = parms.write_file;  % write final run statistics to file
  out_file    = parms.out_file;    % write one line per iteration to file
% -----------------------
% End: control parameters
% -----------------------

% Initialize run statistics

  Stats      = zeros(1,2);
  jSkipped   = 1;          Skipped = 0;
  jTime      = 2;

  outfile    = fopen(strcat('./Statistics/', prob.name,'.out'),'w');

  if self_scaled suffix = 'S'; else suffix = ''; end

  if write_file
    str = strcat('./Statistics/', solver, suffix, '.stats');
    fid = fopen(str, 'a+');
  % fid     = fopen(strcat('./Statistics/Summary_', solver, '.csv'), 'a+');
  end

  trial_itns  = NaN(trials, 1);
  trial_nfs   = NaN(trials, 1);
  trial_times = NaN(trials, 1);

  for t = 1:trials
    tic
    switch(solver)
      case 'bfgsM'
      % Explicit approximate inverse Hessian
        [x,itn,nf,Skipped,status] = bfgsM( prob,outfile);
      case 'bfgsR'
      % BFGS with factored approximate Hessian (Recurrences)
        [x,itn,nf,Skipped,status] = bfgsR( prob,outfile);
      case 'bfgsC'
      % Conjugate factored BFGS
        [x,itn,nf,Skipped,status] = bfgsC( prob,outfile);
      case 'bfgsF'
      % Frobenius condition monitoring
        [x,itn,nf,Skipped,status] = bfgsF( prob,outfile);
      case 'bfgsS'
      % Cholesky factored inverse Hessian
        [x,itn,nf,Skipped,status] = bfgsS( prob,outfile);
      case 'bfgsZ'
      % Conjugate factored BFGS with auto-scaling
        [x,itn,nf,Skipped,status] = bfgsZ( prob,outfile);
      case 'bfgsL'
      % BFGS with LDL' factored approximate Hessian
        [x,itn,nf,Skipped,status] = bfgsL( prob,outfile);
      case 'bfgsN'
      % bfgsM with adaptive scaling
        [x,itn,nf,Skipped,status] = bfgsN(prob,outfile);
      case 'bfgsI'
      % bfgsM satisfying multi step condition
        [x,itn,nf,Skipped,status] = bfgsI(prob,outfile);
      case 'bfgsY'
      % bfgsH with extra interpolation condition
        [x,itn,nf,Skipped,status] = bfgsY(prob,outfile);
      case 'lrh'
        [x,itn,nf,Skipped,status] = lrh(prob,outfile);
      otherwise
        fprintf('ERROR: Unknown solver\n\n')
        fprintf(outfile, '%s\n', 'Unknown solver');
        fclose(outfile);
        x = []; status = []; stats = NaN;
        return
    end
    time = toc;
    trial_itns(t)  = itn;
    trial_nfs(t)   = nf;
    trial_times(t) = time;
    if status > 2
      trial_itns(:)  = itn;
      trial_nfs(:)   = nf;
      trial_times(:) = time;
      break;
    end
  end

  itn  = mean(trial_itns);
  nf   = mean(trial_nfs);
  Time = mean(trial_times);

%--------------------------------------------------------------------------
%  Collect statistics.
%    status = 0;  Optimal
%    status = 1;  Unbounded
%    status = 2;  Near Optimal
%    status = 3;  Steepest descent line search failure
%    status = 4;  Iteration limit
%--------------------------------------------------------------------------
  if itn > 0,
    pSkipped      = 100*Skipped/itn;
  else
    pSkipped      = 0;
  end
  Stats(jSkipped) = pSkipped;
  Stats(jTime  )  = Time;

  if       status == 0
    outcome = 'Optimal';
  elseif   status == 1
    outcome = 'Unbounded Objective';
  elseif   status == 2
    outcome = 'Near Optimal (badly scaled)';
  elseif   status == 3
    outcome = 'Too many iterations';
  elseif   status == 4
    outcome = 'Line search failure';
  end

  if status > 2
    nf   = NaN;  itn = NaN;
  end
  n     = prob.n;

  strP1 = sprintf('\n PROBLEM %-10s--- %-24s:', prob.name,outcome);
  strP2 = sprintf(' Itns  = %5g, Num f =%6g,' ,            itn,nf);
  strP3 = sprintf(' Cpu time   =%6.3f, n     = %4g\n',    Time, n);

  strS1 = sprintf('\n STATS    %-10s:',       prob.name);
  strS2 = sprintf(' Time  =%6.3f, Skipped =  %3i percent\n\n', Time, round(pSkipped) );
  strP  = [ strP1 strP2 strP3 ];   strS       = [ strS1 strS2 ];

  fprintf(strP); fprintf(outfile, strP);
  fprintf(strS); fprintf(outfile, strS);

  if write_file
    fprintf(fid, '%-10s%10d%10d%10d%10.3f  %-s\n', prob.name, n, itn, nf, Time, outcome);
  % fprintf(fid, '%s,%s,%s,%d,%d,%d,%f\n', strtrim(prob.name), ...
  %         solver, outcome, n, itn, nf, Time);
    fclose(fid);
  end

  fclose(outfile);
  if out_file == false
    delete(strcat('./Statistics/', prob.name, '.out'));
  end
end
