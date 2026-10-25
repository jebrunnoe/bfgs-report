function [x,itn,nf,skipped,status] = bfgsM(prob,outfile)
%==========================================================================
%[x,itn,nf,skipped,status] = bfgsM(prob,outfile)
%  finds a local minimizer of a scalar-valued multivariate function using
%  a variant of the BFGS quasi-Newton method with a Wolfe line search.
%
%  prob is a structure that provides the details of the
%  CUTEst problem being solved.
%
%  On successful termination, x is the first point found at
%  which the norm of the gradient is less than sqrt(eps).
%
%  Input arguments are set in the file bfgs_parms.m.
%
%  The meaning of each line of output is as follows
%
%  Itn          The current iteration number.
%
%  Nf           The cumulative number of times that the
%               objective function f(x) has been evaluated.
%
%  Step         The steplength taken along the search direction.
%
%  f(x)         The current value of the objective function f(x).
%
%  Norm g       The current value of the gradient of f(x). In
%               the neighborhood of a local minimizer, Norm g
%               will be small.
%
%  cond H       condition estimate of the approximate Hessian
%==========================================================================

%  This version dated 29 Jun 2022.
%  Philip E. Gill & Jeb Runnoe, University of California, San Diego.
% -------------------------------------------------------------------------
  format compact

% -------------------------------
% Assign local control parameters
% -------------------------------
  parms         = bfgs_parms();

  solver        = parms.solver       ;  % Solver
  self_scaled   = parms.self_scaled  ;  % Self-scaled variant (true/false)
  maxIterations = parms.maxIterations;  % maximum iterates allowed
  tolStny       = parms.tolStat      ;  % stationarity   tolerance
  jfMax         = parms.jfMax        ;  % max functions per line search
  wolfeTols     = parms.wolfeTols    ;  % Wolfe function tolerance
  caution_level = parms.caution_level;  % caution level when accepting updates
  dxMax         = parms.dxMax        ;  % limit on the line-search initial step
  dxInf         = parms.dxInf        ;  % infinite change in x
  minCurv       = parms.minCurv      ;  % Least positive curvature (mu)
  updateTol     = parms.updateTol    ;  % Tolerance for rejecting BFGS update
  lines         = parms.lines        ;  % lines between printing header
  ResetFreq     = parms.ResetFreq    ;  % Reset frequency
  infinity      = parms.infinity     ;  % user-defined +Infinity
  Assert        = parms.Assert       ;  % (0/1) check consistency
  unBoundedf    = parms.unBoundedf   ;  % unbounded objective value
  debugItn      = parms.debugItn     ;  % iteration for setting dbstop
  searchType    = parms.searchType   ;  % line-search type
  count_near    = parms.count_near   ;  % Option to count near-optimal solution
  print_run_det = parms.print_run_det;  % Option to print iteration details
% -----------------------
% End: control parameters
% -----------------------
  if  print_run_det
    if  self_scaled
      ssYN = 'Yes';
    else
      ssYN = 'No';
    end

    if  count_near
      count_nearYN = 'Yes';
    else
      count_nearYN = 'No';
    end

    if     caution_level == 0                % Incautious
      caution = '    None';
    elseif caution_level == 1                % Li-Fukushima rule 1
      caution = '    LF 1';
    elseif caution_level == 2                % Li-Fukushima rule 2
      caution = '    LF 2';
    elseif caution_level == 3                % Powell
      caution = '  Powell';
      %          12345678
    elseif caution_level == 4                % Powell proof
      caution = 'Powell 2';
    end

    fprintf('\nSolver');
    fprintf('\n=======\n');
    fprintf(' Solver type.......... %14s%7s',solver,                  ' ');
    fprintf(' Self-scaled updates.. %14s%7s',ssYN,                    ' ');
    fprintf(' Line search.......... %14s\n' ,searchType                  );

    fprintf(' Max iterations....... %6s%8g%7s'  , ' ', maxIterations, ' ');
    fprintf(' Stationarity tol..... %5s%9.2e%7s', ' ', tolStny,       ' ');
    fprintf(' Count near-optimal... %14s\n'     ,      count_nearYN      );

    fprintf(' Wolfe tolerance...... %5s%9.2e%7s', ' ', wolfeTols(1),  ' ');
    fprintf(' Armijo tolerance..... %5s%9.2e%7s', ' ', wolfeTols(2),  ' ');
    fprintf(' BFGS acceptance test. %14s\n'     ,      caution           );

    fprintf(' Initial search step.. %5s%9.2e%7s', ' ', dxMax,         ' ');
    fprintf(' Maximum search step.. %5s%9.2e\n' , ' ', dxInf             );
  end

% Initialize the counters

  nf       = 1;                 % number of function evaluations
  skipped  = 0;

  if isnan(ResetFreq), ResetFreq = maxIterations + 1; end
  if isnan(debugItn ), debugItn  = maxIterations + 1; end

  QN       =  0;                SD      = 1;

  sqrteps  = eps^(1/2);         cbrteps =  eps^(1/3);
  stepMax  = 10^5;
  MjjMax   = 1.0d+4;            MjjMin  = 1.0d-2;

  x        = prob.x;
  n        = prob.n;
  [f,g]    = prob.obj(x);

  normx    = norm(x);
  normg    = norm(g);         normgMax = 1 + normg;
  fMax     = 1 + abs(f);

  Mjj      = 1/max(1,normg);      Mjj   =  min(max(Mjj,MjjMin),MjjMax);
  M        = diag(Mjj*ones(1,n));

  Info     = '     ';
% header   = '  Itn    Nf    Step     jf    Objective      Norm g   cond H';
  header   = '  Itn    Nf    Step     jf    Objective      Norm g';

  step     =  0;                  iExit = 1; % First line of output
  dType    = SD;

  Optimal  = false;   NearOptimal = false;   Unbounded = false;
  ItnLimit = false;   BadSearch   = false;

  itn      = 0;
% -------------------------------------------------------------------------
% Main loop.
% -------------------------------------------------------------------------
  while (true)
    if      normg    <= tolStny*(1 + abs(f)) || normg <= sqrteps
      Optimal     = true;
    elseif  f        <= unBoundedf
      Unbounded   = true;
    elseif count_near && BadSearch && (normg  <= sqrteps*normgMax ||  ...
                                       abs(f) <= cbrteps*cbrteps*fMax)
      NearOptimal = true;
    end

    ItnLimit  = itn > maxIterations - 1;

    if BadSearch && dType == QN
      Mjj     = 1/max(1,normg);    Mjj   =  min(max(Mjj,MjjMin),MjjMax);
      M       = diag(Mjj*ones(1,n));
      dType   = SD;              Info(5:5) = 'R';
      BadSearch = false;
    end

  %------------------------------------------------------------------------
  % Print iteration details.
  %------------------------------------------------------------------------
    if print_run_det
      nLine    = rem( itn, lines );
      if nLine == 0
        fprintf('\n%s\n', header);  fprintf(outfile, '\n%s\n', header);
      end

      if itn == 0
        cstep = '    --    ';
        cjf   =     '-- ';
      else
        cstep = sprintf ( ' %9.2e', step );
        cjf   = sprintf ( '%2g '  ,   jf );
      end
      if iExit ~= 1
        cExit = sprintf ( '(%1g)', iExit );
      else
        cExit = '   ';
      end
      str1  = sprintf ( '%5g %5g',        itn,      nf );
      str2  = sprintf ( '%s%s%s',               cstep,  cExit, cjf );
      str3  = sprintf ( '%14.7e %9.2e', f, norm(g) );
      str4  = Info;
      str   = [ str1 str2 str3 str4 ];
      fprintf('%s\n', str);    fprintf(outfile, '%s\n', str);
      Info  = '     ';  % Info = ' nMlR'
    end

    if  Optimal || Unbounded || NearOptimal || ItnLimit || BadSearch
      break;
    end

    p         = -M*g;
    gtp       =  g'*p;     normp = norm(p);

    if itn == 0
      stepPrev = 10^(20);
    end

    stepMax     = dxInf/(1+normp);
    if (1 + normx)*dxMax >= stepMax*normp
      stepLimit = stepMax;
    else
      stepLimit = (1 + normx)*dxMax/normp;
    end

%   step        = min([1 stepPrev stepLimit stepMax]);
    step        = min([1 stepLimit stepMax]);
    if step    == stepMax, Info(4:4) = 'l'; end

    while (true)
      if     strcmp(searchType,'Wolfe' )
        [step,xnew,fnew,gnew,jf,iExit] ...
             = wolfeLSfgUC(prob,wolfeTols,jfMax,step,stepMax,f,g,p,x);
      elseif strcmp(searchType,'Armijo')
        [step,xnew,fnew,gnew,jf,iExit] ...
             = armijoLS(prob,step,stepMax,f,g,p,x);
      end
      if (sum(isnan([fnew; gnew])) + sum(isinf([fnew; gnew]))) == 0, break; end
      step   = step/10; stepMax = step;
    end

    nf       = nf + jf;
    stepPrev = step;

    if  iExit > 3
      BadSearch = true; % Line search failure
    else

    % BFGS update.

      d      =  step*p;
      y      =  gnew - g;
      ytd    =  y'*d;
      My     =  M*y;
      ytMy   =  y'*My;
      sigma  =  ytd/(ytMy);

      gamma  = 1;
      if self_scaled

      % Self-scaled update

        Hd    = -step*g;
        dtHd  =  d'*Hd;
        gamma =  dtHd/ytd;
      end

      normg   =  norm(gnew);

      ResetH  = itn > 0  &&  mod( itn, ResetFreq ) == 0;

      if ResetH
        Mjj   = sigma/max(1,normg);         Mjj =  min(max(Mjj,MjjMin),MjjMax);
        M     = diag(Mjj*ones(1,n));
        dType = SD;                   Info(5:5) = 'R';
      else

        if     caution_level == 0                % Incautious
          Update = ytd > 0;
        elseif caution_level == 1                % Li-Fukushima rule 1
          if normg >= 1
            c = 1.0e-2;
          else
            c = 3;
          end
          Update = ytd >  updateTol*(d'*d)*normg^c;
        elseif caution_level == 2                % Li-Fukushima rule 2
          Update = ytd >  updateTol*(d'*d)*normg;
        elseif caution_level == 3                % Powell
          Update = ytd >  updateTol*(d'*d);
        elseif caution_level == 4                % Powell proof
          Update = ytd >  updateTol*(y'*y);
        else
          fprintf('ERROR: caution level unrecognized.\n\n')
          return
        end

        if Update
        %--------------------------------------------------
        % Update the inverse Hessian approximation.
        %--------------------------------------------------
          sd   = (1/ytd)*d;        sMy = (1/ytMy)*My;
          w    = sd - sMy;
          M    = gamma*(M - sMy*My' + (ytMy*w)*w') + sd*d';
        else
          Info(2:2) = 'n';
          skipped   = skipped + 1;
        end
      end
      x     = xnew;  f = fnew;  g = gnew;
      normx = norm(x);
    end

    itn    = itn + 1;
    if itn >= debugItn,
      jj   = 1; % Set break here for debugging
    end
  end

  if      Optimal
    status = 0;
  elseif  Unbounded
    status = 1;
  elseif  NearOptimal
    status = 2;
  elseif  ItnLimit
    status = 3;
  elseif  BadSearch
    status = 4;
  end
end
