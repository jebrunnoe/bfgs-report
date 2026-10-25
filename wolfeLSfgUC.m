function [alfa,x,objf,objg,numf,iExit] ...
         = wolfeLSfgUC(prob,wolfeTols,jfMax,alfa0,alfaMax,objf0,objg0,p,x0)
%==========================================================================
% [alfa,x,objf,objg,numf,iExit]
%            = wolfeLSfgUC(prob,wolfeTols,alfa0,alfaMax,objf0,objg0,p,x0)
%
%  wolfeLSfg  finds a step along the search direction  p,  such that
%  the function  f  is sufficiently reduced, i.e.,
%            f(x0 + alfa*p)  <  f(x).
%
%  On entry,
%  alfa0         is the initial estimate of the step.
%  alfaMax       is the maximum allowable step.
%  wolfeTols(1)  Line search accuracy parameter in the range (0,1).
%                0.001 means very accurate.   0.99 means quite sloppy.
%  objf0,objg0   are the function and gradient of the objective at x.
%  x0            is the base point for the search.
%  p             is the search direction.
%
%  Output parameters...
%
%  alfa          is the final step.
%  x             is the new point  x0 + alfa*p.
%  objf          is the function value at x.
%  objg          is the gradient at x.
%
%  iExit    Result
%  -----    ------
%    1      The search is successful and step < alfaMax.
%    2      The search is successful and step = alfaMax.
%    3      A better point was found but too many functions
%           were needed (no sufficient decrease obtained).
%    4      alfaMax < tolabs (too small to do a search).
%    5      is not used
%
%    6      No useful step.
%           The interval of uncertainty is less than 2*tolabs.
%           The minimizer is very close to step = zero
%           or the gradients are not sufficiently accurate.
%    7      Too many function calls.
%    8      Bad input parameters
%           (alfaMax le toltny  or  objg'p ge 0).
%==========================================================================

%
%  Local variables
%  ---------------
%
%  Bracketed     is false if f and f' have not been evaluated at
%                the far end of the interval of uncertainty.  In this
%                case, the point b will be at alfaMax + tol(alfaMax).
%
%  Cramped       is true if alfaMax is very small (le tolabs).  If the
%                search fails, this indicates that a zero step should
%                be taken.
%
%  Extrap        is true if xw lies outside the interval of
%                uncertainty.  In this case, extra safeguards are
%                applied to allow for instability in the polynomial
%                fit.
%
%  Moved         is true if a better point has been found, i.e.,
%                alfbst gt 0.
%
%  Wset          records whether a second-best point has been
%                determined it will always be true when convergence
%                is tested.
%
%  nSamea        is the number of consecutive times that the
%                left-hand end point of the interval of uncertainty
%                has remained the same.
%
%  nSameb        similarly for the right-hand end.
%
%  a,b,alfaBest  define the current interval of uncertainty.
%                A minimizer lies somewhere in the interval
%                [alfaBest + a, alfaBest + b].
%
%  alfaBest      is the best point so far.  It is always at one end
%                of the interval of uncertainty.  hence we have
%                either  a lt 0,  b = 0  or  a = 0,  b gt 0.
%
%  fbest, gbest  are the values of f, f' at the point alfbst.
%
%  factor        controls the rate at which extrapolated estimates
%                of alfa may expand into the interval of uncertainty.
%                factor is not used if a minimizer has been bracketed
%                (i.e., when the variable Bracketed is true).
%
%  fw, gw        are the values of f, f' at the point alfaBest + xw.
%                they are not defined until Wset is true.
%
%  xTry          is the trial point in the shifted interval (a, b).
%
%  xw            is such that  alphaBest + xw  is the second-best point.
%                It is not defined until  Wset  is true.
%                In some cases,  xw  will replace a previous  xw
%                that has a lower function but has just been excluded
%                from the interval of uncertainty.
%
%  Systems Optimization Laboratory, Stanford University, California.
%  Original version February 1982.  Rev. May 1983.
%  Original f77 version 22-August-1985.
%  14 Sep 1992: Introduced QuitI, QuitF, etc.
%  22 Nov 1995: Altered criterion for reducing the step below tolabs.
%  17 Jul 1997: Removed saved variables for thread-safe version.
%  19 Apr 2000: QuitF only allowed after a move.
%  28 Nov 2019: Switched to the precise Wolfe conditions.
%  01 Dec 2024: Matlab version.
%  09 Feb 2025: Reverse communication removed for UC version.
%
%  This version dated 3 Sept 2026.
%  Philip E. Gill & Jeb Runnoe, University of California, San Diego.
% -------------------------------------------------------------------------
  format compact

  wolfeGtol = wolfeTols(1);
  wolfeFtol = wolfeTols(2);

  maxfevals = jfMax;

% Compute the line search tolerances.

%==========================================================================
% epsaf   is the absolute function precision. If f(x1) and f(x2) are
%         as close as  epsaf,  we cannot safely conclude from the
%         function values alone which of x1 or x2 is a better point.
%
% tolabs  is an estimate of the absolute spacing between points
%         along  p.  This step should produce a perturbation
%         of  epsaf  in the merit function.
%
% tolrel  is an estimate of the relative spacing between points
%         along  p.
%
% toltny  is the minimum allowable absolute spacing between points
%         along  p.
%==========================================================================
  eps0   = eps^(4/5);
  epsrf  = 10^(-9);

  epsaf  = max(epsrf, eps)*(1 + abs(objf0));
  tolax  = eps0;
  tolrx  = eps0;

  xNorm  = norm(x0);
  pnorm  = norm(p);

% Compute the line search tolerances  tolabs, tolrel, toltny, epsaf.

  t      = xNorm*tolrx  +  tolax;
  if  t < pnorm*alfaMax
    tolabs = t/pnorm;
  else
    tolabs = alfaMax;
  end
  tolrel = max(tolrx, eps);

  t      = max(abs(p)./(abs(x0)*tolrx + tolax));

  if  t*tolabs > 1
    toltny = 1 / t;
  else
    toltny = tolabs;
  end

% Compute initial values.

  QuitF     = false;            Improved = false;

  objg0p    = objg0'*p;

  objf      = objf0;            objg     = objg0;
  x         = x0;

  f0        =   0;                              % the shifted objective
  g0        =  (1          - wolfeFtol)*objg0p; % the shifted slope
  targetgL  =  (wolfeGtol  - wolfeFtol)*objg0p;
  targetgU  = -(wolfeGtol  + wolfeFtol)*objg0p;

  numf      =  0;

  alfaBest  =  0;                      fBest  = f0;  gBest = g0;

  BadArgs   = alfaMax <= toltny  ||  g0 >= 0;
  Moved     = false;
  tolA      = tolabs;     tolR   = tolrel;

  Done      = BadArgs;

  if  ~Done
    Bracketed = false;
    Cramped   = alfaMax <= tolA;
    Extrap    = false;
    Wset      = false;
    nSamea    = 0;                     nSameb = 0;

    tol       = tolA;
    tolMax    = tolA  + tolR*alfaMax;
    a         = 0;                     b      = alfaMax + tolMax;
    factor    = 5;
    alfaNew   = alfa0;
    xTry      = alfaNew;
  end

% -------------------------------------------------------------------------
% Main loop.
% a, b and xTry have been shifted by alfaBest.
% -------------------------------------------------------------------------
  while  ~Done

    alfaNew  = alfaBest + xTry;
    if  Bracketed  ||  alfaNew < alfaMax - tolMax

    % The trial step is not near or beyond alfaMax.
    % Check that the function will not be evaluated too close to a or b.
    % It has already been evaluated at both of these points.

      if  xTry <= a + tol  ||  xTry >= b - tol
        if  (a + b)/2 <= 0
          xTry   = - tol;
        else
          xTry   =   tol;
        end
        alfaNew  = alfaBest + xTry;
      end
    else

    % The step is close to, or larger than alfaMax, replace it by
    % alfaMax to force the evaluation of  f  at the boundary.

      Bracketed = true;
      xTry      = alfaMax - alfaBest;
      alfaNew   = alfaMax;
    end

  % -----------------------------------------------------------------------
  % Evaluate the function and gradient at alfaNew.
  % -----------------------------------------------------------------------
    xNew              = x0 + alfaNew*p;
    [objfNew,objgNew] = prob.obj(xNew);
    objgpNew          = objgNew'*p;

    fNew              = objfNew  - objf0 - wolfeFtol*objg0p*alfaNew;
    gNew              = objgpNew         - wolfeFtol*objg0p;

    numf   = numf + 1;  nSamea = nSamea + 1;   nSameb = nSameb + 1;

    if  ~Bracketed
      b    = alfaMax + tolMax - alfaBest;
    end

  % See if the new step is better.  If alfaNew is large enough that
  % fNew can be distinguished numerically from 0,  the function
  % is required to be sufficiently negative.

    Closef     = abs(fNew - fBest) <=  epsaf;
    if  Closef
      Improved = abs(gNew)         <= abs(gBest);
    else
      Improved =     fNew          <      fBest;
    end

    if  Improved
    % ---------------------------------------------------------------------
    % We seem to have an improvement.  The new point becomes the
    % origin and other points are shifted accordingly.
    % ---------------------------------------------------------------------
      fw        = fBest;
      fBest     = fNew;
      gw        = gBest;
      gBest     = gNew;
      alfaBest  = alfaNew;
      Moved     = true;

      a         = a  - xTry;
      b         = b  - xTry;
      xw        = 0  - xTry;
      Wset      = true;
      Extrap    = xw < 0  &&  gBest < 0  ||  xw > 0  &&  gBest > 0;

    % Decrease the length of the interval of uncertainty.

      if  gNew <= 0
        a         = 0;      nSamea = 0;
      else
        b         = 0;      nSameb = 0;
        Bracketed = true;
      end

      objf = objfNew;       objg   = objgNew;

    else
    % ---------------------------------------------------------------------
    % The new function value is no better than the best value so far
    % The origin remains unchanged but the new point may qualify as xw.
    % xTry must be a new bound on the best point.
    % ---------------------------------------------------------------------
      if  xTry <= 0
        a         = xTry;      nSamea = 0;
      else
        b         = xTry;      nSameb = 0;
        Bracketed = true;
      end

    % If xw has not been set or fNew is better than fw, update the
    % points accordingly.

      if  Wset
        Setxw = fNew < fw  ||  ~Extrap;
      else
        Setxw = true;
      end

      if  Setxw
        xw     = xTry;
        fw     = fNew;
        gw     = gNew;
        Wset   = true;
        Extrap = false;
      end
    end

  % -----------------------------------------------------------------------
  % Check the termination criteria.  Wset will always be true.
  % truea  = alfaBest + a;     trueb  = alfaBest + b;
  % -----------------------------------------------------------------------
    tol    = tolA   + tolR*alfaBest;

    Found  = targetgL   <= gBest       &&  gBest <= targetgU;
    QuitF  = numf       >= maxfevals   &&  Moved;
    QuitI  = b - a      <= tol + tol;

    if  QuitI  &&  ~Moved

    % The interval of uncertainty appears to be small enough,
    % but no better point has been found.  Check that changing
    % alfa by b - a changes f by less than epsaf.

      tol   = tol/10;
      tolA  = tol;
      QuitI = tol <= toltny  ||  abs(fw) <= epsaf   &&  gw <= epsaf;
    end

    Done  = QuitF  ||  QuitI  ||  Found;

    if Done, break; end

  % -----------------------------------------------------------------------
  % Proceed with the computation of a trial steplength.
  % The choices are...
  % 1. Parabolic fit using derivatives only, if the f values are
  %    close.
  % 2. Cubic fit for a minimizer, using both f and f'.
  % 3. Damped cubic or parabolic fit if the regular fit appears to
  %    be consistently overestimating the distance to a minimizer.
  % 4. Bisection, geometric bisection, or a step of  tol  if
  %    choices 2 or 3 are unsatisfactory.
  % -----------------------------------------------------------------------
    xmidpt = (a + b)/2;
    s      = 0;               q  = 0;

    if   Closef
    % ---------------------------------------------------------------------
    % Fit a parabola to the two best gradient values.
    % ---------------------------------------------------------------------
      s      = gBest;
      q      = gBest - gw;
    else
    % ---------------------------------------------------------------------
    % Fit cubic through  fBest  and  fw.
    % ---------------------------------------------------------------------
      FitOK  = true;
      rc     = 3*(fBest - fw)/xw + gBest + gw;
      absrc  = abs(rc);
      sc     = sqrt(abs(gBest))*sqrt(abs(gw));

    % Compute  qc =  the square root of  r*r - gBest*gw.
    % The method avoids unnecessary underflow and overflow.

      if    (gw < 0  &&  gBest > 0) || (gw > 0  &&  gBest < 0)
        scale = absrc + sc;
        if  scale == 0
          qc  = 0;
        else
          qc  = scale*sqrt((absrc/scale)^2 + (sc/scale)^2);
        end
      elseif  absrc >= sc
        qc    = sqrt(absrc + sc)*sqrt(absrc - sc);
      else
        FitOK = false;
      end

      if  FitOK         % Compute a minimizer of the fitted cubic.
        if  xw < 0
          qc = - qc;
        end
        s  = gBest -  rc - qc;    q  = gBest - gw - qc - qc;
      end
    end
  % -----------------------------------------------------------------------
  % Construct an artificial interval  (artifa, artifb)  in which
  % the new estimate of a minimizer must lie.  Set a default
  % value of xTry that will be used if the polynomial fit fails.
  % -----------------------------------------------------------------------
    artifa = a;    artifb = b;
    if  ~Bracketed

    % A minimizer has not been bracketed.  Set an artificial
    % upper bound by expanding the interval  xw  by a suitable
    % factor.

      xTry   = - factor*xw;
      artifb =   xTry;
      if  alfaBest + xTry < alfaMax
        factor = 5*factor;
      end
    elseif  Extrap

    % The points are configured for an extrapolation.
    % Set a default value of  xTry  in the interval  (a, b)
    % that will be used if the polynomial fit is rejected.  In
    % the following,  dTry  and  daux  denote the lengths of
    % the intervals  (a, b)  and  (0, xw)  (or  (xw, 0),  if
    % appropriate).  The value of  xTry is the point at which
    % the exponents of  dTry  and  daux  are approximately
    % bisected.

      daux = abs(xw);
      dTry = b - a;
      if  daux >= dTry
        xTry = 5*dTry*(0.1 + dTry/daux)/11;
      else
        xTry = sqrt(daux) * sqrt(dTry)/2;
      end

      if  xw > 0
        xTry = - xTry;
      end

    % Reset the artificial bounds.  If the point computed by
    % extrapolation is rejected,  xTry will remain at the
    % relevant artificial bound.

      if  xTry <= 0
        artifa = xTry;
      else
        artifb = xTry;
      end
    else

    % The points are configured for an interpolation.  The
    % default value xTry bisects the interval of uncertainty.
    % the artificial interval is just (a, b).

      xTry   = xmidpt;
      if  nSamea >= 3  ||  nSameb >= 3

      % If the interpolation appears to be overestimating the
      % distance to a minimizer,  damp the interpolation.

        factor = factor / 5;
        s      = factor * s;
      else
        factor = 1;
      end
    end
  % -----------------------------------------------------------------------
  % The polynomial fits give  (s/q)*xw  as the new step.
  % Reject this step if it lies outside  (artifa, artifb).
  % -----------------------------------------------------------------------
    if  q ~= 0
      if  q < 0
        s = - s;  q = - q;
      end
      if  s*xw >= q*artifa  &&  s*xw <= q*artifb

      % Accept the polynomial fit.

        if  abs(s*xw) >= q*tol
          xTry = (s/q)*xw;
        else
          xTry = 0;
        end
      end
    end
  end % while not Done

% -------------------------------------------------------------------------
% Exit.
% -------------------------------------------------------------------------
  if     BadArgs
    iExit = 8;             % bad arguments
  elseif  Found
    if  alfaBest < alfaMax
      iExit = 1;           % Sufficient decrease
    else
      iExit = 2;           % Suff. Decrease on the boundary
    end
  elseif  Moved
    iExit   = 3;           % Decr at boundary or max funs
  elseif  QuitF
    iExit   = 7;           % No new point after max funs
  elseif  Cramped
    iExit   = 4;           % alfaMax too mall
  else
    iExit   = 6;           % [a,b] too small
  end

% =========================================================================
%  Finish with  x = the best point found so far.
% =========================================================================
  alfa = alfaBest;

  if  Improved
    x  = xNew;
  elseif  alfa > 0
    x  = x0 + alfa*p;
  end
end
