function [prob] = getCUSTomProblem()
   prob      = struct();
   prob.n    = 2;
   prob.m    = 0;
   prob.x    = [ -1.2 1.0 ]';
   prob.name = "Rosenbrock";
   prob.obj  = @(x)custom_obj(x);
end

function [f, g] = custom_obj(x)
  x1 = x(1);                    x2 = x(2);
  a  = x2-x1^2;                 b  = 1-x1;
  f  = 100*a^2 + b^2;       
  g  = [200*a*(-2*x1)-2*b
	200*a            ];
end
