set t pdf
set output "fitting_crit_lam_for_limiting.pdf"


f(x)= l+ a*(x**b)

fit f(x) "critical_lambda_for_limit.dat" u 1:2 via a,b,l

p "critical_lambda_for_limit.dat" u 1:2, f(x)
#0.475972

set output