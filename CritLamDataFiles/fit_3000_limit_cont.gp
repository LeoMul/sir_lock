set t pdf
set output "fitting_3000_limit_cont.pdf"


f(x)=a*exp(-(x-(0.05 + (0.05 - 0.25)/(1 + b^2)))**2/(2*c**2))
n = 3000*3000


fit [0.05:0.2] f(x) "ver0.1.0CriticalLam_ThisFileN3000Size3000to200Lam0.05to0.25_50_NumNet20000_GTba210_GS875629289_SIRS1489264107025_THRk4_LOCKLimitContacts.dat" u ($1):($3/n) via a,b,c 

p "ver0.1.0CriticalLam_ThisFileN3000Size3000to200Lam0.05to0.25_50_NumNet20000_GTba210_GS875629289_SIRS1489264107025_THRk4_LOCKLimitContacts.dat" u ($1):($3/n), f(x)
#0.475972

set output