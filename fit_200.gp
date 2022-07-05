

set t pdf

set output "fitting_200.pdf"
a=2000
b=0.15 
c=1.3
f(x)=a*exp(-(x-b)**2/(2*c**2))

fit [0.14:0.18] f(x) "ver0.1.0CriticalLam_ThisFileN200Size3000to200Lam0.05to0.25_50_NumNet20000_GTba210_GS875629289_SIRS1489264107025_THRk4_LOCKLimitContacts.dat" u 1:3 via a,b,c

p "ver0.1.0CriticalLam_ThisFileN200Size3000to200Lam0.05to0.25_50_NumNet20000_GTba210_GS875629289_SIRS1489264107025_THRk4_LOCKLimitContacts.dat" u 1:3, f(x)
#0.475972

set output