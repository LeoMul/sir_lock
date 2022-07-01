set t pdf
set output "fitting_2800.pdf"

f(x)=a*exp(-(x-b)**2/(2*c**2))

fit [0.1:0.2] f(x) "ver0.1.0CriticalLam_ThisFileN2800Size200to3200Lam0.05to0.25_20_NumNet100000_GTsw0.1_GS875629289_SIRS1489264107025_THRk4.dat" u 1:3 via a,b,c 

p "ver0.1.0CriticalLam_ThisFileN2800Size200to3200Lam0.05to0.25_20_NumNet100000_GTsw0.1_GS875629289_SIRS1489264107025_THRk4.dat" u 1:3, f(x)
#0.475972

set output