set t pdf
set output "fitting.pdf"

f(x)=a*exp(-(x-b)**2/(2*c**2))

fit [0.05:0.25] f(x) "ver0.1.0LamScan_C_N50r0.14t0-1_20SamStep5000_Graphba210_GSeed875629289_SS1489264107025_THRk4_LOCKRandom_LT0.1_RT0.05_COMPtrue.dat" u 1:3 via a,b,c 

p "ver0.1.0LamScan_C_N50r0.14t0-1_20SamStep5000_Graphba210_GSeed875629289_SS1489264107025_THRk4_LOCKRandom_LT0.1_RT0.05_COMPtrue.dat" u 1:3, f(x)
#0.475972

set output