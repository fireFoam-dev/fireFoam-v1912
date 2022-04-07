#Compute mean flame height, starting from t=20 s
f(x) = h
tstartfit=20.
fit [ x=tstartfit:* ] f(x) 'outFlameHeight' via h
theTitle=sprintf("mean = %.2f [m]",h)
set xlabel "Time [s]"
set ylabel "Flame height [m]"
p 'outFlameHeight' u 1:2 w l lw 2 sm bez notitle, f(x) w l lw 2 lc rgb 'red' t theTitle
