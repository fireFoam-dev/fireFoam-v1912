at(file, row, col) = system( sprintf("awk -v row=%d -v col=%d 'NR == row {print $col}' %s", row, col, file) )
hRef = 0.0

#Enthalpy conservation
allPlanes = "0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8  0.9     1.0 1.1 1.2     1.3     1.4     1.5     1.6     1.7     1.8     1.9 2.0     2.1     2.2     2.3     2.4     2.5     2.6     2.7     2.8     2.9     3.0"
fname = "enthalpy_vs_z.xy"
set print fname
#first compute averages at each plane: h hs hc (in columns 7 9 8 respectively)
f(x) = hmean
g(x) = hsmean
h(x) = hcmean
tstart = 0
thisFile(x) = sprintf("postProcessing/plane_z%s/%d/surfaceFieldValue.dat",word(allPlanes,x),tstart)
print "#height [m]    H_mean    Hs_mean    Hc_mean    (Hc_mean(0) - Hc_mean)"
i=1
hcmean0 = -at(thisFile(i),8,8)
hmean = at(thisFile(i),8,7)
hsmean = at(thisFile(i),8,9)
hcmean = at(thisFile(i),8,8)
print word(allPlanes,i)," ", -hmean," ", hsmean," ", -hcmean," ", 0.0
do for [i=2:words(allPlanes)] {
    hmean = at(thisFile(i),8,7)
    hsmean = at(thisFile(i),8,9)
    hcmean = at(thisFile(i),8,8)

    print word(allPlanes,i)," ", hmean," ", hsmean," ", hcmean," ", (hcmean0 - hcmean)
}
set print "-"

set grid
set key on left bottom opaque notitle
p 'enthalpy_vs_z.xy' u 1:($5*1e-3) w l lw 2 t 'Hc(0)-Hc', 'enthalpy_vs_z.xy' u 1:($3*1e-3) w l lw 2 t 'Hs', 'enthalpy_vs_z.xy' u 1:($4*1e-3) w l lw 2 t 'Hc', 'enthalpy_vs_z.xy' u 1:($2*1e-3) w l lw 2 t 'H'
