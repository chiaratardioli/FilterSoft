UD=42164.1697748545

#####—— Keplerian Elements ——#####

reset

set terminal x11 0 title "Keplerian elements" font "helvetica,10"
unset key
NX=1; NY=6
SX=0.9; SY=0.15; vertMargin=0.75 ; horMargin=20
set bmargin vertMargin
set tmargin vertMargin
set lmargin horMargin
set rmargin 0
## to connect each plot without space
set size 1,1
set origin 0,0
set multiplot
set grid 
set autoscale
###—— Mean anomaly
set size SX,SY
set origin 0,0.05;
set xtic auto offset 0,0.2 font "Helvetica,8" 
set ytic font "Helvetica,8"
set xrange [*:*]
set title ""
set xlabel "Time [Years]" 
set label 1 'Mean anomaly' at graph -0.2, graph 0.6
set label 2 '[rad]' at graph -0.15, graph 0.4
plot 'results_kepl.dat' u($1/2/pi/365):($7) ev 1 w l
###—— Resonant angle
#set size SX,SY
#set origin 0,0.05;
#set xtic auto offset 0,0.2 font "Helvetica,8" 
#set ytic font "Helvetica,8"
#set xrange [*:*]
#set title ""
#set xlabel "Time [Years]" 
#set label 1 'Resonant angle' at graph -0.2, graph 0.6
#set label 2 '[rad]' at graph -0.15, graph 0.4
#plot 'results_kepl.dat' u($1/2/pi/365):($8) ev 1 w l
###—— Argument of perigee
#unset xtic
set origin 0,0.05+SY;
set title ""
set xlabel ""
unset label 1 
set label 3 'Arg. of perigee' at graph -0.2, graph 0.6
plot 'results_kepl.dat' u($1/2/pi/365):($6) ev 1 w l
###—— Longitude of asc. node
set origin 0,0.05+SY*2;
set title ""
set xlabel ""
unset label 3
set label 4 'Long. of asc. node' at graph -0.23, graph 0.6
plot 'results_kepl.dat' u($1/2/pi/365):($5) ev 1 w l, 'results_kepl.dat' u($1/2/pi/365):($9) ev 1 w d , 'results_kepl.dat' u($1/2/pi/365):($10) ev 1 w d
###—— Inclination
set origin 0,0.05+SY*3;
set title ""
set xlabel ""
unset label 4
set label 5 'Inclination' at graph -0.18, graph 0.6
plot 'results_kepl.dat' u($1/2/pi/365):($4) ev 1 w l
###—— Eccentricity
set origin 0,0.05+SY*4;
set title ""
set xlabel ""
unset label 2
unset label 5
set label 6 'Eccentricity' at graph -0.21, graph 0.6
plot 'results_kepl.dat' u($1/2/pi/365):($3) ev 1 w l
###—— Semi-major axis
set origin 0,0.05+SY*5;
set title "Keplerian elements"
set xlabel ""
unset format
unset label 5
unset label 6
set label 7 'Semi-major axis' at graph -0.25, graph 0.6
set label 8 '[km]' at graph -0.1485, graph 0.4
plot 'results_kepl.dat' u($1/2/pi/365):($2*UD) ev 1 w l
unset multiplot

reset

set term postscript eps color dashed defaultplex "Helvetica" 10
set output "kepl.eps"
unset key
NX=1; NY=6
DX=0.075; DY=0.05; SX=0.85; SY=0.25
set bmargin 2; set tmargin 0.5; set lmargin DY; set rmargin DY
## to connect each plot without space
set size SX*NX+0.2,SY*NY+0.15
set origin 0,0
set multiplot
set grid lt 3
set autoscale
###—— Mean anomaly
set size SX,SY
set origin DX+0.075,DY;
set xtic auto
set xrange [*:*]
set title ""
set xlabel "Time [Years]"
set ylabel "Mean anomaly [rad]"
plot 'results_kepl.dat' u($1/2/pi/365):($7) ev 1 w l
###—— Resonant angle
#set size SX,SY
#set origin DX+0.075,DY;
#set xtic auto
#set xrange [*:*]
#set title ""
#set xlabel "Time [Years]"
#set ylabel "Resonant angle [rad]"
#plot 'results_kepl.dat' u($1/2/pi/365):($8) ev 1 w l lw 3
###—— Argument of perigee
#unset xtic
set origin DX+0.075,DY+SY#+0.0075;
set title ""
set xlabel ""
set ylabel "Argument of perigee [rad]"
plot 'results_kepl.dat' u($1/2/pi/365):($6) ev 1 w l lw 3
###—— Longitude of asc. node
set origin DX+0.075,DY+SY*2#+0.015;
set title ""
set xlabel ""
set ylabel "Longitude of asc. node [rad]"
plot 'results_kepl.dat' u($1/2/pi/365):($5) ev 1 w l lw 3, 'results_kepl.dat' u($1/2/pi/365):($9) ev 1 w d, 'results_kepl.dat' u($1/2/pi/365):($10) ev 1 w d
###—— Inclination
set origin DX+0.075,DY+SY*3#+0.0225;
set title ""
set xlabel ""
set ylabel "Inclination [rad]"
plot 'results_kepl.dat' u($1/2/pi/365):($4) ev 1 w l
###—— Eccentricity
set origin DX+0.075,DY+SY*4#+0.03;
set title ""
set xlabel ""
set ylabel "Eccentricity"
plot 'results_kepl.dat' u($1/2/pi/365):($3) ev 1 w l
###—— Semi-major axis
set origin DX+0.075,DY+SY*5#+0.0375;
set title "Keplerian elements"
set xlabel ""
set ylabel "Semi-major axis [km]"
plot 'results_kepl.dat' u($1/2/pi/365):($2*UD) ev 1 w l
unset multiplot

#####—— Energy ——#####

reset

#set term postscript eps color dashed defaultplex "Helvetica" 10
#set output "energy.eps"
#unset key
#set grid
#set lmargin 8
#set bmargin 2
#set size 1,1
#set origin 0,0
#set xtic auto 
#set ytic auto 
#set title "Relative error in energy"
#set xlabel "Time [years]" 
#set ylabel "Relative error in energy" 
#plot 'results_cart.dat' u($1/2/pi/365):($10) ev 1 w l 

#reset 

#set terminal x11 2 title "Energy" font "helvetica,10"
#unset key
#set grid
#set lmargin 8
#set bmargin 2
#set size 1,1
#set origin 0,0
#set xtic auto 
#set ytic auto 
#set title "Relative error in energy"
#set xlabel "Time [years]" 
#set ylabel "Relative error in energy" 
#plot 'results_cart.dat' u($1/2/pi/365):($10) ev 1 w l 

#####—— Cartesian coordinates ——#####

reset
vertMargin=0.75 
set bmargin 0.75
set tmargin 0.75
set lmargin 1
set rmargin 1
set terminal x11 2 title "Cartesian coordinates" font "helvetica,10"
unset key
SX=0.95
SY=0.225
vertMargin=0.75 ; horMargin=20
set bmargin vertMargin
set tmargin vertMargin
set lmargin horMargin
set rmargin 0
set size 1,1
set origin 0,0
set multiplot
set grid 
set autoscale
###—— Energy
set size SX,SY
set origin 0,0.05;
set xtic auto offset 0,0.2 font "Helvetica,8" 
set ytic font "Helvetica,8"
#set logscale y
set title ""
set xlabel "Time [Years]" 
set ylabel "Relative error in energy"
plot 'results_cart.dat' u($1/2/pi/365):($10) ev 1 w l
#unset logscale y
###—— Z
#unset xtic
set origin 0,0.05+SY;
set title ""
set xlabel ""
set ylabel "Z [km]"
plot 'results_cart.dat' u($1/2/pi/365):($8*UD) ev 1 w l
###—— Y
set origin 0,0.05+SY*2;
set title ""
set xlabel ""
set ylabel "Y [km]"
plot 'results_cart.dat' u($1/2/pi/365):($7*UD) ev 1 w l
###—— X
set origin 0,0.05+SY*3;
set title "Cartesian coordinates and relative error in energy"
set xlabel ""
set ylabel "X [km]"
plot 'results_cart.dat' u($1/2/pi/365):($6*UD) ev 1 w l
unset multiplot

reset

set term postscript eps color dashed defaultplex "Helvetica" 10
set output "cart.eps"
unset key
SX=0.95; SY=0.225
set bmargin 2
set tmargin 0.5
set lmargin 20
set rmargin 0
set size SX+0.1,SY*4+0.15
set origin 0,0
set multiplot
set grid lt 3
set autoscale	
###—— Energy
set size SX,SY
set origin 0,0.05
set xtic auto
#set logscale y
set title ""
set xlabel "Time [Years]"
set ylabel "Relative error in energy"
plot 'results_cart.dat' u($1/2/pi/365):($10) ev 1 w l
#unset logscale y
###—— Z
#unset xtic
set origin 0,0.05+SY
set title ""
set xlabel ""
set ylabel "Z [km]"
plot 'results_cart.dat' u($1/2/pi/365):($8*UD) ev 1 w l
###—— Y
set origin 0,0.05+SY*2
set title ""
set xlabel ""
set ylabel "Y [km]"
plot 'results_cart.dat' u($1/2/pi/365):($7*UD) ev 1 w l
###—— X
set origin 0,0.05+SY*3
set title "Cartesian coordinates and relative error in energy"
set xlabel ""
set ylabel "X [km]"
plot 'results_cart.dat' u($1/2/pi/365):($6*UD) ev 1 w l
unset multiplot

#####—— Resonant angle ——#####

reset

set term postscript eps color dashed defaultplex "Helvetica" 10
set output "resonant.eps"
unset key
set grid
set lmargin 8
set bmargin 2
set size 1,1
set origin 0,0
set xtic auto 
set ytic auto 
set title "Resonant angle"
set xlabel "Time [years]" 
set ylabel "theta [rad]"
plot 'results_kepl.dat' u($1/2/pi/365):($8) ev 1 w l

reset 

set terminal x11 10 title "Resonant angle" font "helvetica,10"
unset key
set grid
set lmargin 8
set bmargin 2
set size 1,1
set origin 0,0
set xtic auto 
set ytic auto 
set title "Resonant angle"
set xlabel "Time [years]" 
set ylabel "theta [rad]" 
plot 'results_kepl.dat' u($1/2/pi/365):($8) ev 1 w l

## #####—— Cartesian coordinates (3D) ——#####

## reset 

## set term postscript eps color dashed defaultplex "Helvetica" 10
## set output "cart3D.eps"
## unset key
## set grid
## set autoscale
## set lmargin 5
## set size 1,1
## set origin 0,0
## set xtic auto font "Helvetica,8" offset -0.5,-0.5
## set ytic auto font "Helvetica,8" offset 0.5,-0.5
## set ztic auto font "Helvetica,8"
## set title "Orbit"
## set xlabel "x [km]" 2,0
## set ylabel "y [km]" 0,0
## set zlabel "z [km]" 0,0
## splot 'results_cart.dat' u($6*UD):($7*UD):($8*UD) ev 1 w p 

## reset 

## set terminal x11 1 title "Orbit" font "helvetica,10"
## unset key
## set grid
## set autoscale
## set lmargin 5
## set size 1,1
## set origin 0,0
## set xtic auto font "Helvetica,8" offset -0.5,-0.5
## set ytic auto font "Helvetica,8" offset 0.5,-0.5
## set ztic auto font "Helvetica,8"
## set title "Orbit"
## set xlabel "x [km]"
## set ylabel "y [km]"
## set zlabel "z [km]"
## splot 'results_cart.dat' u($6*UD):($7*UD):($8*UD) ev 1 w p 






