#set term png
#set output "configuration.png"
#set terminal png size 2000,2000

set term gif size 700,700 animate delay 10
set output "P5-configuration.gif"


set tmargin at screen 0.9
set bmargin at screen 0.07
set rmargin at screen 0.95
set lmargin at screen 0.07

set key  top left font "Times, 10"
set key spacing 2

set xtics font "Times, 10"  offset 0,0
set ytics font "Times, 10"  offset 0,0

set title "Configuaració Spins" font "Times, 20" offset 0,0

file = "P5-configuration.conf"
L = 32

set xrange [0.5:L+0.5]
set yrange [0.5:L+0.5]

do for [j = 0:4000:10] {
set label 2 sprintf('IMC: %5i',j) at 2,34 left front font 'Times,10'
set label 1 sprintf('T = 1.0') at 28,34 left front font 'Times,10'
plot file index j using 1:2 with points pt 5 ps 2 lc rgb "#990000" t ""
}