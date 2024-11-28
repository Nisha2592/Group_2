# Set output to a PNG file
set terminal pngcairo size 600,400
set output 'benchmark.png'

# Titles and labels
set title "Benchmark Comparison: Optimized vs Optimized-MPI"
set xlabel "natom"
set ylabel "Time (s)"

# Grid and key
set grid
set key top left

# Plot the data
plot "benchmark_data.dat" using 1:2 with linespoints title "Non-MPI" lc rgb "red" lw 2 pt 7, \
     "benchmark_data.dat" using 1:3 with linespoints title "2-procs" lc rgb "blue" lw 2 pt 7, \
     "benchmark_data.dat" using 1:4 with linespoints title "4-procs" lc rgb "green" lw 2 pt 7, \
     "benchmark_data.dat" using 1:5 with linespoints title "6-procs" lc rgb "black" lw 2 pt 7, \
     ""                   using 1:6 with linespoints title "8-procs" lc rgb "purple" lw 2 pt 7, \
     ""                   using 1:7 with linespoints title "10-procs" lc rgb "yellow" lw 2 pt 7
