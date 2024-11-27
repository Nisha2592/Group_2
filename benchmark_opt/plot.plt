# Set output to a PNG file
set terminal pngcairo size 600,400
set output 'benchmark.png'

# Titles and labels
set title "Benchmark Comparison: Non-Optimized vs Optimized"
set xlabel "natom"
set ylabel "Time (s)"

# Grid and key
set grid
set key top left

# Plot the data
plot "benchmark_data.dat" using 1:2 with linespoints title "Non-Optimized" lc rgb "red" lw 2 pt 7, \
     "benchmark_data.dat" using 1:3 with linespoints title "Optimized" lc rgb "blue" lw 2 pt 7
