using Plots
using DelimitedFiles

samples = readdlm("/home/antonio/Documentos/points.txt", ' ', Float64, '\n')

interior = scatter(samples[:, 1], samples[:, 2], samples[:, 3], markeralpha = 0.1, ms=1.1, color=RGB(1, 0, 0), title="Bench", label="samples")

wp_bench = readdlm("/home/antonio/Documentos/gsoc/cgal-public-dev/Barycentric_coordinates_3/benchmark/Barycentric_coordinates_3/build/wp_bench.txt", ' ', Float64, '\n')
dh_bench = readdlm("/home/antonio/Documentos/gsoc/cgal-public-dev/Barycentric_coordinates_3/benchmark/Barycentric_coordinates_3/build/dh_bench.txt", ' ', Float64, '\n')
mv_bench = readdlm("/home/antonio/Documentos/gsoc/cgal-public-dev/Barycentric_coordinates_3/benchmark/Barycentric_coordinates_3/build/mv_bench.txt", ' ', Float64, '\n')
num_bench = readdlm("/home/antonio/Documentos/gsoc/cgal-public-dev/Barycentric_coordinates_3/benchmark/Barycentric_coordinates_3/build/num_bench.txt", ' ', Float64, '\n')

plot(num_bench, [wp_bench, dh_bench, mv_bench], label=["WP" "DH" "MV"], xaxis=("n", :log), yaxis=("time (sec)", :log), title="Log-log scale plot")
