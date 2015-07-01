To compile all of the benchmarks in the user manual, simply run:

> ./run_benchmarks.sh

This will generate 3 sets of files.

Data Files:

- Files are of the form: "_modeldata_<modelname>.dat" , where
  - <modelname> is the base name of the model used
- These simply contain the raw data points from the benchmarks, you do not need to interact with them directly

Table Files:

- Files are of the form: "benchmark_table_##.txt" , where
  - ## is the number of source points used for that trial
- These files contain text for tables that can be pasted directly into the benchmarks section of the user manual

Comparison Figures

- Files are of the form: "benchmark_plot_<modelname>_<runtype>".png , where
  - <modelname> is the base name of the model used
  - <runtype> is one of 'query', 'construction', or 'memory'
- These are images with scatter-plots tracking a specific variable for each model over multiple different number of source points
- These figures should be pasted into the doc/Surface_mesh_shortest_path/fig/ subdirectory

If you want to change the set of models used, simply edit the file "testModels.txt" to add or remove the models you wish to test on (make sure that the model files exist in the data/ directory).
