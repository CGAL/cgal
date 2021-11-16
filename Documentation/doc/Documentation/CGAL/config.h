namespace CGAL{
/// returns the full path of a data file from \cgal.
/// The data files are located in the `data` directory of a \cgal release
/// or in the directory `Data/data` from a git branch checkout.
/// That function uses as prefix the environment variable `CGAL_DATA_DIR` if set, and
/// the value of the macro `CGAL_DATA_DIR` otherwise.
/// When using `cmake` and a standard \cgal setup, linking the target using that function
/// with the target `CGAL::Data` will automatically set the macro `CGAL_DATA_DIR`
///  pointing to the data directory of the \cgal version used.
/// The function will attempt to open the file at the returned location
/// and will print a warning message via `std::cerr` if the file could not be opened.
std::string data_file_path(const std::string& filename);
}
