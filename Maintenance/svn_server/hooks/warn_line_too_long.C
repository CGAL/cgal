// Takes filenames as argument, and prints a warning for each too long line
// (>= 80 characters).
// First argument is the name of the directory.
// Following arguments are file names in this directory.
// Directory must not contain "/test".
// File names must end with ".h", ".hpp", ".C", ".c", ".cpp".
//
// Sylvain Pion. 2001, 2003, 2004.
// $Date$
// $Revision$

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main(int argc, char **argv)
{
  if (argc < 3) {
    std::cerr << "Input file name(s) as argument" << std::endl;
    return -1;
  }

  // Check we are not in test directory.
  std::string dirname (argv[1]);
  if (dirname.find("/test") != std::string::npos)
    return 0;

  // Iterate over all files.
  for (int file = 2; file < argc; ++file) {
    std::string filename(argv[file]);

    // Test if filename ends with ".h", ".hpp", ".C", ".c", ".cpp".
    if (filename.find(".h") == std::string::npos &&
        filename.find(".hpp") == std::string::npos &&
        filename.find(".c") == std::string::npos &&
        filename.find(".C") == std::string::npos &&
        filename.find(".cpp") == std::string::npos)
	continue;

    std::ifstream iFile(filename.c_str());

    if (!iFile.is_open())
      std::cerr << "could not open file " << filename << std::endl;

    unsigned line_number=0;
    std::string line;

    while (getline(iFile, line)) {
      ++line_number;
      // We don't warn if line_number <= 20 because the headers usually
      // emit the warning validly.
      if (line.size() >= 80 && line_number > 20) {
        std::cerr << std::endl << std::endl << std::endl
                  << "WARNING WARNING WARNING : Line " << line_number
	          << " of file " << filename
	          << " is too long (" << line.size() << " characters)"
		  << std::endl << std::endl << std::endl;
      }
    }
  }

  return 0;
}
