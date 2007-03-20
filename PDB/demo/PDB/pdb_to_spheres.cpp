/* Copyright 2004
Stanford University

This file is part of the DSR PDB Library.

The DSR PDB Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

The DSR PDB Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the DSR PDB Library; see the file LICENSE.LGPL.  If not, write to
the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <CGAL/PDB/PDB.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <CGAL/PDB/cgal.h>
#include <iterator>
#include <CGAL/Simple_cartesian.h>

#include <boost/program_options.hpp>

int main(int argc, char *argv[]){
  std::string input_file;
  std::string output_file;
  bool print_help=false;
  bool verbose=false;

  {
    boost::program_options::options_description o("Allowed options"), po, ao;
    o.add_options()
      ("help", boost::program_options::bool_switch(&print_help), "produce help message")
      ("verbose,v", boost::program_options::bool_switch(&verbose),
       "Print error messages from reading the pdb files.");
    po.add_options()
      ("input-pdb", boost::program_options::value< std::string >(&input_file),
       "The input file.")
      ("output-spheres", boost::program_options::value< std::string >(&output_file),
       "The output file.");

    ao.add(o).add(po);

    boost::program_options::positional_options_description p;
    p.add("input-pdb", 1);
    p.add("output-spheres", 1);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(ao).positional(p).run(), vm);
    boost::program_options::notify(vm);


    if (input_file.empty() || output_file.empty() || print_help) {
      std::cout << "Read a PDB file and write spheres:\n# num_spheres \nx y z r\n...\nto the output file.\n";
      std::cout << "usage: " << argv[0] << " input.pdb output.spheres\n";
      std::cout << o << "\n";
      return EXIT_SUCCESS;
    }

  }


  std::ifstream in(input_file.c_str());
  if (!in){
    std::cerr << "Error opening input file " << input_file << std::endl;
    return EXIT_FAILURE;
  }
  CGAL_PDB_NS::PDB inpdb(in, verbose);

  std::ofstream out(output_file.c_str());
  if (!out) {
    std::cerr<< "Error opening output file: " << output_file << std::endl;
    return EXIT_FAILURE;
  }
  typedef CGAL::Simple_cartesian<double> K;

  CGAL_PDB_NS::all_spheres<K>(inpdb.model(0), std::ostream_iterator<K::Sphere_3>(out, "\n"));


  return EXIT_SUCCESS;
}
