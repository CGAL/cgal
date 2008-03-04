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
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <fstream>

int main(int argc, char *argv[]) {
  std::string input_file, output_file;
  bool print_help=false;
  bool verbose=false;

  CGAL_SET_LOG_LEVEL(CGAL::Log::LOTS);

  boost::program_options::options_description o("Allowed options"), po, ao;
  o.add_options()
    ("help", boost::program_options::bool_switch(&print_help),
     "produce help message")
    ("verbose,v", boost::program_options::bool_switch(&verbose),
     "print out verbose messages about reading and writing pdb files");
  po.add_options()
    ("input-pdb", boost::program_options::value< std::string>(&input_file),
     "input file")
    ("output-pdb", boost::program_options::value< std::string>(&output_file),
     "output file");

  ao.add(o).add(po);

  boost::program_options::positional_options_description p;
  p.add("input-pdb", 1);
  p.add("output-pdb", 1);

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
                                options(ao).positional(p).run(), vm);
  boost::program_options::notify(vm);

  std::ifstream fi(input_file.c_str());
  CGAL::PDB::PDB in(fi);
  std::ofstream fo(output_file.c_str());
  in.write(fo);
  return EXIT_SUCCESS;
}
