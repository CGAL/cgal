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
#include <CGAL/PDB/range.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cctype>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <sstream>

/*!
  \example pdb_cat.cpp

  This example shows how to print info about a pdb.
*/



int main(int argc, char *argv[]){
  using namespace CGAL::PDB;
  std::string output_file;
  std::vector<std::string> input_files;
  bool print_help=false;
  bool verbose=false;

  boost::program_options::options_description o("Allowed options"), po, ao;
  o.add_options()
    ("help", boost::program_options::bool_switch(&print_help),
     "produce help message")
    ("verbose,v", boost::program_options::bool_switch(&verbose),
     "print out verbose messages about reading and writing pdb files");
  po.add_options()
    ("input-pdbs", boost::program_options::value< std::vector<std::string> >(&input_files),
     "The input and output files.");

  ao.add(o).add(po);

  boost::program_options::positional_options_description p;
  p.add("input-pdbs", -1);

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
                                options(ao).positional(p).run(), vm);
  boost::program_options::notify(vm);

    
  if (input_files.size() <1 || print_help) {
    std::cout << "Print info about a pdb file.\n";
    std::cout << "usage: " << argv[0] << " input-pdb-0 input-pdb-1 ... \n";
    std::cout << o << "\n";
    return EXIT_FAILURE;
  }

  for (unsigned int i=0; i < input_files.size(); ++i){
    std::ifstream in(input_files[i].c_str());
    if (!in){
      std::cerr << "Error opening input file " << input_files[i] << std::endl;
      continue;
    }
    PDB inpdb(in, verbose);
    std::cout << "File " << input_files[i] 
              << " " << CGAL::PDB::distance(inpdb.models()) << " models" << std::endl;
    
    CGAL_PDB_FOREACH(const PDB::Model_pair& m, inpdb.models()) {
      std::cout << "Model " << m.key() << std::endl;
      CGAL_PDB_FOREACH(const Model::Chain_pair &c,
                       m.model().chains()) {
        std::cout << " Chain " << c.key() << " has "
                  << CGAL::PDB::distance(c.chain().monomers())
                  << " residues" << std::endl;
      }
    }
  }


  return EXIT_SUCCESS;
}
