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
#include <cctype>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <sstream>

/*!
  \example pdb_cat.cpp

  This example shows how to concatenate several PDB files into one.
*/



int main(int argc, char *argv[]){
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

    
  if (input_files.size() <2 || print_help) {
    std::cout << "Concatenate a bunch of pdb files into one pdb file with many models.\n";
    std::cout << "usage: " << argv[0] << " input-pdb-0 input-pdb-1 ... output-pdb\n";
    std::cout << o << "\n";
    return EXIT_FAILURE;
  }

  output_file=input_files.back();
  input_files.pop_back();

  
  CGAL::PDB::PDB outpdb;
  char next_chain='A';
  outpdb.insert(CGAL::PDB::PDB::Model_key(0), CGAL::PDB::Model());
  CGAL::PDB::Model& outm= outpdb.find(CGAL::PDB::PDB::Model_key(0))->model();
  for (unsigned int i=0; i < input_files.size(); ++i){
    std::ifstream in(input_files[i].c_str());
    if (!in){
      std::cerr << "Error opening input file " << input_files[i] << std::endl;
      return EXIT_FAILURE;
    }
    CGAL::PDB::PDB inpdb(in, verbose);
    
    for (CGAL::PDB::PDB::Model_iterator mit= inpdb.models_begin(); 
         mit != inpdb.models_end(); ++mit){
      const CGAL::PDB::Model &m= mit->model();
      for (CGAL::PDB::Model::Chain_const_iterator cit= m.chains_begin();
           cit != m.chains_end(); ++cit) {
        outm.insert(CGAL::PDB::Model::Chain_key(next_chain), cit->chain());
        ++next_chain;
      }
      for (CGAL::PDB::Model::Heterogen_const_iterator cit= m.heterogens_begin(); 
           cit != m.heterogens_end(); ++cit) {
        outm.insert(cit->key(), cit->heterogen());
      }
      
    }
  }
  
  if (output_file.empty()){
    outpdb.write(std::cout);
  } else {
    std::ofstream out(output_file.c_str());
    outpdb.write(out);
  }


  return EXIT_SUCCESS;
}
