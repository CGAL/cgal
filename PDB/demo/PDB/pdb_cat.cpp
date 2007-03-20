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
#include <CGAL/PDB/align.h>

#include <boost/program_options.hpp>

int main(int argc, char *argv[]){
  std::vector<std::string> input_files;
  std::string output_file;
  bool print_help=false;
  bool verbose=false;
  bool align=false;

  {
    boost::program_options::options_description o("Allowed options"), po, ao;
    o.add_options()
      ("help", boost::program_options::bool_switch(&print_help), "produce help message")
      ("verbose,v", boost::program_options::bool_switch(&verbose),
       "Print error messages from reading the pdb files.")
      ("align,a", boost::program_options::bool_switch(&align),
       "Align each protein to the previous.");
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
      return EXIT_SUCCESS;
    }

    output_file=input_files.back();
    input_files.pop_back();
  }


  CGAL_PDB_NS::PDB outpdb;
  CGAL_PDB_NS::Model last;
  for (unsigned int i=0; i < input_files.size(); ++i){
    std::ifstream in(input_files[i].c_str());
    if (!in){
      std::cerr << "Error opening input file " << input_files[i] << std::endl;
      return EXIT_FAILURE;
    }
    CGAL_PDB_NS::PDB inpdb(in, verbose);

    for (CGAL_PDB_NS::PDB::Models_iterator mit= inpdb.models_begin(); mit != inpdb.models_end(); ++mit){
      CGAL_PDB_NS::Model m= *mit;
      m.set_index(outpdb.number_of_models());
      if (align && i != 0) {
	for (unsigned int j= 0; j< m.number_of_chains(); ++j) {
	  CGAL_PDB_NS::Protein &p= m.chain(j);
	  const CGAL_PDB_NS::Protein &bp= last.chain(j);
	  CGAL_PDB_NS::align_second_protein_to_first(bp, p);
	}
      }
      outpdb.new_model(m);
      last= m;
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
