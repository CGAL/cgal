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

#include <CGAL/PDB/Transform.h>

#include <boost/program_options.hpp>


int main(int argc, char *argv[]){
  std::string input_file, output_file;
  bool print_help=false;
  bool verbose=false;
  bool dali=false;
  bool matrix=false;

  boost::program_options::options_description o("Allowed options"), po, ao;
  o.add_options()
    ("help", boost::program_options::bool_switch(&print_help),
     "produce help message")
    ("verbose,v", boost::program_options::bool_switch(&verbose),
     "print out any errors that occur during reading of the pdb file.");

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

  if (input_file.empty() || output_file.empty() || print_help || (dali && matrix)) {
    std::cout << "This program produces a pdb identical to the input but with no sidechains.\n";
    std::cout << "usage: " << argv[0] << " input-pdb output-pdb\n\n";
    std::cout << o << "\n";
    return EXIT_SUCCESS;
  }




  // std::cout << input_file << " " << output_template << " " << split_domain << " " << split_chains << std::endl;



  //= new char[strlen(argv[2]+1000)];


  std::ifstream in(input_file.c_str());
  if (!in) {
    std::cerr << "Error opening output file " << input_file << std::endl;
    return EXIT_FAILURE;
  }

  CGAL_PDB_NS::PDB pdb(in, verbose);

  if (verbose) std::cout << "Input PDB has " << pdb.number_of_models() << " models." << std::endl;

  CGAL_PDB_NS::PDB outpdb;

  for (unsigned int i=0;i< pdb.number_of_models(); ++i){
    CGAL_PDB_NS::Model &m= pdb.model(i);
    CGAL_PDB_NS::Model om;

    std::cout << "Model " << i << " has " << m.number_of_chains() << " chains."<< std::endl;
    for (unsigned int j=0; j< m.number_of_chains(); ++j){
      CGAL_PDB_NS::Protein &p= m.chain(j);
      CGAL_PDB_NS::Protein op;
      for (CGAL_PDB_NS::Protein::Residues_iterator rit= p.residues_begin(); rit != p.residues_end(); ++rit){
	CGAL_PDB_NS::Residue nr(rit->type());
	nr.set_index(rit->index());
	if (rit->has_atom(CGAL_PDB_NS::Residue::AL_N)) {
	  nr.set_atom(CGAL_PDB_NS::Residue::AL_N, rit->atom(CGAL_PDB_NS::Residue::AL_N));
	}
	if (rit->has_atom(CGAL_PDB_NS::Residue::AL_CA)) {
	  nr.set_atom(CGAL_PDB_NS::Residue::AL_CA, rit->atom(CGAL_PDB_NS::Residue::AL_CA));
	}
	if (rit->has_atom(CGAL_PDB_NS::Residue::AL_C)) {
	  nr.set_atom(CGAL_PDB_NS::Residue::AL_C, rit->atom(CGAL_PDB_NS::Residue::AL_C));
	}
	op.new_residue(nr);
      }
      om.new_chain(op);
    }
    outpdb.new_model(om);
  }

  std::ofstream out(output_file.c_str());
  if (!out) {
    std::cerr << "Error opening output file " << output_file << std::endl;
    return EXIT_FAILURE;
  }

  outpdb.write(out);

  //delete[] buf;
  return EXIT_SUCCESS;
}
