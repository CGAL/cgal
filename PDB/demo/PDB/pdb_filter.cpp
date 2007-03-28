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

#include <boost/program_options.hpp>



int main(int argc, char *argv[]){
  std::string input_file, output_file, nearby_file;

  bool print_help=false;
  bool verbose=false;
  double dist_threshold=3;

  boost::program_options::options_description o("Allowed options"), po, ao;
  o.add_options()
    ("help", boost::program_options::bool_switch(&print_help),
     "produce help message")
    ("verbose,v", boost::program_options::bool_switch(&verbose),
     "print out verbose messages about reading and writing pdb files")
    ("distance,d", boost::program_options::value<double>(&dist_threshold),
     "The maximum distance to find neighbors when looking for nearby atoms.")
    ("nearby-file,n", boost::program_options::value<std::string>(&nearby_file),
     "Only output atoms which are near atoms of a similar type in this file.");
  po.add_options()
    ("input-pdb", boost::program_options::value< std::string>(&input_file),
     "input file")
    ("output-pdb", boost::program_options::value< std::string>(&output_file),
     "Output file.");

  ao.add(o).add(po);

  boost::program_options::positional_options_description p;
  p.add("input-pdb", 1);
  p.add("output-pdb", 2);

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				options(ao).positional(p).run(), vm);
  boost::program_options::notify(vm);



  //boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
  //boost::program_options::notify(vm);

  if (input_file.empty()
      || output_file.empty()
      || print_help) {
    std::cout << "This program filters a pdb file, removing some of the atoms and residues.\n";
    std::cout << "useage: " << argv[0]
	      << " [-c] input-pdb output.pdb\n" << std::endl;
    std::cout << o << "\n";
    return EXIT_SUCCESS;
  }

  // std::cout << input_file << " " << output_template << " " << split_domain << " " << split_chains << std::endl;

  std::ifstream in(input_file.c_str());
  if (!in) {
    std::cerr << "Error opening input file " << input_file << std::endl;
    return EXIT_FAILURE;
  }

  //= new char[strlen(argv[2]+1000)];
  CGAL_PDB_NS::PDB pdb(in, verbose);

  std::cout << "Input PDB has " << pdb.number_of_models() << " models." << std::endl;



  if (!nearby_file.empty()){
    std::ifstream ns(input_file.c_str());
    if (!ns) {
      std::cerr << "Error opening input file " << nearby_file << std::endl;
      return EXIT_FAILURE;
    }
    CGAL_PDB_NS::PDB fpdb(ns, verbose);

    CGAL_PDB_NS::Squared_distance sd;

    for (unsigned int i=0; i< pdb.number_of_models(); ++i){
      CGAL_PDB_NS::Model &m= pdb.model(i);
      CGAL_PDB_NS::Model &fm= fpdb.number_of_models()==1?fpdb.model(0):fpdb.model(i);

      for (unsigned int j=0; j< m.number_of_chains(); ++j){
	CGAL_PDB_NS::Protein &c= m.chain(j);
	std::vector<CGAL_PDB_NS::Atom::Index> to_erase;


	for (CGAL_PDB_NS::Protein::Const_atoms_iterator it= c.atoms_begin(); it != c.atoms_end(); ++it){
	  bool found=false;
	  for (unsigned int j=0; j< fm.number_of_chains(); ++j){
	    const CGAL_PDB_NS::Protein &fc= fm.chain(j);
	    for (CGAL_PDB_NS::Protein::Const_atoms_iterator fit= fc.atoms_begin(); fit != fc.atoms_end(); ++fit){
	      if (fit->second.type()== it->second.type()) {
		double dd= sd(fit->second.cartesian_coords(), it->second.cartesian_coords());
		if (dd < dist_threshold*dist_threshold){
		  found=true;
		  break;
		}
	      }
	    }
	    if (found==true) break;
	  }
	  if (!found) {
	    to_erase.push_back(it->second.index());
	  }
	}

	for (unsigned int i=0; i< to_erase.size(); ++i){
	  c.erase_atom(to_erase[i]);
	}

      }
    }
  }




  std::ofstream out(output_file.c_str());
  if (!out) {
    std::cerr << "Error opening output file " << output_file << std::endl;
    return EXIT_FAILURE;
  }
  pdb.write(out);
  //delete[] buf;
  return EXIT_SUCCESS;
}
