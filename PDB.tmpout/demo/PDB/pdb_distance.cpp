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

#include <vector>
#include <iterator>
#include <fstream>
#include <CGAL/PDB/align.h>
#include <CGAL/PDB/distance.h>


#include <boost/program_options.hpp>

/*!
  \example pdb_distance.cc

  This example shows how to compute various types of distance measures
  between proteins. For cRMS computations, the proteins are first
  aligned and then the cRMS is computed (this could instead be done
  with one function call).
*/

int main(int argc, char *argv[]){
  std::string base_file, input_file, output_file;
  bool print_help=false;
  bool crms=false, drms=false;
  int base_model=0;
  bool warn=false;
  bool all_atoms=false;
  bool no_align=false;
  bool max_cdist=false;
  bool max_ddist=false;
  {
    boost::program_options::options_description o("Allowed options"), po, ao;
    o.add_options()
      ("crms,c", boost::program_options::bool_switch(&crms),
       "Output the cRMS between the two pdbs (after alignment).")
      ("drms,d", boost::program_options::bool_switch(&drms),
       "Output the dRMS between the two pdbs.")
      ("all-atoms,a", boost::program_options::bool_switch(&all_atoms),
       "Output the distances between all atoms, not just the C_alphas.")
      ("no-align,n", boost::program_options::bool_switch(&no_align),
       "Do not rigidly align the proteins before computing cRMS.")
      ("verbose,v", boost::program_options::bool_switch(&warn),
       "Warn about errors parsing pdb file.")
      ("max-c-dist,C", boost::program_options::bool_switch(&max_cdist),
       "Output the max distance between any two corresponding atoms.")
      ("max-d-dist,D", boost::program_options::bool_switch(&max_ddist),
       "Output the max distance difference between pairwise distances.")
      ("help", boost::program_options::bool_switch(&print_help), "Produce help message");

    po.add_options()("input-pdb-0",
		     boost::program_options::value< std::string>(&base_file),
		     "First input file.")
      ("input-pdb-1", boost::program_options::value< std::string>(&input_file),
       "Second input file");
    ao.add(o).add(po);

    boost::program_options::positional_options_description p;
    p.add("input-pdb-0", 1);
    p.add("input-pdb-1", 1);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(ao).positional(p).run(), vm);
    boost::program_options::notify(vm);


    if (base_file.empty() || input_file.empty() || print_help) {
      std::cout << "This program computes the distance between two pdb files..\n";
      std::cout << "usage: " << argv[0] << " input-pdb-0 input-pdb-1\n\n";
      std::cout << o << "\n";
      return EXIT_SUCCESS;
    }
  }


  std::ifstream in(input_file.c_str());

  if (!in){
    std::cerr<< "Error opening input file: " << input_file << std::endl;
    return EXIT_FAILURE;
  }
  std::ifstream bin(base_file.c_str());

  if (!bin){
    std::cerr<< "Error opening input file: " << base_file << std::endl;
    return EXIT_FAILURE;
  }

  CGAL_PDB_NS::PDB input(in, warn);
  CGAL_PDB_NS::PDB base(bin, warn);
  CGAL_PDB_NS::Protein base_protein;
  for (unsigned int i=0; i< base.number_of_models(); ++i){
    const CGAL_PDB_NS::Model &m= base.model(i);
    if (base_model ==0 || m.index()==base_model) {
      for (unsigned int j=0; j< m.number_of_chains(); ++j){
	base_protein= m.chain(j);
      }
    }
  }

  if ( crms && !no_align) {
    for (unsigned int i=0; i< input.number_of_models(); ++i){
      CGAL_PDB_NS::Model &m= input.model(i);
      CGAL_PDB_NS::Protein &p= m.chain(0);
      CGAL_PDB_NS::align_second_protein_to_first(base_protein, p);
    }
  }

  if (crms) {
    std::cout << "cRMS:\n";
    for (unsigned int i=0; i< input.number_of_models(); ++i){
      CGAL_PDB_NS::Model &m= input.model(i);
      CGAL_PDB_NS::Protein &p= m.chain(0);
      double d;
      if (all_atoms) {
	d =CGAL_PDB_NS::no_align_cRMS(base_protein, p);
      } else {
	d =CGAL_PDB_NS::no_align_ca_cRMS(base_protein, p);
      }
      if (d < .00001) d=0;
      std::cout << "Model " << i << " is " << d << std::endl;
    }
  }

  if (max_cdist) {
    std::vector<CGAL_PDB_NS::Point> base_coords;
    if (all_atoms) {
      CGAL_PDB_NS::coordinates(base_protein.atoms_begin(), base_protein.atoms_end(), std::back_inserter(base_coords));
    } else {
      CGAL_PDB_NS::ca_coordinates(base_protein.atoms_begin(), base_protein.atoms_end(), std::back_inserter(base_coords));
    }
    std::cout << "Max dist:\n";
    for (unsigned int i=0; i< input.number_of_models(); ++i){
      CGAL_PDB_NS::Model &m= input.model(i);
      CGAL_PDB_NS::Protein &p= m.chain(0);
      std::vector<CGAL_PDB_NS::Point> coords;
      if (all_atoms) {
	CGAL_PDB_NS::coordinates(p.atoms_begin(), p.atoms_end(), std::back_inserter(coords));
      } else {
	CGAL_PDB_NS::ca_coordinates(p.atoms_begin(), p.atoms_end(), std::back_inserter(coords));
      }
      if (coords.size() != base_coords.size()){
	std::cerr<< "Proteins being compared must have the same number of atoms.\n";
	return EXIT_FAILURE;
      }
      double mmax= 0;
      CGAL_PDB_NS::Squared_distance sd;
      for (unsigned int j=0; j< coords.size(); ++j){
	mmax=std::max BOOST_PREVENT_MACRO_SUBSTITUTION(std::sqrt(sd(coords[j], base_coords[j])), mmax);
      }
      std::cout << "Model " << i << " is " << mmax << std::endl;
    }
  }

  if (drms) {
    std::cout << "\ndRMS:\n";
    for (unsigned int i=0; i< input.number_of_models(); ++i){
      CGAL_PDB_NS::Model &m= input.model(i);
      CGAL_PDB_NS::Protein &p= m.chain(0);
      double d;
      if (all_atoms) {
	d=CGAL_PDB_NS::dRMS(base_protein, p);
      } else {
	d=CGAL_PDB_NS::ca_dRMS(base_protein, p);
      }
      if (d < .00001) d=0;
      std::cout << "Model " << i << " is " << d << std::endl;
    }
  }

  if (max_ddist) {
    CGAL_PDB_NS::Matrix base_matrix;
    if (all_atoms) {
      base_matrix=CGAL_PDB_NS::distance_matrix(base_protein);
    } else {
      base_matrix=CGAL_PDB_NS::ca_distance_matrix(base_protein);
    }
    std::cout << "Max dist:\n";
    for (unsigned int i=0; i< input.number_of_models(); ++i){
      CGAL_PDB_NS::Model &m= input.model(i);
      CGAL_PDB_NS::Protein &p= m.chain(0);
      CGAL_PDB_NS::Matrix input_matrix;
      if (all_atoms) {
	input_matrix= CGAL_PDB_NS::distance_matrix(p);
      } else {
	input_matrix = CGAL_PDB_NS::ca_distance_matrix(p);
      }
      if (base_matrix.dim1() != input_matrix.dim1()){
	std::cerr << "Proteins being compared must have the same number of atoms.\n";
	std::cerr << "Base protein had " << base_matrix.dim1()
		  << " atoms and second protein had " << input_matrix.dim1()
		  << " atoms." << std::endl;
	return EXIT_FAILURE;
      }
      double mmax= 0;
      double maxr=0;
      for (int j=0; j< base_matrix.dim1(); ++j){
	for (int k=0; k< j; ++k){
	  mmax=std::max BOOST_PREVENT_MACRO_SUBSTITUTION(std::abs(base_matrix[j][k]-input_matrix[j][k]), mmax);
	  maxr=std::max BOOST_PREVENT_MACRO_SUBSTITUTION(std::abs(base_matrix[j][k]-input_matrix[j][k])
			/std::max BOOST_PREVENT_MACRO_SUBSTITUTION(base_matrix[j][k],input_matrix[j][k]), maxr);
	}
      }
      std::cout << "Model " << i << " is " << mmax << " with ratio " << maxr << std::endl;
    }
  }

  return EXIT_SUCCESS;

}
