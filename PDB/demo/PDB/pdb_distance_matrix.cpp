/* Copyright 2004
   Stanford University

   This file is part of the CGAL PDB Library.

   The CGAL PDB Library is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

   The DSR PDB Library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
   MA 02111-1307, USA. */

#include <CGAL/PDB/PDB.h>
#include <CGAL/PDB/distance.h>
#include <CGAL/PDB/range.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include <Magick++.h>


/*!
  \example pdb_distance_matrix.cc

  This example shows how to write a distance matrix to a file.
*/

int main(int argc, char *argv[]){
  using namespace CGAL::PDB;
  std::string input_file, output_file;
  bool print_help=false;
  bool verbose=false;
  bool all_atoms=false;
  double contact_map_threshold=std::numeric_limits<double>::infinity();

  {
    boost::program_options::options_description o("Allowed options"), po, ao;
    o.add_options()
      ("help", boost::program_options::bool_switch(&print_help), "produce help message")
      ("verbose,v", boost::program_options::bool_switch(&verbose), "print verbose error messages")
      ("all-atoms,a", boost::program_options::bool_switch(&all_atoms),
       "Output the distances between all atoms, not just the C_alphas.")
      ("image,i",boost::program_options::value<std::string>(&output_file), "Output image file.")
      ("contact-threshold,c",boost::program_options::value<double>(&contact_map_threshold), "Output a contact map with the given threshold. (assuming the image argument is passed).");
    po.add_options()
      ("input-pdb", boost::program_options::value< std::string>(&input_file), "input file");

    ao.add(o).add(po);
    
    boost::program_options::positional_options_description p;
    p.add("input-pdb", 1);
    p.add("output-image", 1);
    
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(ao).positional(p).run(), vm);
    boost::program_options::notify(vm);
    
    
    
    //boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    //boost::program_options::notify(vm);    
    
    if (input_file.empty() || print_help) {
      std::cout << "This program writes the distance matrix to an image file.\n";
      std::cout << "usage: " << argv[0] << " input-pdb\n\n";
      std::cout << o << "\n";
      return EXIT_FAILURE;
    }
  }
  
  
  std::ifstream in(input_file.c_str());
  if (!in) {
    std::cerr << "Error opening input file " << input_file << std::endl;
    return EXIT_FAILURE;
  }

  PDB pdb(in, verbose);

  if (verbose) std::cout << "Input PDB has " << CGAL::PDB::distance(pdb.models()) << " models." << std::endl;
  if (CGAL::PDB::distance(pdb.models()) != 1){
    std::cout << "Attempting to write multiple image files: assuming the output argument has a %d in it.\n";
  }

  CGAL_PDB_FOREACH(PDB::Model_pair& m, pdb.models()) {
    Matrix arr;
    if (all_atoms) {
      arr= distance_matrix(make_point_range(make_atom_range(m.model().atoms())));
    } else {
      arr= distance_matrix(make_point_range(make_atom_range(make_ca_range(m.model().atoms()))));
    }

    if (!output_file.empty()) {

      char buf[100000];
      if (CGAL::PDB::distance(pdb.models()) != 1){
	sprintf(buf, output_file.c_str(), m.key().index());
      } else {
	sprintf(buf, output_file.c_str());
      }
      Magick::Geometry geom(arr.dim1(),arr.dim2());


      Magick::Image im(geom, "red");

      if (contact_map_threshold != std::numeric_limits<double>::infinity()){
	for (int j=0; j< arr.dim1(); ++j){
	  for (int k=0; k< arr.dim2(); ++k){
	    if (arr[j][k] < contact_map_threshold){
	      im.pixelColor(j, k, Magick::ColorGray(1));
	    } else {
	      im.pixelColor(j,k, Magick::ColorGray(0));
	    }
	  }
	}
      } else {
	double max=0;
	for (int i=0; i< arr.dim1(); ++i){
	  for (int j=0; j< arr.dim2(); ++j){
	    max=std::max(max, arr[i][j]);
	    //min=std::min(min, arr[i][j]);
	  }
	}
	//std::cout << "Maximum distance is " << max << std::endl;
	for (int i=0; i< arr.dim1(); ++i){
	  for (int j=0; j< arr.dim2(); ++j){
	    im.pixelColor(i,j, Magick::ColorGray(arr[i][j]/max));
	  }
	}
      }
	
      im.write(buf);
    }
  }
  return EXIT_SUCCESS;
  }

