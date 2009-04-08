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
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Simple_cartesian.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cctype>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <sstream>

/*!
  \example pdb_transform.cpp

  This example shows how to tranform all the coordinates of a PDB.
*/


int main(int argc, char *argv[]){
  using namespace CGAL::PDB;
  std::string input_file, output_file;
  bool print_help=false;
  bool verbose=false;

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

    
  if (input_file.empty() || output_file.empty() || print_help) {
    std::cout << "Transform all the cordinates of a pdb. Transform is read from stdin\n";
    std::cout << "usage: " << argv[0] << " input output\n";
    std::cout << o << "\n";
    return EXIT_FAILURE;
  }

  
  
  double c[3][4];
  for (unsigned int i=0; i< 3; ++i) {
    char buf[1000];
    std::cin.getline(buf,1000);
    if (!std::cin) {
      std::cerr << "error reading matrix from standard input.\n";
      return EXIT_FAILURE;
    }
    std::istringstream iss(buf);
    for (unsigned int j=0; j< 4; ++j) {
      iss >> c[i][j];
    }
    if (!iss) {
      std::cerr << "error reading matrix from standard input.\n";
      return EXIT_FAILURE;
    }
  }
  CGAL::Aff_transformation_3<Kernel> 
    t(c[0][0], c[0][1], c[0][2], c[0][3],
      c[1][0], c[1][1], c[1][2], c[1][3],
      c[2][0], c[2][1], c[2][2], c[2][3]);
  Kernel::Point_3 pt(0,0,1);
  // std::cout << t(pt) << std::endl;
 
   
  for (unsigned int i=0; i< 4; ++i) {
    for (unsigned int j=0; j< 4; ++j) {
      std::cout << t.hm(i,j) << " ";
    }
    std::cout << std::endl;
  }
  std::ifstream in(input_file.c_str());
  if (!in) {
    std::cerr << "Error opening input file " << input_file << std::endl;
    return EXIT_FAILURE;
  }

  PDB pdb(in, verbose);

  if (verbose) std::cout << "Input PDB has " << CGAL::PDB::distance(pdb.models()) << " models." << std::endl;
 

  CGAL_PDB_FOREACH(PDB::Model_pair &m, pdb.models()) {
    std::cout << "Model " << m.key() << " has " << CGAL::PDB::distance(m.model().chains()) 
              << " chains."<< std::endl;
    CGAL_PDB_FOREACH(Chain &c, make_chain_range(m.model().chains())) {
      CGAL_PDB_FOREACH(Monomer &r, make_monomer_range(c.monomers())) {
        CGAL_PDB_FOREACH(Atom &a, make_atom_range(r.atoms())) {
	  a.set_point(t(a.point()));
	}
      }
    }
    CGAL_PDB_FOREACH(Heterogen &h, make_heterogen_range(m.model().heterogens())) {
      CGAL_PDB_FOREACH(Atom &a, make_atom_range(h.atoms())) {
        a.set_point(t(a.point()));
      }
    }
  }


  std::ofstream out(output_file.c_str());
  if (!out) {
    std::cerr << "Error opening output file " << output_file << std::endl;
    return EXIT_FAILURE;
  }

  pdb.write(out);

  return EXIT_SUCCESS;
}
