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

// return the cRMS
double align_to_points(const std::vector<CGAL_PDB_NS::Point> &points,
		       char mode, CGAL_PDB_NS::Protein &p){
  std::vector<CGAL_PDB_NS::Point> ppoints;
  switch (mode) {
  case 'b':
    CGAL_PDB_NS::backbone_coordinates(p.atoms_begin(), p.atoms_end(),
				 std::back_inserter(ppoints));
    break;
  case 'c':
    CGAL_PDB_NS::ca_coordinates(p.atoms_begin(), p.atoms_end(),
			   std::back_inserter(ppoints));
    break;
  case 'a':
    CGAL_PDB_NS::coordinates(p.atoms_begin(), p.atoms_end(),
			std::back_inserter(ppoints));
    break;
  };

  assert(points.size() == ppoints.size());

  CGAL_PDB_NS::Transform tr= CGAL_PDB_NS::transform_taking_first_to_second(ppoints, points);
  //std::cout << tr << std::endl;

  for (CGAL_PDB_NS::Protein::Atoms_iterator it= p.atoms_begin(); it != p.atoms_end(); ++it){
    CGAL_PDB_NS::Point np= tr(it->second.cartesian_coords());
    it->second.set_cartesian_coords(np);
  }

  /*
  {

   CGAL_PDB_NS::Transform tr2= CGAL_PDB_NS::compute_transform_taking_first_to_second(points,
									   ppoints);
   std::cout << tr2 << std::endl;

  }

  {
    std::vector<CGAL_PDB_NS::Point> npoints(ppoints.size());
    for (unsigned int i=0; i< points.size(); ++i){
      CGAL_PDB_NS::Point tp = tr(ppoints[i]);
      npoints[i]=tp;
    }
    CGAL_PDB_NS::Transform tr2= CGAL_PDB_NS::compute_transform_taking_first_to_second(npoints,
									    points);
    std::cout << tr2 << std::endl;

  }

  {
    std::vector<CGAL_PDB_NS::Point> npoints(ppoints.size());
    for (unsigned int i=0; i< points.size(); ++i){
      CGAL_PDB_NS::Point tp = tr(points[i]);
      npoints[i]=tp;
    }
    CGAL_PDB_NS::Transform tr2= CGAL_PDB_NS::compute_transform_taking_first_to_second(npoints,
									    ppoints);
    std::cout << tr2 << std::endl;

    }*/

  double dist=0;

  CGAL_PDB_NS::Squared_distance sd;

  for (unsigned int i=0; i< points.size(); ++i){
    CGAL_PDB_NS::Point tp = tr(ppoints[i]);
    dist += std::sqrt(sd(tp, points[i]));
  }
  return dist/ points.size();
}



int main(int argc, char *argv[]){
  std::string base_file, input_file, output_file;
  bool print_help=false;
  bool crms=false;
  char mode='b';
  bool warn=false;
  {
    boost::program_options::options_description o("Allowed options"), po, ao;
    o.add_options()
      ("aligned-pdb,o", boost::program_options::value< std::string>(&output_file),
       "Where to write the result of aligning input-pdb to base-points.")
      ("atoms,a", boost::program_options::value< char>(&mode),
       "Which atoms to use from the protein. Values are b (backbone), c (c-alpha), a (all).")
      ("crms,c", boost::program_options::bool_switch(&crms),
       "Output the cRMS between the two pdbs (after alignment).")
      ("verbose,v", boost::program_options::bool_switch(&warn),
       "Warn about errors parsing pdb file.")
      ("help", boost::program_options::bool_switch(&print_help),
       "Produce help message");
    po.add_options()("base-points", boost::program_options::value< std::string>(&base_file),
		     "Base points")
      ("input-pdb", boost::program_options::value< std::string>(&input_file),
       "input file");
    ao.add(o).add(po);

    boost::program_options::positional_options_description p;
    p.add("base-points", 1);
    p.add("input-pdb", 1);

    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).
				  options(ao).positional(p).run(), vm);
    boost::program_options::notify(vm);


    if (base_file.empty() || input_file.empty() || print_help
	|| (mode != 'b' && mode != 'c' && mode != 'a')) {
      std::cout << "This program aligns aligns a protein with a set of points and"
	" can compute distances between them.\n";
      std::cout << "usage: " << argv[0] << " base-points input-pdb\n\n";
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

  std::vector<CGAL_PDB_NS::Point> points;
  while (bin) {
    char buf[10000];
    bin.getline(buf, 10000);
    if (!bin) break;
    if (buf[0]=='#') continue;
    std::istringstream iss(buf);
    CGAL_PDB_NS::Point pt;
    iss >> pt;
    points.push_back(pt);
  }

  std::cout << "Read " << points.size() << " points.\n";

  if ( !output_file.empty() || crms) {
    for (unsigned int i=0; i< input.number_of_models(); ++i){
      CGAL_PDB_NS::Model &m= input.model(i);
      CGAL_PDB_NS::Protein &p= m.chain(0);
      double cRMS= align_to_points(points, mode, p);
      if (crms) {
	std::cout << cRMS << std::endl;
      }
    }
  }

  if (!output_file.empty()) {
    std::ofstream out(output_file.c_str());
    if (!out){
      std::cerr << "Error opening output file "
		<< output_file << std::endl;
      return EXIT_FAILURE;
    }
    input.write(out);
  }



  return EXIT_SUCCESS;

}
