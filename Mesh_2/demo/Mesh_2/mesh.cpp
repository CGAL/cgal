// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Laurent Rineau

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_area_criteria_2.h>
#include <CGAL/IO/File_poly.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds,
  CGAL::Exact_predicates_tag> Tr;
typedef CGAL::Delaunay_mesh_area_criteria_2<Tr> Meshcriteria;

typedef K::Point_2 Point;

void usage(char** argv)
{
  std::cerr << "Usage: " << std::endl
	    << argv[0] << " [-Q] [-a area] input.poly [output.poly]" << std::endl;
}

int main(int argc, char** argv)
{
  int arg_count = 1;
  bool terminal_output = true;
  Meshcriteria criteria;
  Tr t;

  if(argc < 2)
    {
      usage(argv);
      return 1;
    }

  while(argv[arg_count][0] == '-' && std::string(argv[arg_count]) != "--")
    {
      if(std::string(argv[arg_count]) == "-Q")
	terminal_output = false;
      else if(std::string(argv[arg_count]) == "-a")
	{
	  double area_bound;
	    if( (argc > arg_count+1) &&
		std::istringstream(argv[arg_count+1]) >> area_bound )
	      {
		criteria.set_area_bound(area_bound);
		++arg_count;
	      }
	    else
	      {
		std::cerr << "The area " << argv[arg_count+1]
			  << " is not a double." << std::endl;
		usage(argv);
		return 1;
	      }
	}
      else
	{
	  std::cerr << "Unknown option " << argv[arg_count] << std::endl;
	  usage(argv);
	  return 1;
	}
      ++arg_count;
    }
  if(std::string(argv[arg_count]) == "--")
    ++arg_count;

  if(argc < arg_count+1 || argc > arg_count+2)
    {
      usage(argv);
      return 1;
    }

  std::ifstream input(argv[arg_count]);
  if(input)
    {
      CGAL::read_triangle_poly_file(t, input);
      CGAL::refine_Delaunay_mesh_2(t, criteria);

      if(argc==arg_count+1)
	{
	  if(terminal_output)
	    CGAL::write_triangle_poly_file(t, std::cout);
	}
      else
	{
	  std::ofstream output(argv[arg_count+1]);
	  CGAL::write_triangle_poly_file(t, output);
	}
      if(terminal_output)
	std::cerr
	  << "Mesh points: " << t.number_of_vertices() << std::endl
	  << "Mesh triangles: " << t.number_of_faces () << std::endl;

    }
  else
    {
      std::cerr << "Bad file: " << argv[arg_count] << std::endl;
      usage(argv);
      return 1;
    }

  return 0;
}

