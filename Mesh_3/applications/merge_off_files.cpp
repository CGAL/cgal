// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Laurent Rineau

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stack>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/format.hpp>
#include <cstdlib>
#include <numeric> // for std::accumulate

using std::cout;
using std::cin;
using std::endl;

using std::ifstream;
std::vector<ifstream*> inputs;

int main(int argc, char** argv)
{
  const int number_of_inputs = argc-1;
  inputs.reserve(number_of_inputs);
  for(int i = 1; i < argc; ++i)
  {
    ifstream* input = new ifstream(argv[i]);
    if(!*input)
    {
      std::cerr << "Cannot open file \"" << argv[i] << "\"!\n";
      return EXIT_FAILURE;
    }
    else
      inputs.push_back(input);
  }

  std::vector<unsigned int> n_vertices(number_of_inputs);
  std::vector<unsigned int> n_facets(number_of_inputs);
  std::string header;
  for(int i = 0; i < number_of_inputs; ++i)
  {
    *(inputs[i]) >> header;
    
    if(header != "OFF")
    {
      std::cerr << "In file \"" << argv[i+1]
                << "\", header is \"" << header << "\"\nshould be \"OFF\".\n";
      return EXIT_FAILURE;
    }
  }

  cout << header << endl;
  std::string dummy;

  for(int i = 0; i < number_of_inputs; ++i)
  {
    *(inputs[i]) >> n_vertices[i];
    *(inputs[i]) >> n_facets[i];
    getline(*(inputs[i]), dummy);
  }

  const unsigned int total_n_vertices = 
    std::accumulate(n_vertices.begin(), n_vertices.end(), 0);
  const unsigned int total_n_facets = 
    std::accumulate(n_facets.begin(), n_facets.end(), 0);

  cout << total_n_vertices << " " << total_n_facets << " 0\n";

  for(int i = 0; i < number_of_inputs; ++i)
    for(unsigned int j = 0; j < n_vertices[i]; ++j)
    {
      std::string x, y, z;
      *(inputs[i]) >> x >> y >> z; // read x,y,z and write them to cout
      cout << x << " " << y << " " << z << std::endl;
      getline(*(inputs[i]), dummy); // ignore the end of the line
    }
  
  unsigned int vertex_index_offset = 0;
  for(int i_input = 0; i_input < number_of_inputs; ++i_input)
  {
    for(unsigned int i_facet = 0; i_facet < n_facets[i_input]; ++i_facet)
    {
      // Read a facet, then reindex its vertices.
      int i, j, k;
      *(inputs[i_input]) >> dummy >> i >> j >> k;
      if( dummy != "3" )
      {
        std::cerr << "In file \"" << argv[i_input+1] << "\""
                  << ", in facet #" << i_facet << ", expected \"3\", found \""
                  << dummy << "\"!\n";
        return 1;
      }
      getline(*(inputs[i_input]), dummy); // ignore the end of the line
      i += vertex_index_offset;
      j += vertex_index_offset;
      k += vertex_index_offset;
      cout << "3 " << i << " " << j << " " << k << "\n";
    }
    vertex_index_offset += n_vertices[i_input];
  }
  for(int i = 0; i < number_of_inputs; ++i)
  {
    if(!*(inputs[i]))
      return EXIT_FAILURE;
  }
  if(!cout) 
    return EXIT_FAILURE;
  else
    return EXIT_SUCCESS;
}
