// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Rineau

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stack>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/format.hpp>

using std::cout;
using std::cin;
using std::endl;

int main()
{
  std::string header;
  cin >> header;

  if(header != "OFF")
  {
    std::cerr << "header is \"" << header << "\"\nshould be \"OFF\".\n";
    return 1;
  }

  cout << header << endl;
  unsigned int n_vertices;
  unsigned int n_facets;
  std::string dummy;

  cin >> n_vertices;
  cin >> n_facets;
  getline(cin, dummy);

  // Vector that maps from old vertex index to new vertex index,
  // as some (old) vertices may have identical point coordinates.
  std::vector<int> new_index(n_vertices);

  typedef boost::tuple<double, double, double> Point;
  // Vector that stores the points coordinates (in a Point).
  std::vector<Point> points(n_vertices);

  // Map that retains the mapping from the point coordinates (in a Point)
  // to the nex index.
  typedef std::map<Point, int> Renumber;
  Renumber renumber;

  unsigned int index = 0;
  for(unsigned int i = 0; i < n_vertices; ++i)
  {
    cin >> boost::get<0>(points[index])
        >> boost::get<1>(points[index]) >> boost::get<2>(points[index]);
    Renumber::const_iterator it = renumber.find(points[index]);
    if( it == renumber.end() )
    {
      renumber[points[index]] = index;
      new_index[i] = index;
      ++index;
    }
    else
    {
      new_index[i] = it->second;
    }
  }

  cout << index << " " << n_facets << dummy << endl;
  for(unsigned int i = 0; i < index; ++i)
  {
    cout << boost::get<0>(points[i]) << " "
         << boost::get<1>(points[i]) << " "
         << boost::get<2>(points[i]) << "\n";
  }

  // Vector that stores each facet.
  typedef std::vector<boost::tuple<int, int, int> > Facets;
  Facets facets;

  // For each (oriented) edge, that map stores a vector of adjacent facets,
  // with a boolean that tells if the edge is in the opposite orientation,
  // in the facet.
  // Edges are stores in the direction from the smallest index to the
  // greatest.
  typedef std::pair<int, int> Edge;
  typedef std::map<Edge, std::vector<std::pair<int, bool> > > Edges_map;
  Edges_map edges;

  // "nested function" opposite, that returns the edge, in the opposite
  // direction.
  struct {
    Edge operator()(Edge e) const {
      return std::make_pair(e.second, e.first);
    };
  } opposite;

  for(unsigned int i_facet = 0; i_facet < n_facets; ++i_facet)
  {
    // Read a facet, then reindex its vertices.
    int i, j, k;
    cin >> dummy >> i >> j >> k;
    if( dummy != "3" )
    {
      std::cerr << "In facet #" << i_facet << ", expected \"3\", found \""
                << dummy << "\"!\n";
      return 1;
    }
    i = new_index[i];
    j = new_index[j];
    k = new_index[k];
    facets.push_back(boost::make_tuple(i, j, k));

    // Create the three edges of the facet.
    Edge e[3];
    e[0] = std::make_pair(i, j);
    e[1] = std::make_pair(j, k);
    e[2] = std::make_pair(k, i);
    for(int i_edge = 0; i_edge < 3; ++i_edge)
    {
      if( e[i_edge].first < e[i_edge].second )
        edges[e[i_edge]].push_back(std::make_pair(i_facet, false));
      else
        edges[opposite(e[i_edge])].push_back(std::make_pair(i_facet, true));
    }
  }

  // Map that stores all already passed facet, and retains the orientation
  // of the facet. "true" means that the facet needs to be reoriented.
  std::map<int, bool> oriented_set;

  // Stack of facets indices to be handled.
  std::stack<int> stack;
  int seed_facet_candidate = 0;

  while (oriented_set.size() != n_facets) {
    // find a facet index that is not yet in 'oriented_set'.
    while( oriented_set.find(seed_facet_candidate) != oriented_set.end() )
      ++seed_facet_candidate;
    std::cerr << "Need seed facet: " << seed_facet_candidate << "\n";
    // push it in oriented set
    oriented_set[seed_facet_candidate] = false;
    stack.push(seed_facet_candidate);

    while(! stack.empty() ) {
      const int f = stack.top();
      stack.pop();
      const int i = boost::get<0>(facets[f]);
      const int j = boost::get<1>(facets[f]);
      const int k = boost::get<2>(facets[f]);
      Edge e[3];
      e[0] = std::make_pair(i, j);
      e[1] = std::make_pair(j, k);
      e[2] = std::make_pair(k, i);
      for(int ih = 0 ; ih < 3 ; ++ih)
      {
        bool f_orient = false;
        if(e[ih].first > e[ih].second) {
          f_orient = true;
          e[ih] = opposite(e[ih]);
        }

        Edges_map::iterator edge_it = edges.find(e[ih]);
        if(edge_it->second.size() == 2) { // regular edge
          int fn = edge_it->second[0].first;
          bool fn_orient = edge_it->second[0].second;
          if(fn == f)
          {
            fn = edge_it->second[1].first;
            fn_orient = edge_it->second[1].second;
          }
          if (oriented_set.find(fn) == oriented_set.end())
          {
            if(f_orient == fn_orient)
              oriented_set[fn] = ! oriented_set[f];
            else
              oriented_set[fn] = oriented_set[f];
            stack.push(fn);
          }
        } // end "if the edge is regular"
        else {
          std::cerr << boost::format("Irregular edge: (%1%,%2%)"
                                     ", %3% facets.\n")
            % e[ih].first % e[ih].second
            % edge_it->second.size();
        }
      } // end "for each neighbor of f"
    } // end "stack non empty"
  } // end "oriented_set not full"

  for(unsigned int i_facet = 0; i_facet < n_facets; ++i_facet)
  {
    const int i = boost::get<0>(facets[i_facet]);
    const int j = boost::get<1>(facets[i_facet]);
    const int k = boost::get<2>(facets[i_facet]);
    if(oriented_set[i_facet])
      cout << "3 " << j << " " << i << " " << k << "\n";
    else
      cout << "3 " << i << " " << j << " " << k << "\n";
  }
}

