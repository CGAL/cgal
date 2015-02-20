// Copyright (c) 2014  INRIA Sophia-Antipolis (France)
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin


#ifndef SIMPLICIAL_COMPLEX_H
#define SIMPLICIAL_COMPLEX_H

#include <CGAL/Tangential_complex/config.h>
#include <CGAL/Tangential_complex/utilities.h>

#include <CGAL/basic.h>

namespace CGAL {
namespace Tangential_complex_ {

class Simplicial_complex
{
public:
  typedef std::set<std::size_t>                     Simplex;
  typedef std::set<Simplex>                         Simplex_range;

  void add_simplex(const std::set<std::size_t> &s)
  {
    m_complex.insert(s);
  }

  const Simplex_range &simplex_range() const
  {
    return m_complex;
  }

  void collapse(int max_simplex_dim)
  {
    // We note k = max_simplex_dim - 1
    int k = max_simplex_dim - 1;
    
    typedef Complex::iterator                         Simplex_iterator;
    typedef std::vector<Simplex_iterator>             Simplex_iterator_list;
    typedef std::map<Simplex, Simplex_iterator_list>  Cofaces_map;

    // Create a map associating each non-maximal k-faces to the list of its
    // maximal cofaces
    Cofaces_map cofaces_map;
    for (Complex::const_iterator it_simplex = m_complex.begin(), 
                                 it_simplex_end = m_complex.end() ;
         it_simplex != it_simplex_end ; 
         ++it_simplex)
    {
      if (it_simplex->size() > k + 1)
      {
        std::vector<Simplex> k_faces;
        // Get the k-faces composing the simplex
        combinations(*it_simplex, k + 1, std::back_inserter(k_faces));
        for (const auto &comb : k_faces) // CJTODO C++1
          cofaces_map[comb].push_back(it_simplex);
      }
    }

    // For each non-maximal k-face F, if F has only one maximal coface Cf:
    // - Look for the other k-faces F2 of Cf in the map and:
    //    * if the list contains only Cf, clear the list (we don't remove the
    //      list since it was create troubles with the iterators) and add f2
    //      to the complex
    //    * otherwise, remove Cf from the associated list
    // - Remove Cf from the complex
    for (Cofaces_map::const_iterator it_map_elt = cofaces_map.begin(),
                                     it_map_end = cofaces_map.end() ;
         it_map_elt != it_map_end ;
         ++it_map_elt)
    {
      if (it_map_elt->second.size() == 1)
      {
        std::vector<Simplex> k_faces;
        const Simplex_iterator_list::value_type &it_Cf = 
          *it_map_elt->second.begin();
        CGAL_assertion(it_Cf->size() == max_simplex_dim + 1);
        // Get the k-faces composing the simplex
        combinations(*it_Cf, k + 1, std::back_inserter(k_faces));
        for (const auto &f2 : k_faces) // CJTODO C++1
        {
          if (f2 != it_map_elt->first)
          {
            Cofaces_map::iterator it_comb_in_map = cofaces_map.find(f2);
            if (it_comb_in_map->second.size() == 1)
            {
              it_comb_in_map->second.clear();
              m_complex.insert(f2);
            }
            else // it_comb_in_map->second.size() > 1
            {
              Simplex_iterator_list::iterator it = std::find(
                it_comb_in_map->second.begin(),
                it_comb_in_map->second.end(), 
                it_Cf);
              CGAL_assertion(it != it_comb_in_map->second.end());
              it_comb_in_map->second.erase(it);                
            }
          }
        }
        m_complex.erase(it_Cf);
      }
    }

    // Collapse the lower dimension simplices
    if (k > 1)
      collapse(max_simplex_dim - 1);
  }

  void display_stats() const
  {
    std::cerr << "==========================================================\n";
    std::cerr << "Complex stats:\n";      
    
    if (m_complex.empty())
    {
      std::cerr << "No simplices.\n";
    }
    else
    {
      // Number of simplex for each dimension
      std::map<int, std::size_t> simplex_stats;
    
      typedef std::set<std::size_t>       Simplex;
      typedef std::set<Simplex>           Complex;

      for (Complex::const_iterator it_simplex = m_complex.begin(), 
                                   it_simplex_end = m_complex.end() ;
           it_simplex != it_simplex_end ; 
           ++it_simplex)
      {
        ++simplex_stats[static_cast<int>(it_simplex->size()) - 1];
      }
   
      for (std::map<int, std::size_t>::const_iterator it_map = simplex_stats.begin() ;
           it_map != simplex_stats.end() ; ++it_map)
      {
        std::cerr << "  * " << it_map->first << "-simplices: " 
                  << it_map->second << std::endl;
      }
    }

    std::cerr << "==========================================================\n";
  }

  // CJTODO: this is only based on the fact that each edge should have two
  // cofaces
  bool is_manifold(bool exit_at_the_first_problem = true)
  {
    typedef std::set<std::size_t>               Edge;
    typedef std::map<Edge, std::size_t>         Cofaces_map;
    
    // Counts the number of cofaces of each edge
    
    // Create a map associating each non-maximal k-faces to the list of its
    // maximal cofaces
    Cofaces_map cofaces_map;
    for (Complex::const_iterator it_simplex = m_complex.begin(), 
                                 it_simplex_end = m_complex.end() ;
         it_simplex != it_simplex_end ; 
         ++it_simplex)
    {
      if (it_simplex->size() > 2)
      {
        std::vector<Edge> edges;
        // Get the k-faces composing the simplex
        combinations(*it_simplex, 2, std::back_inserter(edges));
        for (const auto &edge : edges) // CJTODO C++1
        {
          ++cofaces_map[edge];
        }
      }
    }

    bool manifold = true;
    for (Cofaces_map::const_iterator it_map_elt = cofaces_map.begin(),
                                     it_map_end = cofaces_map.end() ;
         it_map_elt != it_map_end ;
         ++it_map_elt)
    {
      if (it_map_elt->second != 2)
      {
        std::cerr << "Found an edge with " << it_map_elt->second << " cofaces\n";
        if (exit_at_the_first_problem)
          return false;
        else
          manifold = false;
      }
    }

    return manifold;
  }

private:
  typedef Simplex_range Complex;

  Complex       m_complex;

}; // /class Simplicial_complex

} // namespace Tangential_complex_
} //namespace CGAL

#endif // SIMPLICIAL_COMPLEX_H
