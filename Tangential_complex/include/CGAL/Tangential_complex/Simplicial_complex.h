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
#include "CGAL/Tangential_complex/console_color.h"

#include <CGAL/basic.h>
#include <CGAL/iterator.h>

// For is_pure_pseudomanifold
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/container/flat_map.hpp>

#include <algorithm>
#include <string>
#include <fstream>

namespace CGAL {
namespace Tangential_complex_ {

class Simplicial_complex
{
public:
  typedef std::set<std::size_t>                     Simplex;
  typedef std::set<Simplex>                         Simplex_set;

  // If perform_checks = true, the function:
  // - won't insert the simplex if it is already in a higher dim simplex
  // - will erase any lower-dim simplices that are faces of the new simplex
  // Returns true if the simplex was added
  bool add_simplex(
    const std::set<std::size_t> &s, bool perform_checks = true)
  {
    if (perform_checks)
    {
      unsigned int num_pts = static_cast<int>(s.size());
      std::vector<Complex::iterator> to_erase;
      bool check_higher_dim_simpl = true;
      for (Complex::iterator it_simplex = m_complex.begin(), 
                             it_simplex_end = m_complex.end() ;
           it_simplex != it_simplex_end ; 
           ++it_simplex)
      {
        // Check if the simplex is not already in a higher dim simplex
        if (check_higher_dim_simpl 
          && it_simplex->size() > num_pts 
          && std::includes(it_simplex->begin(), it_simplex->end(),
                           s.begin(), s.end()))
        {
          // No need to insert it, then
          return false;
        }
        // Check if the simplex includes some lower-dim simplices
        if (it_simplex->size() < num_pts 
          && std::includes(s.begin(), s.end(),
                           it_simplex->begin(), it_simplex->end()))
        {
          to_erase.push_back(it_simplex);
          // We don't need to check higher-sim simplices any more
          check_higher_dim_simpl = false;
        }
      }
      for (std::vector<Complex::iterator>::const_iterator it= to_erase.begin();
           it != to_erase.end() ; ++it)
      {
        m_complex.erase(*it);
      }
    }
    return m_complex.insert(s).second;
  }

  // Ignore the point coordinates
  bool load_simplices_from_OFF(const char *filename)
  {
    std::ifstream in(filename);
    if (!in.is_open())
    {
      std::cerr << "Could not open '" << filename << "'" << std::endl;
      return false;
    }

    std::string line;
    std::getline(in, line); // Skip first line
    std::size_t num_points, num_triangles, dummy;
    in >> num_points;
    in >> num_triangles;
    std::getline(in, line); // Skip the rest of the line
    
    // Skip points
    for (std::size_t i = 0 ; i < num_points ; ++i)
      std::getline(in, line);

    // Read triangles
    for (std::size_t i = 0 ; i < num_triangles ; ++i)
    {
      in >> dummy;
      CGAL_assertion(dummy == 3);
      std::set<std::size_t> tri;
      std::size_t idx;
      for (int j = 0 ; j < 3 ; ++j)
      {
        in >> idx;
        tri.insert(idx);
      }
      std::getline(in, line); // Skip the rest of the line

      add_simplex(tri, false);
    }

    return true;
  }

  const Simplex_set &simplex_range() const
  {
    return m_complex;
  }

  bool empty()
  {
    return m_complex.empty();
  }

  void clear()
  {
    m_complex.clear();
  }

  template <typename Test, typename Output_it>
  void get_simplices_matching_test(Test test, Output_it out)
  {
    for (Complex::const_iterator it_simplex = m_complex.begin(), 
                                 it_simplex_end = m_complex.end() ;
         it_simplex != it_simplex_end ; 
         ++it_simplex)
    {
      if (test(*it_simplex))
        *out++ = *it_simplex;
    }
  }

  // When a simplex S has only one co-face C, we can remove S and C 
  // without changing the topology
  void collapse(int max_simplex_dim, bool quiet = false)
  {
#ifdef CGAL_TC_VERBOSE
    if (!quiet)
      std::cerr << "Collapsing... ";
#endif
    // We note k = max_simplex_dim - 1
    int k = max_simplex_dim - 1;

    typedef Complex::iterator                         Simplex_iterator;
    typedef std::vector<Simplex_iterator>             Simplex_iterator_list;
    typedef std::map<Simplex, Simplex_iterator_list>  Cofaces_map;

    std::size_t num_collapsed_maximal_simplices = 0;
    do 
    {
      num_collapsed_maximal_simplices = 0;
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
      // - Look for the other k-faces F2, F3... of Cf in the map and:
      //    * if the list contains only Cf, clear the list (we don't remove the
      //      list since it creates troubles with the iterators) and add the F2,
      //      F3... to the complex
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
            // Skip F
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
          ++num_collapsed_maximal_simplices;
        }
      }
      // Repeat until no maximal simplex got removed
    } while (num_collapsed_maximal_simplices > 0);

    // Collapse the lower dimension simplices
    if (k > 0)
      collapse(max_simplex_dim - 1, true);

#ifdef CGAL_TC_VERBOSE
    if (!quiet)
      std::cerr << "done.\n";
#endif
  }

  void display_stats() const
  {
    std::cerr << yellow << "Complex stats:\n" << white;

    if (m_complex.empty())
    {
      std::cerr << "  * No simplices.\n";
    }
    else
    {
      // Number of simplex for each dimension
      std::map<int, std::size_t> simplex_stats;

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
                  << it_map->second << "\n";
      }
    }
  }

  // verbose_level = 0, 1 or 2
  bool is_pure_pseudomanifold__do_not_check_if_stars_are_connected(
    int simplex_dim,
    bool allow_borders = false,
    bool exit_at_the_first_problem = false,
    int verbose_level = 0,
    std::size_t *p_num_wrong_dim_simplices = NULL, 
    std::size_t *p_num_wrong_number_of_cofaces = NULL) const
  {
    typedef std::set<std::size_t>               K_1_face;
    typedef std::map<K_1_face, std::size_t>     Cofaces_map;

    std::size_t num_wrong_dim_simplices = 0;
    std::size_t num_wrong_number_of_cofaces = 0;

    // Counts the number of cofaces of each K_1_face

    // Create a map associating each non-maximal k-faces to the list of its
    // maximal cofaces
    Cofaces_map cofaces_map;
    for (Complex::const_iterator it_simplex = m_complex.begin(),
                                 it_simplex_end = m_complex.end() ;
         it_simplex != it_simplex_end ;
         ++it_simplex)
    {
      if (it_simplex->size() != simplex_dim + 1)
      {
        if (verbose_level >= 2)
          std::cerr << "Found a simplex with dim = "
                    << it_simplex->size() - 1 << "\n";            
        ++num_wrong_dim_simplices;
      }
      else
      {
        std::vector<K_1_face> k_1_faces;
        // Get the facets composing the simplex
        combinations(
          *it_simplex, simplex_dim, std::back_inserter(k_1_faces));
        for (const auto &k_1_face : k_1_faces) // CJTODO C++1
        {
          ++cofaces_map[k_1_face];
        }
      }
    }

    for (Cofaces_map::const_iterator it_map_elt = cofaces_map.begin(),
                                     it_map_end = cofaces_map.end() ;
         it_map_elt != it_map_end ;
         ++it_map_elt)
    {
      if (it_map_elt->second != 2
        && (!allow_borders || it_map_elt->second != 1))
      {
        if (verbose_level >= 2)
          std::cerr << "Found a k-1-face with "
                    << it_map_elt->second << " cofaces\n";

        if (exit_at_the_first_problem)
          return false;
        else
          ++num_wrong_number_of_cofaces;
      }
    }

    bool ret = num_wrong_dim_simplices == 0 && num_wrong_number_of_cofaces == 0;

    if (verbose_level >= 1)
    {
      std::cerr << "Pure pseudo-manifold: ";
      if (ret)
      {
        std::cerr << green << "YES" << white << "\n";
      }
      else
      {
        std::cerr << red << "NO" << white << "\n"
          << "  * Number of wrong dimension simplices: "
          << num_wrong_dim_simplices << "\n"
          << "  * Number of wrong number of cofaces: "
          << num_wrong_number_of_cofaces << "\n";
      }
    }

    if (p_num_wrong_dim_simplices)
      *p_num_wrong_dim_simplices = num_wrong_dim_simplices;
    if (p_num_wrong_number_of_cofaces)
      *p_num_wrong_number_of_cofaces = num_wrong_number_of_cofaces;

    return ret;
  }

  template <int K>
  std::size_t num_K_simplices() const
  {
    std::set<Simplex> k_simplices; 

    for (Complex::const_iterator it_simplex = m_complex.begin(),
      it_simplex_end = m_complex.end() ;
      it_simplex != it_simplex_end ;
    ++it_simplex)
    {
      if (it_simplex->size() == K + 1)
      {
        k_simplices.insert(*it_simplex);
      }
      else if (it_simplex->size() > K + 1)
      {
        // Get the k-faces composing the simplex
        combinations(
          *it_simplex, K + 1, std::inserter(k_simplices, k_simplices.begin()));
      }
    }

    return k_simplices.size();
  }

  std::ptrdiff_t euler_characteristic(bool verbose = false) const
  {
    if (verbose)
      std::cerr << "\nComputing Euler characteristic of the complex...\n";
    
    std::size_t num_vertices = num_K_simplices<0>();
    std::size_t num_edges = num_K_simplices<1>();
    std::size_t num_triangles = num_K_simplices<2>();
    
    std::ptrdiff_t ec =
      (std::ptrdiff_t) num_vertices
      - (std::ptrdiff_t) num_edges
      + (std::ptrdiff_t) num_triangles;

    if (verbose)
      std::cerr << "Euler characteristic: V - E + F = "
        << num_vertices << " - " << num_edges << " + " << num_triangles << " = "
        << blue
        << ec
        << white << "\n";

    return ec;
  }

  // CJTODO: ADD COMMENTS
  bool is_pure_pseudomanifold(
    int simplex_dim,
    std::size_t num_vertices,
    bool allow_borders = false,
    bool exit_at_the_first_problem = false, 
    int verbose_level = 0,
    std::size_t *p_num_wrong_dim_simplices = NULL, 
    std::size_t *p_num_wrong_number_of_cofaces = NULL,
    std::size_t *p_num_unconnected_stars = NULL,
    Simplex_set *p_wrong_dim_simplices = NULL,
    Simplex_set *p_wrong_number_of_cofaces_simplices = NULL,
    Simplex_set *p_unconnected_stars_simplices = NULL) const
  {
    // If simplex_dim == 1, we do not need to check if stars are connected
    if (simplex_dim == 1)
    {
      if (p_num_unconnected_stars)
        *p_num_unconnected_stars = 0;
      return is_pure_pseudomanifold__do_not_check_if_stars_are_connected(
        simplex_dim,
        allow_borders,
        exit_at_the_first_problem,
        verbose_level,
        p_num_wrong_dim_simplices,
        p_num_wrong_number_of_cofaces);
    }
    // Associates each vertex (= the index in the vector) 
    // to its star (list of simplices)
    typedef std::vector<std::vector<Complex::const_iterator> > Stars;
    std::size_t num_wrong_dim_simplices = 0;
    std::size_t num_wrong_number_of_cofaces = 0;
    std::size_t num_unconnected_stars = 0;

    // Fills a Stars data structure
    Stars stars;
    stars.resize(num_vertices);
    for (Complex::const_iterator it_simplex = m_complex.begin(), 
                                 it_simplex_end = m_complex.end() ;
         it_simplex != it_simplex_end ; 
         ++it_simplex)
    {
      if (it_simplex->size() != simplex_dim + 1)
      {
        if (verbose_level >= 2)
          std::cerr << "Found a simplex with dim = " 
                    << it_simplex->size() - 1 << "\n";
        ++num_wrong_dim_simplices;
        if (p_wrong_dim_simplices)
          p_wrong_dim_simplices->insert(*it_simplex);
      }
      else
      {
        for (Simplex::const_iterator it_point_idx = it_simplex->begin() ; 
             it_point_idx != it_simplex->end() ; 
             ++it_point_idx)
        {
          stars[*it_point_idx].push_back(it_simplex);
        }
      }
    }

    // Now, for each star, we have a vector of its d-simplices
    // i.e. one index for each d-simplex
    // Boost Graph only deals with indexes, so we also need indexes for the
    // (d-1)-simplices
    std::size_t center_vertex_index = 0;
    for (Stars::const_iterator it_star = stars.begin() ;
         it_star != stars.end() ;
         ++it_star, ++center_vertex_index)
    {
      typedef std::map<Simplex, std::vector<std::size_t> > 
                                                      Dm1_faces_to_adj_D_faces;
      Dm1_faces_to_adj_D_faces dm1_faces_to_adj_d_faces;
      
      for (int i_dsimpl = 0 ; i_dsimpl < it_star->size() ; ++i_dsimpl)
      {
        Simplex dm1_simpl_of_link = *((*it_star)[i_dsimpl]);
        dm1_simpl_of_link.erase(center_vertex_index);
        // Copy it to a vector so that we can use operator[] on it
        std::vector<std::size_t> dm1_simpl_of_link_vec(
          dm1_simpl_of_link.begin(), dm1_simpl_of_link.end());

        CGAL::Combination_enumerator<int> dm2_simplices(
          simplex_dim - 1, 0, simplex_dim);
        for ( ; !dm2_simplices.finished() ; ++dm2_simplices)
        {
          Simplex dm2_simpl;
          for (int j = 0 ; j < simplex_dim - 1 ; ++j)
            dm2_simpl.insert(dm1_simpl_of_link_vec[dm2_simplices[j]]);
          dm1_faces_to_adj_d_faces[dm2_simpl].push_back(i_dsimpl);
        }
      }

      Adj_graph adj_graph;
      std::vector<Graph_vertex> d_faces_descriptors;
      d_faces_descriptors.resize(it_star->size());
      for (int j = 0 ; j < it_star->size() ; ++j)
        d_faces_descriptors[j] = boost::add_vertex(adj_graph);

      Dm1_faces_to_adj_D_faces::const_iterator dm1_to_d_it = 
        dm1_faces_to_adj_d_faces.begin();
      Dm1_faces_to_adj_D_faces::const_iterator dm1_to_d_it_end = 
        dm1_faces_to_adj_d_faces.end();
      for (std::size_t i_km1_face = 0 ; 
           dm1_to_d_it != dm1_to_d_it_end ; 
           ++dm1_to_d_it, ++i_km1_face)
      {
        Graph_vertex km1_gv = boost::add_vertex(adj_graph);

        for (std::vector<std::size_t>::const_iterator kface_it = 
               dm1_to_d_it->second.begin() ;
             kface_it != dm1_to_d_it->second.end() ;
             ++kface_it)
        {
          boost::add_edge(km1_gv, *kface_it, adj_graph);
        }

        if (dm1_to_d_it->second.size() != 2
          && (!allow_borders || dm1_to_d_it->second.size() != 1))
        {
          ++num_wrong_number_of_cofaces;
          if (p_wrong_number_of_cofaces_simplices)
          {
            for (auto idx : dm1_to_d_it->second)
              p_wrong_number_of_cofaces_simplices->insert(*((*it_star)[idx]));
          }
        }
      }

      // What is left is to check the connexity
      bool is_connected = true;
      if (boost::num_vertices(adj_graph) > 0)
      {
        std::vector<int> components(boost::num_vertices(adj_graph));
        is_connected =
          (boost::connected_components(adj_graph, &components[0]) == 1);
      }

      if (!is_connected)
      { 
        if (verbose_level >= 2)
          std::cerr << "Error: star #" << center_vertex_index
                    << " is not connected\n";
        ++num_unconnected_stars;
        if (p_unconnected_stars_simplices)
        {
          for (std::vector<Complex::const_iterator>::const_iterator 
                it_simpl = it_star->begin(),
                it_simpl_end = it_star->end();
               it_simpl != it_simpl_end ;
               ++it_simpl)
          {
            p_unconnected_stars_simplices->insert(**it_simpl);
          }
        }
      }
    }

    // Each one has been counted several times ("simplex_dim" times)
    num_wrong_number_of_cofaces /= simplex_dim;

    bool ret = 
      num_wrong_dim_simplices == 0 
      && num_wrong_number_of_cofaces == 0
      && num_unconnected_stars == 0;

    if (verbose_level >= 1)
    {
      std::cerr << "Pure pseudo-manifold: ";
      if (ret)
      {
        std::cerr << green << "YES" << white << "\n";
      }
      else
      {
        std::cerr << red << "NO" << white << "\n"
                  << "  * Number of wrong dimension simplices: " 
                  << num_wrong_dim_simplices << "\n"
                  << "  * Number of wrong number of cofaces: " 
                  << num_wrong_number_of_cofaces << "\n"
                  << "  * Number of not-connected stars: " 
                  << num_unconnected_stars << "\n";
      }
    }

    if (p_num_wrong_dim_simplices)
      *p_num_wrong_dim_simplices = num_wrong_dim_simplices;
    if (p_num_wrong_number_of_cofaces)
      *p_num_wrong_number_of_cofaces = num_wrong_number_of_cofaces;
    if (p_num_unconnected_stars)
      *p_num_unconnected_stars = num_unconnected_stars;

    return ret;
  }

private:
  typedef Simplex_set Complex;
 
  // graph is an adjacency list
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Adj_graph;
  // map that gives to a certain simplex its node in graph and its dimension
  typedef boost::graph_traits<Adj_graph>::vertex_descriptor Graph_vertex;
  typedef boost::graph_traits<Adj_graph>::edge_descriptor Graph_edge;

  Complex       m_complex;

}; // /class Simplicial_complex

} // namespace Tangential_complex_
} //namespace CGAL

#endif // SIMPLICIAL_COMPLEX_H
