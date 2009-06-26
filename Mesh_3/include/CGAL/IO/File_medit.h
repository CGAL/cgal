// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Laurent RINEAU

#ifndef CGAL_IO_FILE_MEDIT_H
#define CGAL_IO_FILE_MEDIT_H

#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <CGAL/utility.h>

namespace CGAL {

namespace {


//-------------------------------------------------------
// Default property map implementation
//-------------------------------------------------------

#ifdef CGAL_MESH_3_IO_VERBOSE
template<class T>
inline
std::ostream&
operator<<(std::ostream &os, const std::pair<T,T>& pair)
{
  return os << "<" << pair.first << "," << pair.second << ">";
}
#endif


/**
 *
 */
template <typename C3T3>
class Default_facet_index_pmap
{
  typedef typename C3T3::Surface_index Surface_index;
  typedef std::map<Surface_index,int> Surface_map;
  typedef typename C3T3::Facet Facet;
  typedef unsigned int size_type;

public:
  Default_facet_index_pmap(const C3T3& c3t3, const int first_index = 0)
    : r_c3t3_(c3t3)
  {
    typedef typename C3T3::Facet_iterator Facet_iterator;

    int index_counter = first_index + 1;

    for( Facet_iterator facet_it = r_c3t3_.facets_begin();
         facet_it != r_c3t3_.facets_end();
         ++facet_it)
    {
      // Add surface index in internal map if needed
      if ( surface_map_.end() ==
               surface_map_.find(c3t3.surface_index((*facet_it).first,
                                                    (*facet_it).second)) )
      {
        surface_map_.insert(std::make_pair(r_c3t3_.surface_index(*facet_it),
                                           index_counter));
        ++index_counter;
      }
    }

#ifdef CGAL_MESH_3_IO_VERBOSE
    std::cerr << "Nb of surface patches: " << surface_map_.size() << "\n";
    std::cerr << "Surface mapping:\n\t" ;

    typedef typename Surface_map::iterator Surface_map_iterator;
    for ( Surface_map_iterator surf_it = surface_map_.begin() ;
          surf_it != surface_map_.end() ;
          ++surf_it )
    {
      std::cerr << "[" << (*surf_it).first
                << ":" << (*surf_it).second << "] ";
    }
    std::cerr << "\n";
#endif
  }

  int surface_index(const Facet& f) const
  {
    return surface_index(r_c3t3_.surface_index(f));
  }

  int surface_index(const Surface_index& index) const
  {
    typedef typename Surface_map::const_iterator Smi;
    Smi elt_it = surface_map_.find(index);
    if ( elt_it != surface_map_.end() )
      return elt_it->second;
    else
      return -1;
  }

  size_type surface_number() const
  {
    return surface_map_.size();
  }

private:
  const C3T3& r_c3t3_;
  Surface_map surface_map_;
};


template <typename C3T3>
class Default_cell_index_pmap
{
  typedef typename C3T3::Subdomain_index Subdomain_index;
  typedef std::map<Subdomain_index,int> Subdomain_map;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef unsigned int size_type;

public:
  Default_cell_index_pmap(const C3T3& c3t3, const int first_index = 0)
    : r_c3t3_(c3t3)
  {
    typedef typename C3T3::Cell_iterator Cell_iterator;

    int index_counter = first_index + 1;

    for( Cell_iterator cell_it = r_c3t3_.cells_begin();
         cell_it != r_c3t3_.cells_end();
         ++cell_it)
    {
      // Add subdomain index in internal map if needed
      if ( subdomain_map_.end() ==
              subdomain_map_.find(r_c3t3_.subdomain_index(cell_it)) )
      {
        subdomain_map_.insert(std::make_pair(r_c3t3_.subdomain_index(cell_it),
                                             index_counter));
        ++index_counter;
      }
    }

#ifdef CGAL_MESH_3_IO_VERBOSE
    std::cerr << "Nb of subdomains: " << subdomain_map_.size() << "\n";
    std::cerr << "Subdomain mapping:\n\t" ;

    typedef typename Subdomain_map::iterator Subdomain_map_iterator;
    for ( Subdomain_map_iterator sub_it = subdomain_map_.begin() ;
          sub_it != subdomain_map_.end() ;
          ++sub_it )
    {
      std::cerr << "[" << (*sub_it).first << ":" << (*sub_it).second << "] ";
    }
    std::cerr << "\n";
#endif
  }

  int subdomain_index(const Cell_handle& ch) const
  {
    return subdomain_index(r_c3t3_.subdomain_index(ch));
  }

  int subdomain_index(const Subdomain_index& index) const
  {
    typedef typename Subdomain_map::const_iterator Smi;
    Smi elt_it = subdomain_map_.find(index);
    if ( elt_it != subdomain_map_.end() )
      return elt_it->second;
    else
      return -1;
  }

  size_type subdomain_number() const
  {
    return subdomain_map_.size();
  }

private:
  const C3T3& r_c3t3_;
  Subdomain_map subdomain_map_;
};


template <typename C3T3>
class Default_vertex_index_pmap
{
  typedef typename C3T3::Surface_index Surface_index;
  typedef typename C3T3::Subdomain_index Subdomain_index;
  typedef typename C3T3::Index Index;
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef typename C3T3::Facet Facet;

public:
  Default_vertex_index_pmap(const Default_facet_index_pmap<C3T3>& f_pmap,
                            const Default_cell_index_pmap<C3T3>& c_pmap,
                            const C3T3& c3t3)
    : f_pmap_(f_pmap)
    , c_pmap_(c_pmap)
    , r_c3t3_(c3t3)
    , edge_index_(1 + f_pmap_.surface_number() + c_pmap_.subdomain_number()) {}

  int index(const Vertex_handle& vh) const
  {
    switch ( r_c3t3_.in_dimension(vh) )
    {
    case 2:
      {
        // Check if each incident surface facet of vh has the same surface index
        typename std::vector<Facet> facets;
        r_c3t3_.triangulation().finite_incident_facets(
            vh, std::back_inserter(facets));

        if ( facets.begin() == facets.end() )
          return -1;

        // Look for the first surface facet
        typename std::vector<Facet>::iterator it_facet = facets.begin();
        while ( ! r_c3t3_.is_in_complex(*it_facet) )
        {
          if ( ++it_facet == facets.end() )
            return -1;
        }

        Surface_index facet_index = r_c3t3_.surface_index(*it_facet);
        ++it_facet;

        for( ; it_facet != facets.end() ; ++it_facet)
        {
          // If another index is found, return value for edge vertice
          if (   r_c3t3_.is_in_complex(*it_facet)
              && facet_index != r_c3t3_.surface_index(*it_facet) )
            return edge_index_;
        }

        return f_pmap_.surface_index(facet_index);
      }
      break;

    case 3:
      {
        // Returns value of any incident cell
        typename std::vector<Cell_handle> cells;
        r_c3t3_.triangulation().finite_incident_cells(
            vh,std::back_inserter(cells));

        if ( cells.begin() != cells.end() )
          return c_pmap_.subdomain_index(r_c3t3_.subdomain_index(*cells.begin()));
        else
          return -1;
      }
      break;

    default:
      // should not happen
      return -1;
      break;
    }
  }

private:
  Default_facet_index_pmap<C3T3> f_pmap_;
  Default_cell_index_pmap<C3T3> c_pmap_;
  const C3T3& r_c3t3_;
  const unsigned int edge_index_;
};


//template <typename C3T3>
//class Vertex_index_pmap
//{
//  typedef typename C3T3::Surface_index Surface_index;
//  typedef typename C3T3::Subdomain_index Subdomain_index;
//  typedef typename C3T3::Index Index;
//  typedef typename C3T3::Vertex_handle Vertex_handle;
//
//public:
//  Vertex_index_pmap(const Default_facet_index_pmap<C3T3>& f_pmap,
//                    const Default_cell_index_pmap<C3T3>& c_pmap,
//                    const C3T3& c3t3)
//    : f_pmap_(f_pmap)
//    , c_pmap_(c_pmap)
//    , r_c3t3_(c3t3)       { }
//
//  int index(const Vertex_handle& vh) const
//  {
//    switch ( r_c3t3_.in_dimension(vh) )
//    {
//    case 2:
//      return f_pmap_.surface_index(r_c3t3_.index(vh).surface_index());
//      break;
//    case 3:
//      return c_pmap_.subdomain_index(r_c3t3_.index(vh).subdomain_index());
//      break;
//    default:
//      // should not happen
//      return -1;
//      break;
//    }
//  }
//
//private:
//  Default_facet_index_pmap<C3T3> f_pmap_;
//  Default_cell_index_pmap<C3T3> c_pmap_;
//  const C3T3& r_c3t3_;
//};

//-------------------------------------------------------
// Property map access methods
//-------------------------------------------------------
template <typename C3T3>
int get(const Default_vertex_index_pmap<C3T3>& vmap,
        const typename C3T3::Vertex_handle& vh)
{
  return vmap.index(vh);
}

template <typename C3T3>
int get(const Default_facet_index_pmap<C3T3>& fmap,
        const typename C3T3::Facet& f)
{
  return fmap.surface_index(f);
}

template <typename C3T3>
int get(const Default_cell_index_pmap<C3T3>& fmap,
        const typename C3T3::Cell_handle& ch)
{
  return fmap.subdomain_index(ch);
}

//-------------------------------------------------------
// IO functions
//-------------------------------------------------------

template <class C3T3>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3)
{
#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cerr << "Output to medit:\n";
#endif
  output_to_medit(os,
                  c3t3,
                  Default_cell_index_pmap<C3T3>(c3t3));
}

template <class C3T3>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3,
                const Default_cell_index_pmap<C3T3>& cell_pmap)
{
  output_to_medit(
      os,
      c3t3,
      Default_facet_index_pmap<C3T3>(c3t3, cell_pmap.subdomain_number()),
      cell_pmap);
}

template <class C3T3>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3,
                const Default_facet_index_pmap<C3T3>& facet_pmap,
                const Default_cell_index_pmap<C3T3>& cell_pmap)
{
  output_to_medit(
      os,
      c3t3,
      Default_vertex_index_pmap<C3T3>(facet_pmap, cell_pmap, c3t3),
      facet_pmap,
      cell_pmap);
}

template <class C3T3,
          class Vertex_index_property_map,
          class Facet_index_property_map,
          class Cell_index_property_map>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3,
                const Vertex_index_property_map& vertex_pmap,
                const Facet_index_property_map& facet_pmap,
                const Cell_index_property_map& cell_pmap)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Facet_iterator Facet_iterator;
  typedef typename C3T3::Cell_iterator Cell_iterator;

  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Point Point_3;

  const Tr& tr = c3t3.triangulation();

  //-------------------------------------------------------
  // File output
  //-------------------------------------------------------

  //-------------------------------------------------------
  // Header
  //-------------------------------------------------------
  os << std::setprecision(20);

  os << "MeshVersionFormatted 1" << std::endl
     << "Dimension 3" << std::endl;


  //-------------------------------------------------------
  // Vertices
  //-------------------------------------------------------
  os << "Vertices" << std::endl
     << tr.number_of_vertices() << std::endl;

  std::map<Vertex_handle, int> V;
  int inum = 1;
  for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
       vit != tr.finite_vertices_end();
       ++vit)
  {
    V[vit] = inum++;
    Point_3 p = vit->point();
    os << CGAL::to_double(p.x()) << " "
       << CGAL::to_double(p.y()) << " "
       << CGAL::to_double(p.z()) << " "
       << get(vertex_pmap,vit)
       << std::endl;
  }

  //-------------------------------------------------------
  // Facets
  //-------------------------------------------------------
  os << "Triangles" << std::endl
     << c3t3.number_of_facets() << std::endl;

  for( Facet_iterator fit = c3t3.facets_begin();
       fit != c3t3.facets_end();
       ++fit)
  {
    for (int i=0; i<4; i++)
    {
      if (i != fit->second)
      {
        const Vertex_handle& vh = (*fit).first->vertex(i);
        os << V[vh] << " ";
      }
    }
    os << get(facet_pmap, *fit) //ref
       << std::endl;
  }

  //-------------------------------------------------------
  // Tetrahedra
  //-------------------------------------------------------
  os << "Tetrahedra" << std::endl
     << c3t3.number_of_cells() << std::endl;

  for( Cell_iterator cit = c3t3.cells_begin() ;
       cit != c3t3.cells_end() ;
       ++cit )
  {
    for (int i=0; i<4; i++)
      os << V[cit->vertex(i)] << " ";

    os << get(cell_pmap, cit)
       << std::endl;
  }

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End" << std::endl;

} // end output_to_medit(...)

} // end of anomymous namespace for output_to_medit(..)




namespace Mesh_3 {

namespace details {

  class Debug {
    bool active;
    std::ostream* out;
    const std::string header;

    class Debug_aux {
      Debug* debug;
    public:
      Debug_aux(Debug* d) : debug(d) {}

      template <typename T>
      Debug_aux& operator<<(const T& t)
      {
        debug->apply(t, false);
        return *this;
      }

      operator bool() const
      {
        return false;
      }
    };

    Debug_aux aux;

  public:
    Debug(bool debug,
          std::ostream* debug_str = &std::cout,
          const std::string header_string = "")
      : active(debug), out(debug_str), header(header_string),
        aux(this)
    {
    }

    template <class T>
    void apply(T& t, bool with_header = true)
    {
      if(active)
      {
        if(with_header)
          *out << header;
        *out << t;
      }
    }

    template <typename T>
    Debug_aux& operator<<(const T& t)
    {
      apply(t);
      return aux;
    }

    operator bool() const
    {
      return false;
    }
  }; // end internal class Debug

//   template <typename T>
//   Debug_aux& operator<<(Debug& d, const T& t)
//   {
//     d.apply(t);
//     return aux;
//   }
} // end namespace details
} // end namespace Mesh_3


} // end namespace CGAL

#endif // CGAL_IO_FILE_MEDIT_H
