// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb

#ifndef CGAL_IO_FILE_MEDIT_H
#define CGAL_IO_FILE_MEDIT_H

#include <CGAL/license/SMDS_3.h>

#include <CGAL/SMDS_3/Mesh_complex_3_in_triangulation_3_fwd.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>

#include <CGAL/utility.h>
#include <CGAL/basic.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/unordered_map.hpp>

#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include <type_traits>

namespace CGAL {

namespace SMDS_3 {

//-------------------------------------------------------
// Needed in verbose mode
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

// -----------------------------------
// Rebind_cell_pmap
// -----------------------------------
template <typename C3T3>
class Rebind_cell_pmap
{
  typedef typename C3T3::Subdomain_index Subdomain_index;
  typedef std::map<Subdomain_index,int> Subdomain_map;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef unsigned int size_type;

public:
  Rebind_cell_pmap(const C3T3& c3t3)
    : r_c3t3_(c3t3)
  {
    int index_counter = 1;

    for( Cell_handle cell_it : r_c3t3_.cells_in_complex())
    {
      // Add subdomain index in internal map if needed
      std::pair<typename Subdomain_map::iterator, bool> is_insert_successful =
          subdomain_map_.insert(std::make_pair(r_c3t3_.subdomain_index(cell_it),
                                               index_counter));

      if(is_insert_successful.second)
        ++index_counter;
    }

    // Rebind indices in alphanumeric order
    index_counter = 1;
    for ( typename Subdomain_map::iterator mit = subdomain_map_.begin() ;
          mit != subdomain_map_.end() ;
          ++mit )
    {
      mit->second = index_counter++;
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

  size_type subdomain_number() const
  {
    return subdomain_map_.size();
  }

private:
  int subdomain_index(const Subdomain_index& index) const
  {
    typedef typename Subdomain_map::const_iterator Smi;
    Smi elt_it = subdomain_map_.find(index);
    if ( elt_it != subdomain_map_.end() )
      return elt_it->second;
    else
      return 0;
  }

private:
  const C3T3& r_c3t3_;
  Subdomain_map subdomain_map_;
};

// Accessor
template <typename C3T3>
int
get(const Rebind_cell_pmap<C3T3>& cmap,
    const typename C3T3::Cell_handle& ch)
{
  return cmap.subdomain_index(ch);
}

template <typename C3T3>
unsigned int get_size(const Rebind_cell_pmap<C3T3>& cmap)
{
  return cmap.subdomain_number();
}


// -----------------------------------
// No_rebind_cell_pmap
// -----------------------------------
template <typename C3T3>
class No_rebind_cell_pmap
{
  typedef typename C3T3::Subdomain_index Subdomain_index;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef unsigned int size_type;

public:
  No_rebind_cell_pmap(const C3T3& c3t3)
    : r_c3t3_(c3t3) {}

  int subdomain_index(const Cell_handle& ch) const
  {
    return static_cast<int>(r_c3t3_.subdomain_index(ch));
  }

  size_type subdomain_number() const
  {
    std::set<Subdomain_index> subdomain_set;

    for( Cell_handle cell_it : r_c3t3_.cells_in_complex())
    {
      // Add subdomain index in set
      subdomain_set.insert(subdomain_index(cell_it));
    }

    return subdomain_set.size();
  }

private:
  const C3T3& r_c3t3_;
};

// Accessor
template <typename C3T3>
int
get(const No_rebind_cell_pmap<C3T3>& cmap,
    const typename C3T3::Cell_handle& ch)
{
  return cmap.subdomain_index(ch);
}


// -----------------------------------
// Rebind_facet_pmap
// -----------------------------------
template <typename C3T3, typename Cell_pmap>
class Rebind_facet_pmap
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef std::map<Surface_patch_index,int> Surface_map;
  typedef typename C3T3::Facet Facet;
  typedef unsigned int size_type;

public:
  Rebind_facet_pmap(const C3T3& c3t3, const Cell_pmap& cell_pmap)
    : r_c3t3_(c3t3)
    , cell_pmap_(cell_pmap)
  {
    int first_index = 1;
    int index_counter = first_index;

    for( Facet facet_it : r_c3t3_.facets_in_complex())
    {
      // Add surface index in internal map if needed
      std::pair<typename Surface_map::iterator, bool> is_insert_successful =
          surface_map_.insert(std::make_pair(r_c3t3_.surface_patch_index(facet_it),
                                             index_counter));
      if(is_insert_successful.second)
        ++index_counter;
    }

    // Find cell_pmap_ unused indices
    std::set<int> cell_label_set;

    for( typename C3T3::Cell_handle cell_it : r_c3t3_.cells_in_complex())
    {
      // Add subdomain index in set
      cell_label_set.insert(get(cell_pmap_, cell_it));
    }

    // Rebind indices
    index_counter = get_first_unused_label(cell_label_set,first_index);
    for ( typename Surface_map::iterator mit = surface_map_.begin() ;
         mit != surface_map_.end() ;
         ++mit )
    {
      mit->second = index_counter++;
      index_counter = get_first_unused_label(cell_label_set,index_counter);
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
    return surface_index(r_c3t3_.surface_patch_index(f));
  }

  size_type surface_number() const
  {
    return surface_map_.size();
  }

private:
  int surface_index(const Surface_patch_index& index) const
  {
    typedef typename Surface_map::const_iterator Smi;
    Smi elt_it = surface_map_.find(index);
    if ( elt_it != surface_map_.end() )
      return elt_it->second;
    else
      return -1;
  }

  int get_first_unused_label(const std::set<int>& label_set,
                             int search_start) const
  {
    while ( label_set.end() != label_set.find(search_start) )
      ++search_start;

    return search_start;
  }

private:
  const C3T3& r_c3t3_;
  const Cell_pmap& cell_pmap_;
  Surface_map surface_map_;
};


// Accessors
template <typename C3T3, typename Cell_pmap>
int
get(const Rebind_facet_pmap<C3T3,Cell_pmap>& fmap,
    const typename C3T3::Facet& f)
{
  return fmap.surface_index(f);
}

template <typename C3T3, typename Cell_pmap>
unsigned int
get_size(const Rebind_facet_pmap<C3T3,Cell_pmap>& fmap,
         const typename C3T3::Facet& f)
{
  return fmap.surface_number(f);
}


// -----------------------------------
// No_rebind_facet_pmap
// -----------------------------------
template <typename C3T3, typename Cell_pmap>
class No_rebind_facet_pmap
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef typename C3T3::Facet Facet;
  typedef unsigned int size_type;

public:
  No_rebind_facet_pmap(const C3T3& c3t3, const Cell_pmap& /*cell_pmap*/)
    : r_c3t3_(c3t3) {}

  int surface_index(const Facet& f) const
  {
    return static_cast<int>(r_c3t3_.surface_patch_index(f));
  }

private:
  const C3T3& r_c3t3_;
};


// Accessors
template <typename C3T3, typename Cell_pmap>
int
get(const No_rebind_facet_pmap<C3T3,Cell_pmap>& fmap,
    const typename C3T3::Facet& f)
{
return fmap.surface_index(f);
}

// -----------------------------------
// No_rebind_facet_pmap_first
// -----------------------------------
template <typename C3T3, typename Cell_pmap>
class No_rebind_facet_pmap_first
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef typename C3T3::Facet Facet;
  typedef unsigned int size_type;

public:
  No_rebind_facet_pmap_first(const C3T3& c3t3, const Cell_pmap& /*cell_pmap*/)
    : r_c3t3_(c3t3) {}

  int surface_index(const Facet& f) const
  {
    return static_cast<int>(r_c3t3_.surface_patch_index(f).first);
  }

private:
  const C3T3& r_c3t3_;
};


// Accessors
template <typename C3T3, typename Cell_pmap>
int
get(const No_rebind_facet_pmap_first<C3T3,Cell_pmap>& fmap,
  const typename C3T3::Facet& f)
{
  return fmap.surface_index(f);
}


// -----------------------------------
// No_rebind_facet_pmap_second
// -----------------------------------
template <typename C3T3, typename Cell_pmap>
class No_rebind_facet_pmap_second
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef typename C3T3::Facet Facet;
  typedef unsigned int size_type;

public:
  No_rebind_facet_pmap_second(const C3T3& c3t3, const Cell_pmap& /*cell_pmap*/)
  : r_c3t3_(c3t3) {}

  int surface_index(const Facet& f) const
  {
    return static_cast<int>(r_c3t3_.surface_patch_index(f).second);
  }

private:
  const C3T3& r_c3t3_;
};


// Accessors
template <typename C3T3, typename Cell_pmap>
int
get(const No_rebind_facet_pmap_second<C3T3,Cell_pmap>& fmap,
    const typename C3T3::Facet& f)
{
  return fmap.surface_index(f);
}



// -----------------------------------
// No_patch_facet_pmap_first
// -----------------------------------
template <typename C3T3, typename Cell_pmap>
class No_patch_facet_pmap_first
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef typename C3T3::Facet Facet;
  typedef typename C3T3::Cell_handle Cell_handle;

public:
  No_patch_facet_pmap_first(const C3T3&, const Cell_pmap& cell_pmap)
    : cell_pmap_(cell_pmap) { }

  int surface_index(const Facet& f) const
  {
    Cell_handle c1 = f.first;
    Cell_handle c2 = c1->neighbor(f.second);

    int label1 = get(cell_pmap_,c1);
    int label2 = get(cell_pmap_,c2);

    if ( 0 == label1 || -1 == label1 )
      label1 = label2;
    if ( 0 == label2 || -1 == label2 )
      label2 = label1;

    return (std::min)(label1,label2);
  }

private:
  const Cell_pmap& cell_pmap_;
};

// Accessors
template <typename C3T3, typename Cell_pmap>
int
get(const No_patch_facet_pmap_first<C3T3,Cell_pmap>& fmap,
    const typename C3T3::Facet& f)
{
  return fmap.surface_index(f);
}

// -----------------------------------
// No_patch_facet_pmap_second
// -----------------------------------
template <typename C3T3, typename Cell_pmap>
class No_patch_facet_pmap_second
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef typename C3T3::Facet Facet;
  typedef typename C3T3::Cell_handle Cell_handle;

public:
  No_patch_facet_pmap_second(const C3T3&, const Cell_pmap& cell_pmap)
    : cell_pmap_(cell_pmap) { }

  int surface_index(const Facet& f) const
  {
    Cell_handle c1 = f.first;
    Cell_handle c2 = c1->neighbor(f.second);

    int label1 = get(cell_pmap_,c1);
    int label2 = get(cell_pmap_,c2);

    if ( 0 == label1 || -1 == label1 )
      label1 = label2;
    if ( 0 == label2 || -1 == label2 )
      label2 = label1;

    return (std::max)(label1,label2);
  }

private:
  const Cell_pmap& cell_pmap_;
};

// Accessors
template <typename C3T3, typename Cell_pmap>
int
get(const No_patch_facet_pmap_second<C3T3,Cell_pmap>& fmap,
    const typename C3T3::Facet& f)
{
  return fmap.surface_index(f);
}


// -----------------------------------
// Default_vertex_index_pmap
// -----------------------------------
template <typename C3T3, typename Cell_pmap, typename Facet_pmap>
class Default_vertex_pmap
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef typename C3T3::Subdomain_index Subdomain_index;
  typedef typename C3T3::Index Index;
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef typename C3T3::Facet Facet;

public:
  Default_vertex_pmap(const C3T3& c3t3,
                      const Cell_pmap& c_pmap,
                      const Facet_pmap& f_pmap)
    : c_pmap_(c_pmap)
    , f_pmap_(f_pmap)
    , r_c3t3_(c3t3)
    , edge_index_(0) {}

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

        Surface_patch_index facet_index = r_c3t3_.surface_patch_index(*it_facet);
        Facet facet = *it_facet;
        ++it_facet;

        for( ; it_facet != facets.end() ; ++it_facet)
        {
          // If another index is found, return value for edge vertice
          if (   r_c3t3_.is_in_complex(*it_facet)
              && !( facet_index == r_c3t3_.surface_patch_index(*it_facet) ) )
            return edge_index_;
        }

        return get(f_pmap_,facet);
      }
      break;

    case 3:
      {
        // Returns value of any incident cell
        typename std::vector<Cell_handle> cells;
        r_c3t3_.triangulation().finite_incident_cells(
            vh,std::back_inserter(cells));

        if ( cells.begin() != cells.end() )
          return get(c_pmap_, *cells.begin());
        else
          return -1;
      }
      break;

    default:
      // must not happen
      return -1;
      break;
    }
  }

private:
  const Cell_pmap& c_pmap_;
  const Facet_pmap& f_pmap_;
  const C3T3& r_c3t3_;
  const unsigned int edge_index_;
};

template <typename C3T3, typename Cell_pmap, typename Facet_pmap>
int
get(const Default_vertex_pmap<C3T3,Cell_pmap,Facet_pmap>& vmap,
    const typename C3T3::Vertex_handle& vh)
{
  return vmap.index(vh);
}


// -----------------------------------
// Null pmap
// -----------------------------------
template <typename C3T3, typename Cell_pmap>
struct Null_facet_pmap
{
  Null_facet_pmap(const C3T3&, const Cell_pmap&) {}
};

template <typename C3T3, typename Cell_pmap>
int get(const Null_facet_pmap<C3T3,Cell_pmap>&,
        const typename C3T3::Facet&)
{
  return 0;
}

template <typename C3T3, typename Cell_pmap, typename Facet_pmap>
struct Null_vertex_pmap
{
  Null_vertex_pmap(const C3T3&, const Cell_pmap&, const Facet_pmap&) {}
};

template <typename C3T3, typename Cell_pmap, typename Facet_pmap>
int get(const Null_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>&,
        const typename C3T3::Vertex_handle&)
{
  return 0;
}


// -----------------------------------
// Generator
// -----------------------------------
template <typename C3T3, bool rebind, bool no_patch>
struct Medit_pmap_generator{};


template <typename C3T3>
struct Medit_pmap_generator<C3T3, true, false>
{
  typedef Rebind_cell_pmap<C3T3>                            Cell_pmap;
  typedef Rebind_facet_pmap<C3T3, Cell_pmap>                Facet_pmap;
  typedef Null_facet_pmap<C3T3, Cell_pmap>                  Facet_pmap_twice;
  typedef Default_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>  Vertex_pmap;

  bool print_twice() { return false; }
};


template <typename C3T3>
struct Medit_pmap_generator<C3T3, true, true>
{
  typedef Rebind_cell_pmap<C3T3>                            Cell_pmap;
  typedef No_patch_facet_pmap_first<C3T3,Cell_pmap>         Facet_pmap;
  typedef No_patch_facet_pmap_second<C3T3,Cell_pmap>        Facet_pmap_twice;
  typedef Default_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>  Vertex_pmap;

  bool print_twice() { return true; }
};


template <typename C3T3>
struct Medit_pmap_generator<C3T3, false, true>
{
  typedef No_rebind_cell_pmap<C3T3>                         Cell_pmap;
  typedef No_patch_facet_pmap_first<C3T3,Cell_pmap>         Facet_pmap;
  typedef No_patch_facet_pmap_second<C3T3,Cell_pmap>        Facet_pmap_twice;
  typedef Default_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>  Vertex_pmap;

  bool print_twice() { return true; }
};

template <typename C3T3>
struct Medit_pmap_generator<C3T3, false, false>
{
  typedef No_rebind_cell_pmap<C3T3>                         Cell_pmap;
  typedef Rebind_facet_pmap<C3T3,Cell_pmap>                 Facet_pmap;
  typedef Null_facet_pmap<C3T3, Cell_pmap>                  Facet_pmap_twice;
  typedef Null_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>     Vertex_pmap;

  bool print_twice() { return false; }
};


//-------------------------------------------------------
// IO functions
//-------------------------------------------------------



template <class C3T3, bool rebind, bool no_patch>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3,
                const bool all_vertices,
                const bool all_cells)
{
#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cerr << "Output to medit:\n";
#endif

  typedef Medit_pmap_generator<C3T3,rebind,no_patch> Generator;
  typedef typename Generator::Cell_pmap Cell_pmap;
  typedef typename Generator::Facet_pmap Facet_pmap;
  typedef typename Generator::Facet_pmap_twice Facet_pmap_twice;
  typedef typename Generator::Vertex_pmap Vertex_pmap;

  Cell_pmap cell_pmap(c3t3);
  Facet_pmap facet_pmap(c3t3,cell_pmap);
  Facet_pmap_twice facet_pmap_twice(c3t3,cell_pmap);
  Vertex_pmap vertex_pmap(c3t3,cell_pmap,facet_pmap);

  output_to_medit(os,
                  c3t3,
                  vertex_pmap,
                  facet_pmap,
                  cell_pmap,
                  facet_pmap_twice,
                  Generator().print_twice(),
                  all_vertices,
                  all_cells);

#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cerr << "done.\n";
#endif
}



template <class C3T3,
          class Vertex_index_property_map,
          class Facet_index_property_map,
          class Facet_index_property_map_twice,
          class Cell_index_property_map>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3,
                const Vertex_index_property_map& vertex_pmap,
                const Facet_index_property_map& facet_pmap,
                const Cell_index_property_map& cell_pmap,
                const Facet_index_property_map_twice& facet_twice_pmap = Facet_index_property_map_twice(),
                const bool print_each_facet_twice = false,
                const bool all_vertices = true,
                const bool all_cells = false)
{
  typedef typename C3T3::Triangulation Tr;

  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle   Cell_handle;
  typedef typename Tr::Point Point; //can be weighted or not

  const Tr& tr = c3t3.triangulation();

  //-------------------------------------------------------
  // File output
  //-------------------------------------------------------

  //-------------------------------------------------------
  // Header
  //-------------------------------------------------------
  os << std::setprecision(17);

  os << "MeshVersionFormatted 1\n"
     << "Dimension 3\n";
  os << "# CGAL::Mesh_complex_3_in_triangulation_3\n";

  //-------------------------------------------------------
  // Vertices
  //-------------------------------------------------------

  std::unordered_map<Vertex_handle, int> V;
  int inum = 1;
  if (all_vertices || all_cells)
  {
    os << "Vertices\n" << tr.number_of_vertices() << '\n';

    for (typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
         vit != tr.finite_vertices_end();
         ++vit)
    {
      V[vit] = inum++;
      const Point& p = tr.point(vit);
      os << CGAL::to_double(p.x()) << ' '
        << CGAL::to_double(p.y()) << ' '
        << CGAL::to_double(p.z()) << ' '
        << get(vertex_pmap, vit)
        << '\n';
    }
  }
  else
  {
    std::ostringstream oss;
    for (Cell_handle c : c3t3.cells_in_complex())
    {
      for (int i = 0; i < 4; ++i)
      {
        Vertex_handle vit = c->vertex(i);
        if (V.find(vit) == V.end())
        {
          V[vit] = inum++;
          const Point& p = tr.point(vit);
          oss << CGAL::to_double(p.x()) << ' '
            << CGAL::to_double(p.y()) << ' '
            << CGAL::to_double(p.z()) << ' '
            << get(vertex_pmap, vit)
            << '\n';
        }
      }
    }
    os << "Vertices\n" << V.size() << "\n";
    os << oss.str();
  }

  //-------------------------------------------------------
  // Facets
  //-------------------------------------------------------
  typename C3T3::size_type number_of_triangles
    = std::distance(c3t3.facets_in_complex_begin(),
                    c3t3.facets_in_complex_end());

  if ( print_each_facet_twice )
    number_of_triangles += number_of_triangles;

  os << "Triangles\n"
     << number_of_triangles << '\n';

  for(typename C3T3::Facet f : c3t3.facets_in_complex())
  {
    // Apply priority among subdomains, to get consistent facet orientation per subdomain-pair interface.
    if ( print_each_facet_twice )
    {
      // NOTE: We mirror a facet when needed to make it consistent with No_patch_facet_pmap_first/second.
      if (f.first->subdomain_index() > f.first->neighbor(f.second)->subdomain_index())
        f = tr.mirror_facet(f);
    }

    // Get facet vertices in CCW order.
    Vertex_handle vh1 = f.first->vertex((f.second + 1) % 4);
    Vertex_handle vh2 = f.first->vertex((f.second + 2) % 4);
    Vertex_handle vh3 = f.first->vertex((f.second + 3) % 4);

    // Facet orientation also depends on parity.
    if (f.second % 2 != 0)
      std::swap(vh2, vh3);

    os << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3] << ' ';
    os << get(facet_pmap, f) << '\n';

    // Print triangle again if needed, with opposite orientation
    if ( print_each_facet_twice )
    {
      os << V[vh3] << ' ' << V[vh2] << ' ' << V[vh1] << ' ';
      os << get(facet_twice_pmap, f) << '\n';
    }
  }

  //-------------------------------------------------------
  // Tetrahedra
  //-------------------------------------------------------
  typename C3T3::size_type number_of_cells
    = all_cells
    ? c3t3.triangulation().number_of_finite_cells()
    : std::distance(c3t3.cells_in_complex_begin(), c3t3.cells_in_complex_end());;
  os << "Tetrahedra\n"
     << number_of_cells << '\n';

  if (all_cells)
  {
    for (auto cit = c3t3.triangulation().finite_cells_begin();
         cit != c3t3.triangulation().finite_cells_end();
         ++cit)
    {
      for (int i = 0; i < 4; i++)
        os << V[cit->vertex(i)] << ' ';

      os << get(cell_pmap, cit) << '\n';
    }
  }
  else
  {
    for (Cell_handle cit : c3t3.cells_in_complex())
    {
      for (int i = 0; i < 4; i++)
        os << V[cit->vertex(i)] << ' ';

      os << get(cell_pmap, cit) << '\n';
    }
  }

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End\n";

} // end output_to_medit(...)

} // end namespace Mesh_3

namespace IO {

/**
 * @ingroup PkgSMDS3IOFunctions
 * @deprecated This function is deprecated. Users should instead use `CGAL::IO::write_MEDIT()`
 * @brief outputs a mesh complex to the medit (`.mesh`) file format.
        See \cgalCite{frey:inria-00069921} for a comprehensive description of this file format.
 * @param os the output stream
 * @param c3t3 the mesh complex
 * @param rebind if `true`, labels of cells are rebinded into `[1..nb_of_labels]`
 * @param show_patches if `true`, patches are labeled with different labels than
 *                     cells. If `false`, each surface facet is written twice,
 *                     using the label of each adjacent cell.
 * \see \ref IOStreamMedit
 */
template <class C3T3>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3,
                bool rebind,      // = false,
                bool show_patches // = false
#ifndef DOXYGEN_RUNNING
              , bool all_vertices // = true
              , bool all_cells    // = false
#endif
)
{
  if ( rebind )
  {
    if ( show_patches )
      CGAL::SMDS_3::output_to_medit<C3T3,true,false>(os, c3t3,
        all_vertices, all_cells);
    else
      CGAL::SMDS_3::output_to_medit<C3T3,true,true>(os, c3t3,
        all_vertices, all_cells);
  }
  else
  {
    if ( show_patches )
      CGAL::SMDS_3::output_to_medit<C3T3,false,false>(os, c3t3,
        all_vertices, all_cells);
    else
      CGAL::SMDS_3::output_to_medit<C3T3,false,true>(os, c3t3,
        all_vertices, all_cells);
  }
}

/**
 * @ingroup PkgSMDS3IOFunctions
 * @brief outputs a mesh complex to the medit (`.mesh`) file format.
 *      See \cgalCite{frey:inria-00069921} for a comprehensive description of this file format.
 * @tparam T3 can be instantiated with any 3D triangulation of \cgal provided that its
 *  vertex and cell base class are models of the concepts `SimplicialMeshVertexBase_3`
 * and `SimplicialMeshCellBase_3`, respectively.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param os the output stream
 * @param t3 the triangulation
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 * \cgalParamNBegin{all_cells}
 *   \cgalParamDescription{If `true`, all the cells in `t3` are written in `os`,
 *        whether they belong to the complex or not.
 *        Otherwise, only the cells `c` for which
 *        `c->subdomain_index() != Subdomain_index()` are written.}
 *   \cgalParamType{Boolean}
 *   \cgalParamDefault{`true`}
 *   \cgalParamExtra{This parameter must be set to `true` for the file to be readable by `read_MEDIT()`.}
 * \cgalParamNEnd
 *
 * \cgalParamNBegin{all_vertices}
 *   \cgalParamDescription{If `true`, all the finite vertices in `t3` are written in `os`.
 *                         Otherwise, only the vertices that belong to a cell `c` for which
 *                         `c->subdomain_index() != Subdomain_index()` are written}
 *   \cgalParamType{Boolean}
 *   \cgalParamDefault{`true`}
 *   \cgalParamExtra{If `all_cells` is `true`, the value of this parameter is ignored and
                     all vertices are written in `os`. It must be
 *                   set to `true` for the file to be readable by `read_MEDIT()`.}
 * \cgalParamNEnd
 *
 * \cgalParamNBegin{rebind_labels}
 *  \cgalParamDescription{If `true`, labels of cells are rebinded into `[1..nb_of_labels]`}
 *  \cgalParamType{Boolean}
 *  \cgalParamDefault{`false`}
 * \cgalParamNEnd
 *
 * \cgalParamNBegin{show_patches}
 *  \cgalParamDescription{If `true`, patches are labeled with different labels than
 *                        cells. If `false`, each surface facet is written twice,
 *                        using the label of each adjacent cell.}
 *  \cgalParamType{Boolean}
 *  \cgalParamDefault{`true`}
 * \cgalParamNEnd
 * \cgalNamedParamsEnd
 * \see \ref IOStreamMedit
 */
template<typename T3, typename NamedParameters = parameters::Default_named_parameters>
void write_MEDIT(std::ostream& os,
                 const T3& t3,
                 const NamedParameters& np = parameters::default_values())
{
  CGAL::Mesh_complex_3_in_triangulation_3<T3, int, int> c3t3;
  c3t3.triangulation() = t3;
  c3t3.rescan_after_load_of_triangulation();

  using parameters::get_parameter;
  using parameters::choose_parameter;

  bool rebind = choose_parameter(get_parameter(np, internal_np::rebind_labels), false);;
  bool show_patches = choose_parameter(get_parameter(np, internal_np::show_patches), true);
  bool all_c = choose_parameter(get_parameter(np, internal_np::all_cells), true);
  bool all_v = all_c || choose_parameter(get_parameter(np, internal_np::all_vertices), true);

  output_to_medit(os, c3t3, rebind, show_patches, all_v, all_c);
}

/**
 * @ingroup PkgSMDS3IOFunctions
 * @brief outputs a mesh complex to the medit (`.mesh`) file format.
 *      See \cgalCite{frey:inria-00069921} for a comprehensive description of this file format.
 * @tparam T3 can be instantiated with any 3D triangulation of \cgal provided that its
 *  vertex and cell base class are models of the concepts `MeshVertexBase_3` and `MeshCellBase_3`, respectively.
 * @tparam CornerIndex is the type of the indices for corners
 * @tparam CurveIndex is the type of the indices for curves
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

 * @param os the output stream
 * @param c3t3 the mesh complex
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 * \cgalParamNBegin{all_cells}
 *   \cgalParamDescription{If `true`, all the cells in `t3` are written in `os`,
 *        whether they belong to the complex or not.
 *        Otherwise, only the cells `c` for which
 *        `c->subdomain_index() != Subdomain_index()` are written.}
 *   \cgalParamType{Boolean}
 *   \cgalParamDefault{`true`}
 *   \cgalParamExtra{If the complex does not form a topological sphere,
 *        this parameter must be set to `true` for the file to be readable by `read_MEDIT()`.
 *        Otherwise the underlying triangulation data structure will not be valid.}
 * \cgalParamNEnd
 *
 * \cgalParamNBegin{all_vertices}
 *   \cgalParamDescription{If `true`, all the vertices in `t3` are written in `os`.
 *                         Otherwise, only the vertices that belong to a cell `c` for which
 *                         `c->subdomain_index() != Subdomain_index()` are written}
 *   \cgalParamType{Boolean}
 *   \cgalParamDefault{`true`}
 *   \cgalParamExtra{If `all_cells` is `true`, the value of this parameter is ignored
          and all vertices are written in `os`. If the complex does not
          form a topological sphere, it must be
 *        set to `true` for the file to be readable by `read_MEDIT()`.
 *        Otherwise the underlying triangulation data structure will not be valid.}
 * \cgalParamNEnd
 *
 * \cgalParamNBegin{rebind_labels}
 *  \cgalParamDescription{If `true`, labels of cells are rebinded into `[1..nb_of_labels]`}
 *  \cgalParamType{Boolean}
 *  \cgalParamDefault{`false`}
 * \cgalParamNEnd
 *
 * \cgalParamNBegin{show_patches}
 *  \cgalParamDescription{If `true`, patches are labeled with different labels than
 *                        cells. If `false`, each surface facet is written twice,
 *                        using the label of each adjacent cell.}
 *  \cgalParamType{Boolean}
 *  \cgalParamDefault{`true`}
 * \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \see \ref IOStreamMedit
 */
template<typename T3,
         typename CornerIndex,
         typename CurveIndex,
         typename NamedParameters = parameters::Default_named_parameters>
void write_MEDIT(std::ostream& os,
  const CGAL::Mesh_complex_3_in_triangulation_3<T3, CornerIndex, CurveIndex>& c3t3,
  const NamedParameters& np = parameters::default_values())
{
  return write_MEDIT(os, c3t3.triangulation(), np);
}

/**
 * @ingroup PkgSMDS3IOFunctions
 * @brief reads a mesh complex written in the medit (`.mesh`) file format.
 *   See \cgalCite{frey:inria-00069921} for a comprehensive description of this file format.
 * @tparam T3 can be instantiated with any 3D triangulation of \cgal provided that its
 *  vertex and cell base class are models of the concepts `MeshVertexBase_3` and `MeshCellBase_3`,
 *  respectively.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param in the input stream
 * @param t3 the triangulation
 * @param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{indicates whether output warnings and error messages should be printed or not.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{allow_non_manifold}
 *     \cgalParamDescription{allows the construction of a triangulation with non-manifold edges
 *       and non manifold vertices. The triangulation is invalid if this situation is met,
 *       so it should be used only in advanced cases, and the triangulation will be hardly usable.}
 *     \cgalParamType{bool}
 *     \cgalParamDefault{false}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @returns `true` if the connectivity of the triangulation could be built consistently
 * from \p in,
 * and `false` if the triangulation is empty, or if the connectivity
 * of \p t3 could not be built.
 * If `false` is returned, \p t3 is empty when the function returns.
 *
 * This function reads the data about vertices, surface facets, and
 * triangulation cells from `in`, and builds a valid `T3` from it.
 *
 * Note that a valid 3D triangulation of \cgal must have a valid
 * data structure (see `TriangulationDataStructure_3 `),
 * positively oriented cells,
 * and cover the geometric convex hull of all points in `t3`.
 */
template<typename T3, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_MEDIT(std::istream& in,
                T3& t3,
                const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // Default non_manifold value is true if the triangulation periodic, false otherwise
  const bool non_manifold = choose_parameter(get_parameter(np, internal_np::allow_non_manifold),
                                             std::is_same<typename T3::Periodic_tag, Tag_true>::value);
  const bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);

  bool b = CGAL::SMDS_3::build_triangulation_from_file(in, t3, verbose, false /*replace_domain_0*/, non_manifold);
  if(!b)
    t3.clear();
  return b;
}

} // namespace IO

#ifndef CGAL_NO_DEPRECATED_CODE
using IO::output_to_medit;
#endif

} // end namespace CGAL

#endif // CGAL_IO_FILE_MEDIT_H
