// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Laurent RINEAU, Stephane Tayeb

#ifndef CGAL_IO_FILE_MEDIT_H
#define CGAL_IO_FILE_MEDIT_H

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <CGAL/utility.h>

namespace CGAL {

namespace Mesh_3 {


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
// Rebin_cell_pmap
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
    typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

    int first_index = 0;
    int index_counter = first_index + 1;

    for( Cell_iterator cell_it = r_c3t3_.cells_in_complex_begin();
         cell_it != r_c3t3_.cells_in_complex_end();
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
    
    // Rebind indices in alphanumeric order
    index_counter = first_index + 1;
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
      return -1;
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
    typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
    std::set<Subdomain_index> subdomain_set;
    
    for( Cell_iterator cell_it = r_c3t3_.cells_in_complex_begin();
        cell_it != r_c3t3_.cells_in_complex_end();
        ++cell_it)
    {
      // Add subdomain index in set if new
      if ( subdomain_set.end() == subdomain_set.find(subdomain_index(cell_it)) )
      {
        subdomain_set.insert(subdomain_index(cell_it));
      }
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
    typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
    
    int first_index = 1;
    int index_counter = first_index;
    
    for( Facet_iterator facet_it = r_c3t3_.facets_in_complex_begin();
        facet_it != r_c3t3_.facets_in_complex_end();
        ++facet_it)
    {
      // Add surface index in internal map if needed
      if ( surface_map_.end() ==
          surface_map_.find(c3t3.surface_patch_index((*facet_it).first,
                                                     (*facet_it).second)) )
      {
        surface_map_.insert(std::make_pair(r_c3t3_.surface_patch_index(*facet_it),
                                           index_counter));
        ++index_counter;
      }
    }
    
    // Find cell_pmap_ unused indices
    typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
    std::set<int> cell_label_set;
    
    for( Cell_iterator cell_it = r_c3t3_.cells_in_complex_begin();
        cell_it != r_c3t3_.cells_in_complex_end();
        ++cell_it)
    {
      // Add subdomain index in set if new
      if ( cell_label_set.end()
          == cell_label_set.find(get(cell_pmap_,cell_it)) )
      {
        cell_label_set.insert(get(cell_pmap_,cell_it));
      }
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
              && facet_index != r_c3t3_.surface_patch_index(*it_facet) )
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
      // should not happen
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
                const C3T3& c3t3)
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
                  Generator().print_twice());
  
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
                const bool print_each_facet_twice = false)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename C3T3::Facets_in_complex_iterator Facet_iterator;
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;

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
       << get(vertex_pmap, vit)
       << std::endl;
  }

  //-------------------------------------------------------
  // Facets
  //-------------------------------------------------------
  typename C3T3::size_type number_of_triangles = c3t3.number_of_facets_in_complex();
  
  if ( print_each_facet_twice )
    number_of_triangles += number_of_triangles;
  
  os << "Triangles" << std::endl
     << number_of_triangles << std::endl;

  for( Facet_iterator fit = c3t3.facets_in_complex_begin();
       fit != c3t3.facets_in_complex_end();
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
    os << get(facet_pmap, *fit) << std::endl;
    
    // Print triangle again if needed
    if ( print_each_facet_twice )
    {
      for (int i=0; i<4; i++)
      {
        if (i != fit->second)
        {
          const Vertex_handle& vh = (*fit).first->vertex(i);
          os << V[vh] << " ";
        }
      }
      os << get(facet_twice_pmap, *fit) << std::endl;
    }
  }

  //-------------------------------------------------------
  // Tetrahedra
  //-------------------------------------------------------
  os << "Tetrahedra" << std::endl
     << c3t3.number_of_cells_in_complex() << std::endl;

  for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit )
  {
    for (int i=0; i<4; i++)
      os << V[cit->vertex(i)] << " ";

    os << get(cell_pmap, cit) << std::endl;
  }

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End" << std::endl;

} // end output_to_medit(...)

} // end namespace Mesh_3

  

  
/**
 * @brief outputs mesh to medit format
 * @param os the stream
 * @param c3t3 the mesh
 * @param rebind if true, labels of cells are rebinded into [1..nb_of_labels]
 * @param show_patches if true, patches are labeled with different labels than
 * cells. If false, each surface facet is written twice, using label of
 * each adjacent cell.
 */
template <class C3T3>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3,
                bool rebind = false,
                bool show_patches = false) 
{
  if ( rebind )
  {
    if ( show_patches )
      Mesh_3::output_to_medit<C3T3,true,false>(os,c3t3);
    else
      Mesh_3::output_to_medit<C3T3,true,true>(os,c3t3);
  }
  else
  {
    if ( show_patches )
      Mesh_3::output_to_medit<C3T3,false,false>(os,c3t3);
    else
      Mesh_3::output_to_medit<C3T3,false,true>(os,c3t3);
  }
}

} // end namespace CGAL

#endif // CGAL_IO_FILE_MEDIT_H
