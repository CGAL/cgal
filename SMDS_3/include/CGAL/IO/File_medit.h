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
// Renumber_subdomains_pmap
// -----------------------------------
template <typename C3T3>
class Renumber_subdomains_pmap
{
  typedef typename C3T3::Subdomain_index Subdomain_index;
  typedef std::map<Subdomain_index,int> Subdomain_map;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef unsigned int size_type;

public:
  Renumber_subdomains_pmap(const C3T3& c3t3)
    : r_c3t3_(c3t3)
  {
    int index_counter = 1;

    for( Cell_handle c : r_c3t3_.cells_in_complex())
    {
      // Add subdomain index in internal map if needed
      auto [_, is_insert_successful] =
          subdomain_map_.insert(std::make_pair(r_c3t3_.subdomain_index(c),
                                               index_counter));

      if(is_insert_successful)
        ++index_counter;
    }

    // Renumber indices in alphanumeric order
    index_counter = 1;
    for ( auto& [_, index_ref] : subdomain_map_)
    {
      index_ref = index_counter++;
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

  friend int get(const Renumber_subdomains_pmap& cmap, const Cell_handle& ch)
  {
    return cmap.subdomain_index(ch);
  }

  friend unsigned int get_size(const Renumber_subdomains_pmap& cmap)
  {
    return cmap.subdomain_number();
  }

private:
  int subdomain_index(const Subdomain_index& index) const
  {
    auto elt_it = subdomain_map_.find(index);
    if ( elt_it != subdomain_map_.end() )
      return elt_it->second;
    else
      return 0;
  }

private:
  const C3T3& r_c3t3_;
  Subdomain_map subdomain_map_;
};


// -----------------------------------
// Use_subdomain_indices
// -----------------------------------
template <typename C3T3>
class Use_subdomain_indices
{
  typedef typename C3T3::Subdomain_index Subdomain_index;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef unsigned int size_type;

public:
  Use_subdomain_indices(const C3T3& c3t3)
    : r_c3t3_(c3t3) {}

  int subdomain_index(const Cell_handle& ch) const
  {
    return static_cast<int>(r_c3t3_.subdomain_index(ch));
  }

  size_type subdomain_number() const
  {
    std::set<Subdomain_index> subdomain_set;

    for( Cell_handle c : r_c3t3_.cells_in_complex())
    {
      // Add subdomain index in set
      subdomain_set.insert(subdomain_index(c));
    }

    return subdomain_set.size();
  }

  friend int get(const Use_subdomain_indices& cmap, const Cell_handle& ch)
  {
    return cmap.subdomain_index(ch);
  }

private:
  const C3T3& r_c3t3_;
};


// -----------------------------------
// Renumber_surface_patches_pmap
// -----------------------------------
template <typename C3T3, typename Cell_pmap>
class Renumber_surface_patches_pmap
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef std::map<Surface_patch_index,int> Surface_map;
  typedef typename C3T3::Facet Facet;
  typedef unsigned int size_type;

public:
  Renumber_surface_patches_pmap(const C3T3& c3t3, const Cell_pmap& cell_pmap)
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

    for( auto c : r_c3t3_.cells_in_complex())
    {
      // Add subdomain index in set
      cell_label_set.insert(get(cell_pmap_, c));
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

  friend int get(const Renumber_surface_patches_pmap& fmap, const Facet& f)
  {
    return fmap.surface_index(f);
  }

  friend unsigned int get_size(const Renumber_surface_patches_pmap& fmap)
  {
    return fmap.surface_number();
  }

private:
  const C3T3& r_c3t3_;
  const Cell_pmap& cell_pmap_;
  Surface_map surface_map_;
};


// -----------------------------------
// Use_cell_indices_pmap
// -----------------------------------
template <typename C3T3, typename Cell_pmap, int zero_or_one>
class Use_cell_indices_pmap
{
  using Facet = typename C3T3::Facet;

public:
  Use_cell_indices_pmap(const C3T3&, const Cell_pmap& cell_pmap)
    : cell_pmap_(cell_pmap) { }

  int surface_index(const Facet& f) const
  {
    auto c1 = f.first;
    auto c2 = c1->neighbor(f.second);

    int label1 = get(cell_pmap_, c1);
    int label2 = get(cell_pmap_, c2);

    if ( 0 == label1 || -1 == label1 )
      label1 = label2;
    if ( 0 == label2 || -1 == label2 )
      label2 = label1;

    return std::get<zero_or_one>(std::minmax(label1,label2));
  }

  friend auto get(const Use_cell_indices_pmap& fmap, const Facet& f)
  {
    return fmap.surface_index(f);
  }

private:
  const Cell_pmap& cell_pmap_;
};


// -----------------------------------
// Default_vertex_pmap
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

  friend int get(const Default_vertex_pmap& vmap, const Vertex_handle& vh)
  {
    return vmap.index(vh);
  }

private:
  const Cell_pmap& c_pmap_;
  const Facet_pmap& f_pmap_;
  const C3T3& r_c3t3_;
  const unsigned int edge_index_;
};


// -----------------------------------
// Null pmap
// -----------------------------------
struct Null_pmap
{
  template <typename ...Args> Null_pmap(Args&&...) {}

  template <typename T>
  friend int get(const Null_pmap&, const T&)
  {
    return 0;
  }
};


// -----------------------------------
// Generator
// -----------------------------------

enum Renumber_subdomain_indices : bool { RENUMBER_SUBDOMAINS = true, USE_SUBDOMAIN_INDICES = false };
enum Facet_indices : bool { USE_CELL_INDICES = true, RENUMBER_SURFACE_PATCH_INDICES = false };

template <typename, Renumber_subdomain_indices, Facet_indices>
struct Medit_pmap_generator;

template <typename C3T3>
struct Medit_pmap_generator<C3T3, RENUMBER_SUBDOMAINS, RENUMBER_SURFACE_PATCH_INDICES>
{
  typedef Renumber_subdomains_pmap<C3T3>                    Cell_pmap;
  typedef Renumber_surface_patches_pmap<C3T3, Cell_pmap>    Facet_pmap;
  typedef Null_pmap                                         Facet_pmap_twice;
  typedef Default_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>  Vertex_pmap;

  bool print_twice() { return false; }
};


template <typename C3T3>
struct Medit_pmap_generator<C3T3, RENUMBER_SUBDOMAINS, USE_CELL_INDICES>
{
  typedef Renumber_subdomains_pmap<C3T3>                    Cell_pmap;
  typedef Use_cell_indices_pmap<C3T3, Cell_pmap, 0>         Facet_pmap;
  typedef Use_cell_indices_pmap<C3T3, Cell_pmap, 1>         Facet_pmap_twice;
  typedef Default_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>  Vertex_pmap;

  bool print_twice() { return true; }
};


template <typename C3T3>
struct Medit_pmap_generator<C3T3, USE_SUBDOMAIN_INDICES, USE_CELL_INDICES>
{
  typedef Use_subdomain_indices<C3T3>                       Cell_pmap;
  typedef Use_cell_indices_pmap<C3T3, Cell_pmap, 0>         Facet_pmap;
  typedef Use_cell_indices_pmap<C3T3, Cell_pmap, 1>         Facet_pmap_twice;
  typedef Default_vertex_pmap<C3T3, Cell_pmap, Facet_pmap>  Vertex_pmap;

  bool print_twice() { return true; }
};

template <typename C3T3>
struct Medit_pmap_generator<C3T3, USE_SUBDOMAIN_INDICES, RENUMBER_SURFACE_PATCH_INDICES>
{
  typedef Use_subdomain_indices<C3T3>                       Cell_pmap;
  typedef Renumber_surface_patches_pmap<C3T3, Cell_pmap>    Facet_pmap;
  typedef Null_pmap                                         Facet_pmap_twice;
  typedef Null_pmap                                         Vertex_pmap;

  bool print_twice() { return false; }
};


//-------------------------------------------------------
// IO functions
//-------------------------------------------------------



template <class Tr,
          class Vertices_range,
          class Facets_range,
          class Cells_range,
          class Vertex_index_property_map,
          class Facet_index_property_map,
          class Facet_index_property_map_twice = Null_pmap,
          class Cell_index_property_map>
void
output_to_medit(std::ostream& os,
                const Tr& tr,
                const Vertices_range& vertices,
                const Facets_range& facets,
                const Cells_range& cells,
                const Vertex_index_property_map& vertex_pmap,
                const Facet_index_property_map& facet_pmap,
                const Cell_index_property_map& cell_pmap,
                const Facet_index_property_map_twice& facet_twice_pmap = {},
                const bool print_each_facet_twice = false)
{
  using std::size;
  using Vertex_handle = typename Tr::Vertex_handle;

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

  std::ostringstream oss;
  oss.precision(os.precision());
  for(auto v: vertices) {
    auto& v_num = V[v];
    if(v_num != 0) return;
    v_num = inum++;
    const auto& p = tr.point(v);
    oss << CGAL::to_double(p.x()) << ' '
        << CGAL::to_double(p.y()) << ' '
        << CGAL::to_double(p.z()) << ' '
        << get(vertex_pmap, v)
        << '\n';
  }
  os << "Vertices\n" << V.size() << "\n";
  os << oss.str();

  //-------------------------------------------------------
  // Facets
  //-------------------------------------------------------
  auto number_of_triangles = size(facets);

  if ( print_each_facet_twice )
    number_of_triangles += number_of_triangles;

  os << "Triangles\n"
     << number_of_triangles << '\n';

  for (auto f : facets) {
    auto [c, index] = f;
    // Apply priority among subdomains, to get consistent facet orientation per subdomain-pair interface.
    if (print_each_facet_twice) {
      auto mirror_facet = tr.mirror_facet(f);
      [[maybe_unused]] auto [c2, _] = mirror_facet;
      // NOTE: We mirror a facet when needed to make it consistent with Use_cell_indices_pmap.
      if (get(cell_pmap, c) > get(cell_pmap, c2)) {
        f = mirror_facet;
      }
    }

    // Get facet vertices in CCW order.
    auto [vh1, vh2, vh3] = tr.vertices(f);

    os << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3] << ' ';
    os << get(facet_pmap, f) << '\n';

    // Print triangle again if needed, with opposite orientation
    if (print_each_facet_twice) {
      os << V[vh3] << ' ' << V[vh2] << ' ' << V[vh1] << ' ';
      os << get(facet_twice_pmap, f) << '\n';
    }
  }

  //-------------------------------------------------------
  // Tetrahedra
  //-------------------------------------------------------
  os << "Tetrahedra\n"
     << size(cells) << '\n';
  for (const auto& c : cells) {
    for (auto v : tr.vertices(c))
      os << V[v] << ' ';
    os << get(cell_pmap, c) << '\n';
  }

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End\n";

} // end output_to_medit(...)


template <class Tr,
          class Vertices_range,
          class Facets_range,
          class Cells_range>
void
output_T3_to_medit(std::ostream& os,
                const Tr& tr,
                const Vertices_range& vertices,
                const Facets_range& facets,
                const Cells_range& cells)
{
  using std::size;
  using Vertex_handle = typename Tr::Vertex_handle;

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

  std::ostringstream oss;
  oss.precision(os.precision());
  for(auto v: vertices) {
    auto& v_num = V[v];
    if(v_num != 0) return;
    v_num = inum++;
    const auto& p = tr.point(v);
    oss << CGAL::to_double(p.x()) << ' '
        << CGAL::to_double(p.y()) << ' '
        << CGAL::to_double(p.z()) << ' '
        << 0
        << '\n';
  }
  os << "Vertices\n" << V.size() << "\n";
  os << oss.str();


  #if 1
  //-------------------------------------------------------
  // Facets
  //-------------------------------------------------------
  auto number_of_triangles = size(facets);
 bool print_each_facet_twice = false;

  if ( print_each_facet_twice )
    number_of_triangles += number_of_triangles;

  os << "Triangles\n"
     << number_of_triangles << '\n';

  for (auto f : facets) {
    // auto [c, index] = f;
    // Apply priority among subdomains, to get consistent facet orientation per subdomain-pair interface.
    if (print_each_facet_twice) {
      auto mirror_facet = tr.mirror_facet(f);
      [[maybe_unused]] auto [c2, _] = mirror_facet;
      // NOTE: We mirror a facet when needed to make it consistent with Use_cell_indices_pmap.
      if (true) /* AF ???? (get(cell_pmap, c) > get(cell_pmap, c2))*/ {
        f = mirror_facet;
      }
    }

    // Get facet vertices in CCW order.
    auto [vh1, vh2, vh3] = tr.vertices(f);

    os << V[vh1] << ' ' << V[vh2] << ' ' << V[vh3] << ' ';
    os << 1 << '\n';

    // Print triangle again if needed, with opposite orientation
    if (print_each_facet_twice) {
      os << V[vh3] << ' ' << V[vh2] << ' ' << V[vh1] << ' ';
      os << 1 << '\n';
    }
  }
#endif
  //-------------------------------------------------------
  // Tetrahedra
  //-------------------------------------------------------
  os << "Tetrahedra\n"
     << size(cells) << '\n';
  for (const auto& c : cells) {
    for (auto v : tr.vertices(c))
      os << V[v] << ' ';
    os << 1 << '\n';
  }

  //-------------------------------------------------------
  // End
  //-------------------------------------------------------
  os << "End\n";

} // end output_T3_to_medit(...)




template <class C3T3, Renumber_subdomain_indices renumber_subdomain_indices, Facet_indices no_patch>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3,
                const bool all_vertices,
                const bool all_cells)
{
#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cerr << "Output to medit:\n";
#endif

  typedef Medit_pmap_generator<C3T3, renumber_subdomain_indices, no_patch> Generator;
  typedef typename Generator::Cell_pmap Cell_pmap;
  typedef typename Generator::Facet_pmap Facet_pmap;
  typedef typename Generator::Facet_pmap_twice Facet_pmap_twice;
  typedef typename Generator::Vertex_pmap Vertex_pmap;

  Cell_pmap cell_pmap(c3t3);
  Facet_pmap facet_pmap(c3t3,cell_pmap);
  Facet_pmap_twice facet_pmap_twice(c3t3,cell_pmap);
  Vertex_pmap vertex_pmap(c3t3,cell_pmap,facet_pmap);

  const auto& tr = c3t3.triangulation();

  auto all_vertices_range = tr.finite_vertex_handles();
  auto all_cells_range = tr.finite_cell_handles();
  auto cells_in_complex_range = c3t3.cells_in_complex();

  auto output_to_medit = [&](const auto& vertices, const auto& cells) {
    CGAL::SMDS_3::output_to_medit(os, tr,
                                  vertices,
                                  c3t3.facets_in_complex(),
                                  cells,
                                  vertex_pmap, facet_pmap, cell_pmap, facet_pmap_twice,
                                  Generator().print_twice());
  };

  if(false == all_vertices && false == all_cells) {
    std::set<typename C3T3::Vertex_handle> vertices;
    for(auto c : cells_in_complex_range) {
      for(auto v: tr.vertices(c)) {
        vertices.insert(v);
      }
    }
    output_to_medit(vertices, cells_in_complex_range);
  } else {
    if(all_cells) {
      output_to_medit(all_vertices_range, all_cells_range);
    } else {
      // here, necessarily `all_vertices == true`
      output_to_medit(all_vertices_range, cells_in_complex_range);
    }
  }

#ifdef CGAL_MESH_3_IO_VERBOSE
  std::cerr << "done.\n";
#endif
}


} // end namespace Mesh_3

namespace IO {

/**
 * @ingroup PkgSMDS3IOFunctions
 * @deprecated This function is deprecated. Users should instead use `CGAL::IO::write_MEDIT()`
 * @brief outputs a mesh complex to the medit (`.mesh`) file format.
        See \cgalCite{frey:inria-00069921} for a comprehensive description of this file format.
 * @param os the output stream
 * @param c3t3 the mesh complex
 * @param renumber_subdomain_indices if `true`, labels of cells are renumbered into `[1..nb_of_labels]`
 * @param show_patches if `true`, patches are labeled with different labels than
 *                     cells. If `false`, each surface facet is written twice,
 *                     using the label of each adjacent cell.
 * \see \ref IOStreamMedit
 */
template <class C3T3>
void
output_to_medit(std::ostream& os,
                const C3T3& c3t3,
                bool renumber_subdomain_indices, // = false,
                bool show_patches // = false
#ifndef DOXYGEN_RUNNING
              , bool all_vertices // = true
              , bool all_cells    // = false
#endif
)
{
  using namespace CGAL::SMDS_3;
  if ( renumber_subdomain_indices )
  {
    if ( show_patches )
      CGAL::SMDS_3::output_to_medit<C3T3,RENUMBER_SUBDOMAINS,RENUMBER_SURFACE_PATCH_INDICES>(
          os, c3t3, all_vertices, all_cells);
    else
      CGAL::SMDS_3::output_to_medit<C3T3,RENUMBER_SUBDOMAINS,USE_CELL_INDICES>(os, c3t3,
        all_vertices, all_cells);
  }
  else
  {
    if ( show_patches )
      CGAL::SMDS_3::output_to_medit<C3T3,USE_SUBDOMAIN_INDICES,RENUMBER_SURFACE_PATCH_INDICES>(
          os, c3t3, all_vertices, all_cells);
    else
      CGAL::SMDS_3::output_to_medit<C3T3,USE_SUBDOMAIN_INDICES,USE_CELL_INDICES>(os, c3t3,
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
 *  \cgalParamDescription{If `true`, labels of cells are rebound into `[1..nb_of_labels]`}
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

  bool renumber_subdomain_indices = choose_parameter(get_parameter(np, internal_np::rebind_labels), false);;
  bool show_patches = choose_parameter(get_parameter(np, internal_np::show_patches), true);
  bool all_c = choose_parameter(get_parameter(np, internal_np::all_cells), true);
  bool all_v = all_c || choose_parameter(get_parameter(np, internal_np::all_vertices), true);

  output_to_medit(os, c3t3, renumber_subdomain_indices, show_patches, all_v, all_c);
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
 *  \cgalParamDescription{If `true`, labels of cells are rebound into `[1..nb_of_labels]`}
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
