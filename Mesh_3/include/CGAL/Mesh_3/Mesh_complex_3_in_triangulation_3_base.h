// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description : Implements class Mesh_complex_3_in_triangulation_3.
//******************************************************************************

#ifndef CGAL_MESH_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_BASE_H
#define CGAL_MESH_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_BASE_H

#include <CGAL/Mesh_3/utilities.h>
#include <CGAL/iterator.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/Bbox_3.h>
#include <iostream>
#include <fstream>

namespace CGAL {
namespace Mesh_3 {

/**
 * @class Mesh_complex_3_in_triangulation_3_base
 * @brief A data-structure to represent and maintain a 3D complex embedded
 * in a 3D triangulation.
 */
template<typename Tr>
class Mesh_complex_3_in_triangulation_3_base
{
  typedef Mesh_complex_3_in_triangulation_3_base<Tr> Self;

public:
  // Triangulation types
  typedef Tr                            Triangulation;
  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Cell_handle      Cell_handle;
  typedef typename Tr::Facet            Facet;
  typedef typename Tr::Edge             Edge;
  typedef typename Tr::size_type        size_type;

  // Indices types
  typedef typename Tr::Cell::Subdomain_index      Subdomain_index;
  typedef typename Tr::Cell::Surface_patch_index  Surface_patch_index;
  typedef typename Tr::Vertex::Index              Index;

  //-------------------------------------------------------
  // Constructors / Destructors
  //-------------------------------------------------------
  /**
   * @brief Constructor
   * Builds an empty 3D complex.
   */
  Mesh_complex_3_in_triangulation_3_base()
    : number_of_facets_(0)
    , tr_()
    , number_of_cells_(0)    {}
  
  /// Copy constructor
  Mesh_complex_3_in_triangulation_3_base(const Self& rhs)
    : number_of_facets_(rhs.number_of_facets_)
    , tr_(rhs.tr_)
    , number_of_cells_(rhs.number_of_cells_)    {}

  /// Destructor
  ~Mesh_complex_3_in_triangulation_3_base() {}

  void clear() {
    number_of_cells_ = 0;
    number_of_facets_ = 0;
    tr_.clear();
  }
  
  /// Assignment operator
  Self& operator=(Self rhs)
  {
    swap(rhs);
    return *this;
  }

  /// Returns the reference to the triangulation
  Triangulation& triangulation() { return tr_; }
  /// Returns a const reference to the triangulation
  const Triangulation& triangulation() const { return tr_; }


  /// Adds facet \c facet to the 2D complex, with surface index \c index
  void add_to_complex(const Facet& facet, const Surface_patch_index& index)
  {
    add_to_complex(facet.first, facet.second, index);
  }

  /// Adds facet(\c cell, \c i) to the 2D complex, with surface index \c index
  void add_to_complex(const Cell_handle& cell,
                      const int i,
                      const Surface_patch_index& index);

  /// Removes facet \c facet from 2D complex
  void remove_from_complex(const Facet& facet);

  /// Removes facet(\c cell, \c i) from 2D complex
  void remove_from_complex(const Cell_handle& c, const int i) {
    remove_from_complex(Facet(c, i));
  }

  /// Sets surface index of facet \c facet to \c index
  void set_surface_patch_index(const Facet& f, const Surface_patch_index& index)
  {
    set_surface_patch_index(f.first, f.second, index);
  }

  /// Sets surface index of facet(\c cell, \c i) to \c index
  void set_surface_patch_index(const Cell_handle& cell,
                         const int i,
                         const Surface_patch_index& index)
  {
    cell->set_surface_patch_index(i, index);
  }

  /// Returns true if facet \c facet is in complex
  bool is_in_complex(const Facet& facet) const
  {
    return is_in_complex(facet.first, facet.second);
  }

  /// Returns true if facet (\c cell, \c i) is in 2D complex
  bool is_in_complex(const Cell_handle& cell, const int i) const
  {
    return ( cell->is_facet_on_surface(i) );
  }

  /// Returns surface index of facet \c f
  Surface_patch_index surface_patch_index(const Facet& f) const
  {
    return surface_patch_index(f.first,f.second);
  }

  /// Returns surface index of facet(\c cell, \c i)
  Surface_patch_index surface_patch_index(const Cell_handle& cell,
                                          const int i) const
  {
    return cell->surface_patch_index(i);
  }

  /// Returns the number of surface facets of the mesh
  size_type number_of_facets_in_complex() const { return number_of_facets_; }

  /// Adds cell \c cell to the 3D complex, with subdomain index \c index
  void add_to_complex(const Cell_handle& cell, const Subdomain_index& index)
  {
    CGAL_precondition(index != Subdomain_index());

    if ( ! is_in_complex(cell) )
    {
      set_subdomain_index(cell, index);
      ++number_of_cells_;
    }
  }

  /// Removes cell \c cell from the 3D complex
  void remove_from_complex(const Cell_handle& cell)
  {
    if ( is_in_complex(cell) )
    {
      set_subdomain_index(cell, Subdomain_index());
      --number_of_cells_;
    }
  }


  /// Sets subdomain index of cell \c cell to \c index
  void set_subdomain_index(const Cell_handle& cell,
                           const Subdomain_index& index)
  {
    cell->set_subdomain_index(index);
  }

  /// Sets index of vertex \c vertex to \c index
  void set_index(const Vertex_handle& vertex, const Index& index)
  {
    vertex->set_index(index);
  }
  
  /// Sets dimension of vertex \c vertex to \c dimension
  void set_dimension(const Vertex_handle& vertex, int dimension)
  {
    vertex->set_dimension(dimension);
  }

  /// Returns the number of cells which belongs to the 3D complex
  size_type number_of_cells_in_complex() const
  {
    return number_of_cells_;
  }

  /// Returns \c true if cell \c cell belongs to the 3D complex
  bool is_in_complex(const Cell_handle& cell) const
  {
    return ( subdomain_index(cell) != Subdomain_index() );
  }

  /// Returns the subdomain index of cell \c cell
  Subdomain_index subdomain_index(const Cell_handle& cell) const
  {
    return cell->subdomain_index();
  }

  /// Returns the dimension of the lowest dimensional face of the input 3D
  /// complex that contains the vertex
  int in_dimension(const Vertex_handle& v) const { return v->in_dimension(); }

  /// Returns the index of vertex \c v
  Index index(const Vertex_handle& v) const { return v->index(); }
  
  /// Outputs the mesh to medit
  void output_to_medit(std::ofstream& os,
                       bool rebind = true,
                       bool show_patches = false) const
  {
    // Call global function
    CGAL::output_to_medit(os,*this,rebind,show_patches);
  }

  //-------------------------------------------------------
  // Undocumented features
  //-------------------------------------------------------
  /**
   * @brief insert \c [first,last[ in the triangulation (with dimension 2)
   * @param first the iterator on the first point to insert
   * @param last the iterator past the last point to insert
   *
   * InputIterator value type must be \c std::pair<Tr::Point,Index>
   */
  template <typename InputIterator>
  void insert_surface_points(InputIterator first, InputIterator last)
  {
    while ( first != last )
    {
      Vertex_handle vertex = tr_.insert((*first).first);
      vertex->set_index((*first).second);
      vertex->set_dimension(2);
      ++first;
    }
  }

  /**
   * @brief insert \c [first,last[ in the triangulation (with dimension 2 and
   * index \c default_index)
   * @param first the iterator on the first point to insert
   * @param last the iterator past the last point to insert
   * @param default_index the index to be used to insert points
   *
   * InputIterator value type must be \c Tr::Point
   */
  template <typename InputIterator>
  void insert_surface_points(InputIterator first,
                             InputIterator last,
                             const Index& default_index)
  {
    while ( first != last )
    {
      Vertex_handle vertex = tr_.insert(*first);
      vertex->set_index(default_index);
      vertex->set_dimension(2);
      ++first;
    }
  }
  
  /// Swaps this & rhs
  void swap(Self& rhs)
  {
    std::swap(rhs.number_of_facets_, number_of_facets_);
    tr_.swap(rhs.tr_);
    std::swap(rhs.number_of_cells_, number_of_cells_);
  }
  
  /// Returns bbox
  Bbox_3 bbox() const;
  
  //-------------------------------------------------------
  // Traversal
  //-------------------------------------------------------
private:
  typedef Mesh_3::internal::Iterator_not_in_complex<Self> Iterator_not_in_complex;
  
  class Facet_iterator_not_in_complex
  {
    const Self* c3t3_;
    Surface_patch_index index_; //need by SWIG: should be const Surface_patch_index
  public:
    Facet_iterator_not_in_complex(){} //need by SWIG
    Facet_iterator_not_in_complex(const Self& c3t3,
                                  const Surface_patch_index& index = Surface_patch_index())
      : c3t3_(&c3t3)
      , index_(index) { }
    
    template <typename Iterator>
    bool operator()(Iterator it) const
    { 
      if ( index_ == Surface_patch_index() ) { return ! c3t3_->is_in_complex(*it); }
      else { return c3t3_->surface_patch_index(*it) != index_;  }
    }
  };
  
  /**
   * @class Cell_not_in_complex
   * @brief A class to filter cells which do not belong to the complex
   */
  class Cell_not_in_complex
  {
    const Self* r_self_;
    Subdomain_index index_;//needed by SWIG, should be const Subdomain_index
  public:
    Cell_not_in_complex(){}//needed by SWIG
    Cell_not_in_complex(const Self& self,
                        const Subdomain_index& index = Subdomain_index())
      : r_self_(&self)
      , index_(index) { }

    bool operator()(Cell_handle ch) const
    { 
      if ( index_ == Subdomain_index() ) { return !r_self_->is_in_complex(ch); }
      else { return r_self_->subdomain_index(ch) != index_; }
    }
  }; // end class Cell_not_in_complex
  
public:
  /// Iterator type to visit the facets of the 2D complex.
  typedef Filter_iterator<
    typename Triangulation::Finite_facets_iterator,
    Facet_iterator_not_in_complex >               Facets_in_complex_iterator;

  /// Returns a Facets_in_complex_iterator to the first facet of the 2D complex
  Facets_in_complex_iterator facets_in_complex_begin() const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
                                 Facet_iterator_not_in_complex(*this),
                                 tr_.finite_facets_begin());
  }
  
  /// Returns a Facets_in_complex_iterator to the first facet of the 2D complex
  Facets_in_complex_iterator
  facets_in_complex_begin(const Surface_patch_index& index) const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
                                 Facet_iterator_not_in_complex(*this,index),
                                 tr_.finite_facets_begin());
  }

  /// Returns past-the-end iterator on facet of the 2D complex
  Facets_in_complex_iterator facets_in_complex_end() const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
                                 Facet_iterator_not_in_complex(*this));
  }

  /**
   * @class Cells_in_complex_iterator
   * @brief Iterator type to visit the cells of triangulation belonging
   * to the 3D complex
   *
   * This class is usefull to ensure that Cells_in_complex_iterator is convertible
   * to Cell_handle
   */
  class Cells_in_complex_iterator :
    public Filter_iterator<typename Triangulation::Finite_cells_iterator,
                           Cell_not_in_complex>
  {
  private:
    typedef typename Triangulation::Finite_cells_iterator Tr_iterator;
    typedef Filter_iterator<typename Triangulation::Finite_cells_iterator,
                            Cell_not_in_complex> Base;
    typedef Cells_in_complex_iterator Self;

  public:
    Cells_in_complex_iterator() : Base() { }
    Cells_in_complex_iterator(Base i) : Base(i) { }

    Self& operator++() { Base::operator++(); return *this; }
    Self& operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

    operator Cell_handle() const { return Cell_handle(this->base()); }
  }; // end class Cells_in_complex_iterator


  /// Returns a \c Cells_in_complex_iterator to the first cell of the 3D complex
  Cells_in_complex_iterator cells_in_complex_begin() const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this),
                                 tr_.finite_cells_begin());
  }
  
  /// Returns a \c Cells_in_complex_iterator to the first cell of the 3D complex
  Cells_in_complex_iterator
  cells_in_complex_begin(const Subdomain_index& index) const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this,index),
                                 tr_.finite_cells_begin());
  }

  /// Returns the past-the-end iterator for the cells of the 3D complex
  Cells_in_complex_iterator cells_in_complex_end() const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this));
  }
  
  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;

  void set_surface_index(const Facet& f, const Surface_index& index)
  { set_surface_patch_index(f, index); }
  
  void set_surface_index(const Cell_handle& c, const int i, const Surface_index& index)
  { set_surface_patch_index(c,i,index); }
  
  Surface_index surface_index(const Facet& f) const
  { return surface_patch_index(f); }
  
  Surface_index surface_index(const Cell_handle& c, const int i) const
  { return surface_patch_index(c,i); }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  
#ifndef CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS
  typedef Facets_in_complex_iterator  Facet_iterator;
  typedef Cells_in_complex_iterator   Cell_iterator;
  
  Facet_iterator facets_begin() const
  { return facets_in_complex_begin(); }
  
  Facet_iterator facets_end() const
  { return facets_in_complex_end(); }
  
  Cell_iterator cells_begin() const
  { return cells_in_complex_begin(); }
  
  Cell_iterator cells_end() const
  { return cells_in_complex_end(); }
  
  size_type number_of_facets() const
  { return number_of_facets_in_complex(); }

  size_type number_of_cells() const
  { return number_of_cells_in_complex(); }
#endif // CGAL_MESH_3_NO_DEPRECATED_C3T3_ITERATORS
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------

public:
  template <typename Tr2>
  friend
  std::istream & 
  operator>> (std::istream& is, 
              Mesh_complex_3_in_triangulation_3_base<Tr2> &c3t3);
private:
  // Private date members
  size_type number_of_facets_;
  Triangulation tr_;
  size_type number_of_cells_;

};  // end class Mesh_complex_3_in_triangulation_3_base


template <typename Tr>
void
Mesh_complex_3_in_triangulation_3_base<Tr>::add_to_complex(
    const Cell_handle& cell,
    const int i,
    const Surface_patch_index& index)
{
  CGAL_precondition(index != Surface_patch_index());

  if ( ! is_in_complex(cell,i) )
  {
    Facet mirror = tr_.mirror_facet(std::make_pair(cell,i));
    set_surface_patch_index(cell, i, index);
    set_surface_patch_index(mirror.first, mirror.second, index);
    ++number_of_facets_;
  }
}


template <typename Tr>
void
Mesh_complex_3_in_triangulation_3_base<Tr>::remove_from_complex(const Facet& facet)
{
  if ( is_in_complex(facet) )
  {
    Facet mirror = tr_.mirror_facet(facet);
    set_surface_patch_index(facet.first, facet.second, Surface_patch_index());
    set_surface_patch_index(mirror.first, mirror.second, Surface_patch_index());
    --number_of_facets_;
  }
}
  

// -----------------------------------
// Undocumented
// -----------------------------------
template <typename Tr>
Bbox_3
Mesh_complex_3_in_triangulation_3_base<Tr>::
bbox() const
{
  if ( 0 == triangulation().number_of_vertices() )
  {
    return Bbox_3();
  }
  
  typename Tr::Finite_vertices_iterator vit = tr_.finite_vertices_begin();
  Bbox_3 result = (vit++)->point().bbox();
  
  for(typename Tr::Finite_vertices_iterator end = tr_.finite_vertices_end();
      vit != end ; ++vit)
  {
    result = result + vit->point().bbox();
  }
  
  return result;
}

template < class Tr>
std::ostream & 
operator<< (std::ostream& os, 
            const Mesh_complex_3_in_triangulation_3_base<Tr> &c3t3)
{
  return os << c3t3.triangulation();
}


template < class Tr>
std::istream & 
operator>> (std::istream& is, 
            Mesh_complex_3_in_triangulation_3_base<Tr> &c3t3)
{
  c3t3.clear();
  is >> c3t3.triangulation();

  if(!is) {
    c3t3.clear();
    return is;
  }

  for(typename Tr::Finite_facets_iterator 
        fit = c3t3.triangulation().finite_facets_begin(),
        end = c3t3.triangulation().finite_facets_end();
      fit != end; ++fit) 
  {
    if ( c3t3.is_in_complex(*fit) ) {
      ++c3t3.number_of_facets_;
    }
  }

  for(typename Tr::Finite_cells_iterator 
        cit = c3t3.triangulation().finite_cells_begin(),
        end = c3t3.triangulation().finite_cells_end();
      cit != end; ++cit) 
  {
    if ( c3t3.is_in_complex(cit) ) {
      ++c3t3.number_of_cells_;
    }
  }

  return is;
}

}  // end namespace Mesh_3
}  // end namespace CGAL

#endif // CGAL_MESH_3_MESH_COMPLEX_3_IN_TRIANGULATION_3_BASE_H
