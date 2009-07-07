// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description : Implements class Mesh_complex_3_in_triangulation_3.
//******************************************************************************

#ifndef MESH_COMPLEX_3_IN_TRIANGULATION_3_H
#define MESH_COMPLEX_3_IN_TRIANGULATION_3_H


#include <CGAL/iterator.h>
#include <CGAL/IO/File_medit.h>
#include <iostream>
#include <fstream>

namespace CGAL {

/**
 * @class Mesh_complex_3_in_triangulation_3
 * @brief A data-structure to represent and maintain a 3D complex embedded
 * in a 3D triangulation.
 */
template<typename Tr>
class Mesh_complex_3_in_triangulation_3
{
  typedef Mesh_complex_3_in_triangulation_3<Tr> Self;

public:
  //-------------------------------------------------------
  // MeshComplex_2InTriangulation3 types
  //-------------------------------------------------------
  // Triangulation types
  typedef Tr                            Triangulation;
  typedef typename Tr::Vertex_handle    Vertex_handle;
  typedef typename Tr::Cell_handle      Cell_handle;
  typedef typename Tr::Facet            Facet;
  typedef typename Tr::Edge             Edge;
  typedef typename Tr::size_type        size_type;

  //-------------------------------------------------------
  // MeshComplex_3InTriangulation3 types
  //-------------------------------------------------------
  // Indices types
  typedef typename Tr::Cell::Subdomain_index  Subdomain_index;
  typedef typename Tr::Cell::Surface_index    Surface_index;
  typedef typename Tr::Vertex::Index          Index;


  //-------------------------------------------------------
  // Constructors / Destructors
  //-------------------------------------------------------
  /**
   * @brief Constructor
   * Builds an empty 3D complex.
   */
  Mesh_complex_3_in_triangulation_3()
    : number_of_facets_(0)
    , tr_()
    , number_of_cells_(0)    { };
  
  /// Copy constructor
  Mesh_complex_3_in_triangulation_3(const Self& rhs)
    : number_of_facets_(rhs.number_of_facets_)
    , tr_(rhs.tr_)
    , number_of_cells_(rhs.number_of_cells_)    { };

  /// Destructor
  ~Mesh_complex_3_in_triangulation_3() { };
  
  /// Assignment operator
  Self& operator=(Self rhs)
  {
    swap(rhs);
    return *this;
  }

  //-------------------------------------------------------
  // MeshComplex_2InTriangulation3 interface
  //-------------------------------------------------------
  /// Returns the reference to the triangulation
  Triangulation& triangulation() { return tr_; }
  /// Returns a const reference to the triangulation
  const Triangulation& triangulation() const { return tr_; }


  /// Adds facet \c facet to the 2D complex, with surface index \c index
  void add_to_complex(const Facet& facet, const Surface_index& index)
  {
    add_to_complex(facet.first, facet.second, index);
  }

  /// Adds facet(\c cell, \c i) to the 2D complex, with surface index \c index
  void add_to_complex(const Cell_handle& cell,
                      const int i,
                      const Surface_index& index);

  /// Removes facet \c facet from 2D complex
  void remove_from_complex(const Facet& facet);

  /// Sets surface index of facet \c facet to \c index
  void set_surface_index(const Facet& f, const Surface_index& index)
  {
    set_surface_index(f.first, f.second, index);
  }

  /// Sets surface index of facet(\c cell, \c i) to \c index
  void set_surface_index(const Cell_handle& cell,
                         const int i,
                         const Surface_index& index)
  {
    cell->set_surface_index(i, index);
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
  Surface_index surface_index(const Facet& f) const
  {
    return surface_index(f.first,f.second);
  }

  /// Returns surface index of facet(\c cell, \c i)
  Surface_index surface_index(const Cell_handle& cell, const int i) const
  {
    return cell->surface_index(i);
  }

  /// Returns the number of surface facets of the mesh
  size_type number_of_facets() const { return number_of_facets_; }

private:
  /// The number of surface facets
  size_type number_of_facets_;

  //-------------------------------------------------------
  // MeshComplex_3InTriangulation3 interface
  //-------------------------------------------------------
public:
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
  size_type number_of_cells() const
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
  void output_to_medit(std::ofstream& os) const
  {
    // Call global function
    CGAL::output_to_medit(os,*this);
  }

  //-------------------------------------------------------
  // Undocumented features
  //-------------------------------------------------------
  /**
   * @brief returns the number of vertices (any dimension) of the mesh
   * @return the number of vertices
   */
  size_type number_of_vertices() const
  {
    return tr_.number_of_vertices();
  }

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

  //-------------------------------------------------------
  // MeshComplex_2InTriangulation3 traversal
  //-------------------------------------------------------
  /**
   * @class Iterator_not_in_complex
   * @brief A class to filter elements which do not belong to the complex
   */
  class Iterator_not_in_complex
  {
    const Self& r_self_;
  public:
    Iterator_not_in_complex(const Self& self) : r_self_(self) { }

    template <typename Iterator>
    bool operator()(Iterator it) const { return ! r_self_.is_in_complex(*it); }
  }; // end class Iterator_not_in_complex


  /// Iterator type to visit the facets of the 2D complex.
  typedef Filter_iterator<typename Triangulation::Finite_facets_iterator,
                          Iterator_not_in_complex> Facet_iterator;

  /// Returns a Facet_iterator to the first facet of the 2D complex
  Facet_iterator facets_begin() const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
                                 Iterator_not_in_complex(*this),
                                 tr_.finite_facets_begin());
  }

  /// Returns past-the-end iterator on facet of the 2D complex
  Facet_iterator facets_end() const
  {
    return CGAL::filter_iterator(tr_.finite_facets_end(),
                                 Iterator_not_in_complex(*this));
  }


  //-------------------------------------------------------
  // MeshComplex_3InTriangulation3 traversal
  //-------------------------------------------------------
  /**
   * @class Cell_not_in_complex
   * @brief A class to filter cells which do not belong to the complex
   */
  class Cell_not_in_complex
  {
    const Self& r_self_;
  public:
    Cell_not_in_complex(const Self& self) : r_self_(self) { }
    bool operator()(Cell_handle ch) const { return !r_self_.is_in_complex(ch); }
  }; // end class Cell_not_in_complex

  /**
   * @class Cell_iterator
   * @brief Iterator type to visit the cells of triangulation belonging
   * to the 3D complex
   *
   * This class is usefull to ensure that Cell_iterator is convertible
   * to Cell_handle
   */
  class Cell_iterator :
    public Filter_iterator<typename Triangulation::Finite_cells_iterator,
                           Cell_not_in_complex>
  {
  private:
    typedef typename Triangulation::Finite_cells_iterator Tr_iterator;
    typedef Filter_iterator<typename Triangulation::Finite_cells_iterator,
                            Cell_not_in_complex> Base;
    typedef Cell_iterator Self;

  public:
    Cell_iterator(Base i) : Base(i) { }

    Self& operator++() { Base::operator++(); return *this; }
    Self& operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }

    operator Cell_handle() { return Cell_handle(this->base()); }
  }; // end class Cell_iterator


  /// Returns a \c Cell_iterator to the first cell of the 3D complex
  Cell_iterator cells_begin() const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this),
                                 tr_.finite_cells_begin());
  }

  /// Returns the past-the-end iterator for the cells of the 3D complex
  Cell_iterator cells_end() const
  {
    return CGAL::filter_iterator(tr_.finite_cells_end(),
                                 Cell_not_in_complex(*this));
  }


private:
  // Private date members
  Triangulation tr_;
  size_type number_of_cells_;

};  // end class Mesh_complex_3_in_triangulation_3


template <typename Tr>
void
Mesh_complex_3_in_triangulation_3<Tr>::add_to_complex(
    const Cell_handle& cell,
    const int i,
    const Surface_index& index)
{
  CGAL_precondition(index != Surface_index());

  if ( ! is_in_complex(cell,i) )
  {
    Facet mirror = tr_.mirror_facet(std::make_pair(cell,i));
    set_surface_index(cell, i, index);
    set_surface_index(mirror.first, mirror.second, index);
    ++number_of_facets_;
  }
}


template <typename Tr>
void
Mesh_complex_3_in_triangulation_3<Tr>::remove_from_complex(const Facet& facet)
{
  if ( is_in_complex(facet) )
  {
    Facet mirror = tr_.mirror_facet(facet);
    set_surface_index(facet.first, facet.second, Surface_index());
    set_surface_index(mirror.first, mirror.second, Surface_index());
    --number_of_facets_;
  }
}

}  // end namespace CGAL

#endif // MESH_COMPLEX_3_IN_TRIANGULATION_3_H
