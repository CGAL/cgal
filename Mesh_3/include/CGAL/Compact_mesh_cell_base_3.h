// Copyright (c) 2006-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008,2011 GeometryFactory Sarl (France)
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
// $URL: https://scm.gforge.inria.fr/svn/cgal/branches/features/Mesh_3-experimental-GF/Mesh_3/include/CGAL/Compact_mesh_cell_base_3.h $
// $Id: Compact_mesh_cell_base_3.h 70288 2012-07-05 10:09:48Z jtournoi $
//
//
// Author(s)     : Laurent Rineau, Stephane Tayeb, Andreas Fabri


#ifndef CGAL_COMPACT_MESH_CELL_BASE_3_H
#define CGAL_COMPACT_MESH_CELL_BASE_3_H

#include <CGAL/Mesh_3/config.h>

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/internal/Dummy_tds_3.h>

#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Mesh_3/io_signature.h>

namespace CGAL {
  
// Class Compact_mesh_cell_base_3
// Cell base class used in 3D meshing process.
// Adds information to Cb about the cell of the input complex containing it
template< class GT,
          class MD,
          class TDS = void > 
class Compact_mesh_cell_base_3
{
  typedef typename GT::FT FT;
  
public:
  typedef TDS                          Triangulation_data_structure;
  typedef typename TDS::Vertex_handle  Vertex_handle;
  typedef typename TDS::Cell_handle    Cell_handle;
  typedef typename TDS::Vertex         Vertex;
  typedef typename TDS::Cell           Cell;
  typedef typename TDS::Cell_data      TDS_data;


  template <typename TDS2>
  struct Rebind_TDS {
    typedef Compact_mesh_cell_base_3<GT, MD, TDS2> Other;
  };


  // Index Type
  typedef typename MD::Subdomain_index      Subdomain_index;
  typedef typename MD::Surface_patch_index  Surface_patch_index;
  typedef typename MD::Index                Index;
  
  typedef GT                   Geom_traits;
  typedef typename GT::Point_3 Point;
  
  typedef Point*           Point_container;
  typedef Point*           Point_iterator;
  typedef const Point*     Point_const_iterator;

public:
  void invalidate_circumcenter()
  {
    if (circumcenter_) {
      delete circumcenter_;
      circumcenter_ = NULL;
    }
  }
public:
  // Constructors
  Compact_mesh_cell_base_3()
    : surface_index_table_()
    , surface_center_table_()
    , circumcenter_(NULL)
#ifdef CGAL_INTRUSIVE_LIST
    , next_intrusive_()
    , previous_intrusive_()
#endif
    , surface_center_index_table_()
    , sliver_value_(FT(0.))
    , subdomain_index_()
    , bits_(0)
    , sliver_cache_validity_(false)
  {
  }

  Compact_mesh_cell_base_3(const Compact_mesh_cell_base_3& rhs) 
    : circumcenter_(NULL)
#ifdef CGAL_INTRUSIVE_LIST
    , next_intrusive_(rhs.next_intrusive_)
    , previous_intrusive_(rhs.previous_intrusive_)
#endif
    , sliver_value_(rhs.sliver_value_)
    , subdomain_index_(rhs.subdomain_index_)
    , bits_(0)
    , sliver_cache_validity_(false)
  {
    for(int i=0; i <4; i++){
      surface_index_table_[i] = rhs.surface_index_table_[i];
      surface_center_table_[i]= rhs.surface_center_table_[i];
      surface_center_index_table_[i] = rhs.surface_center_index_table_[i];
    }
  }
  
  Compact_mesh_cell_base_3 (Vertex_handle v0,
                            Vertex_handle v1,
                            Vertex_handle v2,
                            Vertex_handle v3)
    : surface_index_table_()
    , surface_center_table_()
    , circumcenter_(NULL)
#ifdef CGAL_INTRUSIVE_LIST
    , next_intrusive_()
    , previous_intrusive_()
#endif
    , surface_center_index_table_()
    , sliver_value_(FT(0.))
    , subdomain_index_() 
    , bits_(0) 
    , sliver_cache_validity_(false)
  {
    set_vertices(v0, v1, v2, v3);
  }


  Compact_mesh_cell_base_3 (Vertex_handle v0,
                            Vertex_handle v1,
                            Vertex_handle v2,
                            Vertex_handle v3,
                            Cell_handle n0,
                            Cell_handle n1,
                            Cell_handle n2,
                            Cell_handle n3)
    : surface_index_table_()
    , surface_center_table_()
    , circumcenter_(NULL)
#ifdef CGAL_INTRUSIVE_LIST
    , next_intrusive_()
    , previous_intrusive_()
#endif
    , surface_center_index_table_()
    , sliver_value_(FT(0.))
    , subdomain_index_()
    , bits_(0)
    , sliver_cache_validity_(false)
  {
    set_neighbors(n0, n1, n2, n3);
    set_vertices(v0, v1, v2, v3);
  }

  ~Compact_mesh_cell_base_3()
  {
    if(circumcenter_ != NULL){
      delete circumcenter_;
    }
  }

  // ACCESS FUNCTIONS

  Vertex_handle vertex(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3 );
    return V[i];
  }

  bool has_vertex(Vertex_handle v) const
  {
    return (V[0] == v) || (V[1] == v) || (V[2]== v) || (V[3]== v);
  }

  bool has_vertex(Vertex_handle v, int & i) const
  {
    if (v == V[0]) { i = 0; return true; }
    if (v == V[1]) { i = 1; return true; }
    if (v == V[2]) { i = 2; return true; }
    if (v == V[3]) { i = 3; return true; }
    return false;
  }

  int index(Vertex_handle v) const
  {
    if (v == V[0]) { return 0; }
    if (v == V[1]) { return 1; }
    if (v == V[2]) { return 2; }
    CGAL_triangulation_assertion( v == V[3] );
    return 3;
  }

  Cell_handle neighbor(int i) const
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3);
    return N[i];
  }

  bool has_neighbor(Cell_handle n) const
  {
    return (N[0] == n) || (N[1] == n) || (N[2] == n) || (N[3] == n);
  }

  bool has_neighbor(Cell_handle n, int & i) const
  {
    if(n == N[0]){ i = 0; return true; }
    if(n == N[1]){ i = 1; return true; }
    if(n == N[2]){ i = 2; return true; }
    if(n == N[3]){ i = 3; return true; }
    return false;
  }

  int index(Cell_handle n) const
  {
    if (n == N[0]) return 0;
    if (n == N[1]) return 1;
    if (n == N[2]) return 2;
    CGAL_triangulation_assertion( n == N[3] );
    return 3;
  }

  // SETTING


  void set_neighbor(int i, Cell_handle n)
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3);
    CGAL_triangulation_precondition( this != &*n );
    N[i] = n;
  }


  void set_neighbors()
  {
    N[0] = N[1] = N[2] = N[3] = Cell_handle();
  }

  void set_neighbors(Cell_handle n0, Cell_handle n1,
                     Cell_handle n2, Cell_handle n3)
  {
    CGAL_triangulation_precondition( this != &*n0 );
    CGAL_triangulation_precondition( this != &*n1 );
    CGAL_triangulation_precondition( this != &*n2 );
    CGAL_triangulation_precondition( this != &*n3 );
    N[0] = n0;
    N[1] = n1;
    N[2] = n2;
    N[3] = n3;
  }

  // CHECKING

  // the following trivial is_valid allows
  // the user of derived cell base classes
  // to add their own purpose checking
  bool is_valid(bool = false, int = 0) const
  { return true; }

  // For use by Compact_container.
  void * for_compact_container() const { return N[0].for_compact_container(); }
  void * & for_compact_container()     { return N[0].for_compact_container(); }

  // TDS internal data access functions.
        TDS_data& tds_data()       { return _tds_data; }
  const TDS_data& tds_data() const { return _tds_data; }


  Point_iterator hidden_points_begin() const { return hidden_points_end(); }
  Point_iterator hidden_points_end() const { return NULL; }
  void hide_point (const Point &) const { }


  // We must override the functions that modify the vertices.
  // And if the point inside a vertex is modified, we fail,
  // but there's not much we can do for this now.
  void set_vertex(int i, Vertex_handle v)
  {
    CGAL_triangulation_precondition( i >= 0 && i <= 3);
    invalidate_circumcenter();
    V[i] = v;
  }

  void set_vertices()
  {
    invalidate_circumcenter();
    V[0] = V[1] = V[2] = V[3] = Vertex_handle();
  }

  void set_vertices(Vertex_handle v0, Vertex_handle v1,
                    Vertex_handle v2, Vertex_handle v3)
  {
    invalidate_circumcenter();
    V[0] = v0;
    V[1] = v1;
    V[2] = v2;
    V[3] = v3;
  }

  const Point &
  circumcenter(const Geom_traits& gt = Geom_traits()) const
  {
    if (circumcenter_ == NULL) {
#ifndef CGAL_REGULAR_TRIANGULATION_3_USE_CIRCUMCENTER_CACHE
      circumcenter_ = new Point(gt.construct_circumcenter_3_object()
                                (this->vertex(0)->point(),
                                 this->vertex(1)->point(),
                                 this->vertex(2)->point(),
                                 this->vertex(3)->point()));
#else // CGAL_REGULAR_TRIANGULATION_3_USE_CIRCUMCENTER_CACHE
      circumcenter_ = new Point(gt.construct_weighted_circumcenter_3_object()
                                (this->vertex(0)->point(),
                                 this->vertex(1)->point(),
                                 this->vertex(2)->point(),
                                 this->vertex(3)->point()));
#endif // CGAL_REGULAR_TRIANGULATION_3_USE_CIRCUMCENTER_CACHE
    }
    return *circumcenter_;
  }


  // Returns the index of the cell of the input complex that contains the cell
  Subdomain_index subdomain_index() const { return subdomain_index_; }
  
  // Sets the index of the cell of the input complex that contains the cell
  void set_subdomain_index(const Subdomain_index& index)
  { subdomain_index_ = index; }  
  
  void set_sliver_value(const FT& value)
  { 
    sliver_cache_validity_ = true;
    sliver_value_ = value;
  }
  const FT& sliver_value() const { return sliver_value_; }
  bool is_cache_valid() const { return sliver_cache_validity_; }
  void reset_cache_validity() const { sliver_cache_validity_ = false;  }

  /// Set surface index of \c facet to \c index
  void set_surface_patch_index(const int facet, const Surface_patch_index& index)
  {
    CGAL_precondition(facet>=0 && facet<4);
    surface_index_table_[facet] = index;
  }

  /// Returns surface index of facet \c facet
  Surface_patch_index surface_patch_index(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return surface_index_table_[facet];
  }

  /// Marks \c facet as visited
  void set_facet_visited (const int facet)
  {
    CGAL_precondition(facet>=0 && facet <4);
    bits_ |= (1 << facet);
  }

  /// Marks \c facet as not visited
  void reset_visited (const int facet)
  {
    CGAL_precondition(facet>=0 && facet<4);
    bits_ &= (15 & ~(1 << facet));
  }

  /// Returns \c true if \c facet is marked as visited
  bool is_facet_visited (const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return ( (bits_ & (1 << facet)) != 0 );
  }

  /// Sets surface center of \c facet to \c point
  void set_facet_surface_center(const int facet, const Point& point)
  {
    CGAL_precondition(facet>=0 && facet<4);
    surface_center_table_[facet] = point;
  }

  /// Returns surface center of \c facet
  Point get_facet_surface_center(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return surface_center_table_[facet];
  }

  /// Sets surface center index of \c facet to \c index
  void set_facet_surface_center_index(const int facet, const Index& index)
  {
    CGAL_precondition(facet>=0 && facet<4);
    surface_center_index_table_[facet] = index;
  }

  /// Returns surface center of \c facet
  Index get_facet_surface_center_index(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return surface_center_index_table_[facet];
  }

  /// Returns true if facet lies on a surface patch
  bool is_facet_on_surface(const int facet) const
  {
    CGAL_precondition(facet>=0 && facet<4);
    return ( Surface_patch_index() != surface_index_table_[facet]);
  }
  
  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;
  
  void set_surface_index(const int facet, const Surface_index& index)
  { set_surface_patch_index(facet,index); }
  
  /// Returns surface index of facet \c facet
  Surface_index surface_index(const int facet) const
  { return surface_patch_index(facet); }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------

  static
  std::string io_signature()
  {
    return 
      Get_io_signature<Subdomain_index>()() + "+" +
      Get_io_signature<Regular_triangulation_cell_base_3<Geom_traits> >()()
      + "+(" + Get_io_signature<Surface_patch_index>()() + ")[4]";
  }

#ifdef CGAL_INTRUSIVE_LIST
public:
  Cell_handle next_intrusive() const { return next_intrusive_; }
  void set_next_intrusive(Cell_handle c)
  { 
    next_intrusive_ = c; 
  }
   
  Cell_handle previous_intrusive() const { return previous_intrusive_; }
  void set_previous_intrusive(Cell_handle c)
  { 
    previous_intrusive_ = c; 
  }
#endif // CGAL_INTRUSIVE_LIST


private:


  /// Stores surface_index for each facet of the cell
  CGAL::cpp11::array<Surface_patch_index, 4> surface_index_table_;
  /// Stores surface center of each facet of the cell
  CGAL::cpp11::array<Point, 4> surface_center_table_;
  /// Stores surface center index of each facet of the cell

  CGAL::cpp11::array<Cell_handle, 4> N;
  CGAL::cpp11::array<Vertex_handle, 4> V;

  mutable Point * circumcenter_;

#ifdef CGAL_INTRUSIVE_LIST
  Cell_handle next_intrusive_, previous_intrusive_;
#endif

  CGAL::cpp11::array<Index, 4> surface_center_index_table_;
  /// Stores visited facets (4 first bits)

  //  Point_container _hidden;

  FT sliver_value_;

  // The index of the cell of the input complex that contains me
  Subdomain_index subdomain_index_;

  TDS_data      _tds_data;
  char bits_;
  mutable bool sliver_cache_validity_;


};  // end class Compact_mesh_cell_base_3



template < class GT, class MT, class Cb >
std::istream&
operator>>(std::istream &is,
           Compact_mesh_cell_base_3<GT, MT, Cb> &c)
{
  typename Compact_mesh_cell_base_3<GT, MT, Cb>::Subdomain_index index;
  if(is_ascii(is))
    is >> index;
  else
    read(is, index);
  if(is) {
    c.set_subdomain_index(index);
    for(int i = 0; i < 4; ++i)
    {
      typename Compact_mesh_cell_base_3<GT, MT, Cb>::Surface_patch_index i2;
      if(is_ascii(is))
        is >> i2;
      else
      {
        read(is, i2);
      }
      c.set_surface_patch_index(i, i2);
    }
  }
  return is;
}

template < class GT, class MT, class Cb >
std::ostream&
operator<<(std::ostream &os,
           const Compact_mesh_cell_base_3<GT, MT, Cb> &c)
{
  if(is_ascii(os))
     os << c.subdomain_index();
  else
    write(os, c.subdomain_index());
  for(int i = 0; i < 4; ++i)
  {
    if(is_ascii(os))
      os << ' ' << c.surface_patch_index(i);
    else
      write(os, c.surface_patch_index(i));
  }
  return os;
}


// Specialization for void.
template <typename GT, typename MD>
class Compact_mesh_cell_base_3<GT, MD, void>
{
public:
  typedef internal::Dummy_tds_3                         Triangulation_data_structure;
  typedef Triangulation_data_structure::Vertex_handle   Vertex_handle;
  typedef Triangulation_data_structure::Cell_handle     Cell_handle;
  template <typename TDS2>
  struct Rebind_TDS { typedef Compact_mesh_cell_base_3<GT, MD, TDS2> Other; };
};

}  // end namespace CGAL


#endif // CGAL_COMPACT_MESH_CELL_BASE_3_H
