// Copyright (c) 2007  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H
#define CGAL_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H

#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Lightweight_vector_3.h>
#include <CGAL/property_map.h>
#include <CGAL/surface_reconstruction_points_assertions.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>
#include <CGAL/centroid.h>

CGAL_BEGIN_NAMESPACE


/// The Reconstruction_vertex_base_3 class is the default
/// vertex class of the Reconstruction_triangulation_3 class.
///
/// It provides the interface requested by the Poisson_reconstruction_function class:
/// - Each vertex stores a normal vector.
/// - A vertex is either an input point or a Steiner point added by Delaunay refinement.
/// - In order to solve a linear system over the triangulation, a vertex may be constrained
///   or not (i.e. may contribute to the right or left member of the linear system),
///   and has a unique index.
///
/// @heading Parameters:
/// @param Gt   Geometric traits class / Point_3 is a typedef to Point_with_normal_3.
/// @param Cb   Vertex base class, model of TriangulationVertexBase_3.

template < typename Gt,
           typename Vb = Triangulation_vertex_base_3<Gt> >
class Reconstruction_vertex_base_3 : public Vb
{
// Public types
public:

  /// Geometric traits class / Point_3 is a typedef to Point_with_normal_3.
  typedef Gt Geom_traits;

  // Repeat Triangulation_vertex_base_3 public types
  /// @cond SKIP_IN_MANUAL
  typedef typename Vb::Cell_handle Cell_handle;
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other                       Vb2;
    typedef Reconstruction_vertex_base_3<Geom_traits, Vb2> Other;
  };
  /// @endcond

  // Geometric types
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Vector_3 Vector;           ///< typedef to Vector_3
  typedef typename Geom_traits::Point_3 Point;             ///< typedef to Point_with_normal_3
  typedef typename Geom_traits::Point_3 Point_with_normal; ///< typedef to Point_with_normal_3

// data members
private:

  // TODO: reduce memory footprint
  FT m_f; // value of the implicit function // float precise enough?
  bool m_constrained; // is vertex constrained? // combine constrained and type
  unsigned char m_type; // INPUT or STEINER
  unsigned int m_index; // index in matrix (to be stored outside)

// Public methods
public:

  Reconstruction_vertex_base_3()
    : Vb()
  {
    m_f = (FT)0.0;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
  }

  Reconstruction_vertex_base_3(const Point_with_normal& p)
    : Vb(p)
  {
    m_f = 0.0f;
    m_type = 0;
    m_constrained = false;
    m_index = 0;

  }

  Reconstruction_vertex_base_3(const Point_with_normal& p, Cell_handle c)
    : Vb(p,c)
  {
    m_f = (FT)0.0;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
  }

  Reconstruction_vertex_base_3(Cell_handle c)
    : Vb(c)
  {
    m_f = (FT)0.0;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
  }

  /// Is vertex constrained, i.e.
  /// does it contribute to the right or left member of the linear system?
  /// Default value is false.
  bool  constrained() const { return m_constrained; }
  bool& constrained()       { return m_constrained; }

  /// Gets/sets the value of the implicit function.
  /// Default value is 0.0.
  FT  f() const { return m_f; }
  FT& f()       { return m_f; }

  /// Gets/sets the type = INPUT or STEINER.
  unsigned char  type() const { return m_type; }
  unsigned char& type()       { return m_type; }

  /// Gets/sets the index in matrix.
  unsigned int  index() const { return m_index; }
  unsigned int& index()       { return m_index; }

  /// Gets/sets normal vector.
  /// Default value is null vector.
  const Vector& normal() const { return this->point().normal(); }
  Vector&       normal()       { return this->point().normal(); }

// Private methods
private:

    /// Copy constructor and operator =() are not implemented.
    Reconstruction_vertex_base_3(const Reconstruction_vertex_base_3& toCopy);
    Reconstruction_vertex_base_3& operator =(const Reconstruction_vertex_base_3& toCopy);

}; // end of Reconstruction_vertex_base_3


/// Helper class:
/// Reconstruction_triangulation_default_geom_traits_3
/// changes in a geometric traits class the Point_3 type to
/// Point_with_normal_3<BaseGt>.
///
/// @heading Parameters:
/// @param BaseGt   Geometric traits class.
template <class BaseGt>
struct Reconstruction_triangulation_default_geom_traits_3 : public BaseGt
{
  typedef Point_with_normal_3<BaseGt> Point_3;
};


/// The Reconstruction_triangulation_3 class
/// provides the interface requested by the Poisson_reconstruction_function class:
/// - Each vertex stores a normal vector.
/// - A vertex is either an input point or a Steiner point added by Delaunay refinement.
/// - In order to solve a linear system over the triangulation, a vertex may be constrained
///   or not (i.e. may contribute to the right or left member of the linear system),
///   and has a unique index.
/// The vertex class must derive from Reconstruction_vertex_base_3.
///
/// @heading Parameters:
/// @param BaseGt   Geometric traits class.
/// @param Gt       Geometric traits class / Point_3 is a typedef to Point_with_normal_3<BaseGt>.
/// @param Tds      Model of TriangulationDataStructure_3. The vertex class
///                 must derive from Reconstruction_vertex_base_3.

template <class BaseGt,
          class Gt = Reconstruction_triangulation_default_geom_traits_3<BaseGt>,
          class Tds = Triangulation_data_structure_3<Reconstruction_vertex_base_3<Gt> > >
class Reconstruction_triangulation_3 : public Delaunay_triangulation_3<Gt,Tds>
{
// Private types
private:

  // Base class
  typedef Delaunay_triangulation_3<Gt,Tds>  Base;

  // Auxiliary class to build an iterator over normals.
  template <class Node>
  struct Project_normal {
    typedef Node                  argument_type;
    typedef typename Node::Vector Vector;
    typedef Vector                result_type;
    Vector&       operator()(Node& x)       const { return x.normal(); }
    const Vector& operator()(const Node& x) const { return x.normal(); }
  };

  // Auxiliary class to build an iterator over input points.
  class Is_steiner_point
  {
  public:
      typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;

      bool operator()(const Finite_vertices_iterator& v) const
      {
        return (v->type() == Reconstruction_triangulation_3::STEINER);
      }
  };

// Public types
public:

  /// Geometric traits class / Point_3 is a typedef to Point_with_normal_3<BaseGt>.
  typedef Gt  Geom_traits;

  // Repeat base class' types
  /// @cond SKIP_IN_MANUAL
  typedef Tds Triangulation_data_structure;
  typedef typename Base::Segment      Segment;
  typedef typename Base::Triangle     Triangle;
  typedef typename Base::Tetrahedron  Tetrahedron;
  typedef typename Base::Line         Line;
  typedef typename Base::Ray          Ray;
  typedef typename Base::Object       Object;
  typedef typename Base::Cell_handle   Cell_handle;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell   Cell;
  typedef typename Base::Vertex Vertex;
  typedef typename Base::Facet  Facet;
  typedef typename Base::Edge   Edge;
  typedef typename Base::Cell_circulator  Cell_circulator;
  typedef typename Base::Facet_circulator Facet_circulator;
  typedef typename Base::Cell_iterator    Cell_iterator;
  typedef typename Base::Facet_iterator   Facet_iterator;
  typedef typename Base::Edge_iterator    Edge_iterator;
  typedef typename Base::Vertex_iterator  Vertex_iterator;
  typedef typename Base::Point_iterator Point_iterator;
  typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Base::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Base::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Base::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Base::All_cells_iterator       All_cells_iterator;
  typedef typename Base::Locate_type Locate_type;
  /// @endcond

  // Geometric types
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Vector_3 Vector; ///< typedef to Vector_3<BaseGt>
  typedef typename Geom_traits::Point_3 Point;  ///< typedef to Point_with_normal_3<BaseGt>
  typedef typename Geom_traits::Point_3 Point_with_normal; ///< Point_with_normal_3<BaseGt>
  typedef typename Geom_traits::Sphere_3 Sphere;

  /// Iterator over all normals.
  typedef Iterator_project<Finite_vertices_iterator,
                           Project_normal<Vertex> >  Normal_iterator;

  /// Point type
  enum Point_type {
    INPUT,    ///< Input point.
    STEINER   ///< Steiner point created by Delaunay refinement.
  };

  /// Iterator over input vertices.
  typedef Filter_iterator<Finite_vertices_iterator, Is_steiner_point>
                                                    Input_vertices_iterator;

  /// Iterator over input points.
  typedef Iterator_project<Input_vertices_iterator,
                           Project_point<Vertex> >  Input_point_iterator;

// Public methods
public:

  /// Default constructor.
  Reconstruction_triangulation_3()
  {
  }

  // Default copy constructor and operator =() are fine.

  // Repeat base class' public methods used below
  /// @cond SKIP_IN_MANUAL
  Base::points_begin;
  Base::points_end;
  Base::number_of_vertices;
  Base::finite_vertices_begin;
  Base::finite_vertices_end;
  Base::geom_traits;
  /// @endcond

  /// Gets first iterator over finite vertices normals.
  Normal_iterator normals_begin()
  {
      return Normal_iterator(finite_vertices_begin());
  }
  /// Gets past-the-end iterator over finite vertices normals.
  Normal_iterator normals_end()
  {
      return Normal_iterator(finite_vertices_end());
  }

  /// Gets first iterator over input vertices.
  Input_vertices_iterator input_vertices_begin() const
  {
      return Input_vertices_iterator(finite_vertices_end(), Is_steiner_point(),
                                     finite_vertices_begin());
  }
  /// Gets past-the-end iterator over input vertices.
  Input_vertices_iterator input_vertices_end() const
  {
      return Input_vertices_iterator(finite_vertices_end(), Is_steiner_point());
  }

  /// Gets iterator over the first input point.
  Input_point_iterator input_points_begin() const
  {
      return Input_point_iterator(input_vertices_begin());
  }
  /// Gets past-the-end iterator over the input points.
  Input_point_iterator input_points_end() const
  {
      return Input_point_iterator(input_vertices_end());
  }

  /// Gets the bounding sphere of all points.
  Sphere bounding_sphere() const
  {
    Min_sphere_d< CGAL::Optimisation_d_traits_3<Gt> >
      ms3(points_begin(), points_end()); // for all points
    return Sphere(ms3.center(), ms3.squared_radius());
  }

  /// Gets the bounding sphere of input points.
  Sphere input_points_bounding_sphere() const
  {
    Min_sphere_d< CGAL::Optimisation_d_traits_3<Gt> >
      input_points_ms3(input_points_begin(), input_points_end()); // for input points
    return Sphere(input_points_ms3.center(), input_points_ms3.squared_radius());
  }

  /// Insert point in the triangulation.
  /// Default type is INPUT.
  Vertex_handle insert(const Point_with_normal& p,
                       Point_type type = INPUT,
                       Cell_handle start = Cell_handle())
  {
    Vertex_handle v = Base::insert(p, start);
    v->type() = type;
    return v;
  }

  /// Insert points in the triangulation using a spatial sort.
  /// Default type is INPUT.
  ///
  /// @commentheading Precondition:
  /// InputIterator value_type must be convertible to Point_with_normal.
  ///
  /// @param first Iterator over first point to add.
  /// @param beyond Past-the-end iterator to add.
  /// @return the number of inserted points.
  template < class InputIterator >
  int insert(InputIterator first, InputIterator beyond,
             Point_type type = INPUT)
  {
    int n = number_of_vertices();

    // spatial sorting
    std::vector<Point_with_normal> points (first, beyond);
    std::random_shuffle (points.begin(), points.end());
    spatial_sort (points.begin(), points.end(), geom_traits());

    Cell_handle hint;
    for (typename std::vector<Point_with_normal>::const_iterator p = points.begin();
         p != points.end(); ++p)
    {
      Vertex_handle v = insert(*p, type, hint);
      hint = v->cell();
    }

    return number_of_vertices() - n;
  }

  /// Delaunay refinement callback:
  /// insert STEINER point in the triangulation.
  template <class CellIt>
  Vertex_handle
  insert_in_hole(const Point_with_normal & p, CellIt cell_begin, CellIt cell_end,
	         Cell_handle begin, int i,
                 Point_type type = STEINER)
  {
      Vertex_handle v = Base::insert_in_hole(p, cell_begin, cell_end, begin, i);
      v->type() = type;
      return v;
  }

  /// Index unconstrained vertices following the order of Finite_vertices_iterator.
  /// @return the number of unconstrained vertices.
  unsigned int index_unconstrained_vertices()
  {
    unsigned int index = 0;
    for (Finite_vertices_iterator v = finite_vertices_begin();
         v != finite_vertices_end();
         v++)
    {
      if(!v->constrained())
        v->index() = index++;
    }
    return index;
  }

}; // end of Reconstruction_triangulation_3


CGAL_END_NAMESPACE

#endif // CGAL_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H
