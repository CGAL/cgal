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
#include <CGAL/surface_reconstruction_assertions.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/centroid.h>

CGAL_BEGIN_NAMESPACE


/// The Implicit_fct_delaunay_triangulation_cell_base_3 class is the default
/// cell class of the Implicit_fct_delaunay_triangulation_3 class.
/// It provides the interface requested by the Poisson_implicit_function class.
///
/// @heading Is Model for the Concepts:
/// Model of the ImplicitFctDelaunayTriangulationCellBase_3 concept.
///
/// @heading Parameters:
/// @param Gt   Geometric traits class / Point_3 is a model of PointWithNormal_3.
/// @param Cb   Cell base class, model of TriangulationCellBase_3.

template < typename Gt,
           typename Cb = Triangulation_cell_base_3<Gt> >
class Implicit_fct_delaunay_triangulation_cell_base_3 : public Cb
{
// Public types
public:

  // Repeat Triangulation_cell_base_3 public types
  /// @cond SKIP_IN_MANUAL
  typedef Gt Geom_traits; ///< Geometric traits class / Point_3 is a model of PointWithNormal_3.
  typedef typename Cb::Cell_handle Cell_handle;
  typedef typename Cb::Vertex_handle Vertex_handle;
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other                    Cb2;
    typedef Implicit_fct_delaunay_triangulation_cell_base_3<Geom_traits,Cb2> Other;
  };
  /// @endcond

// Public methods
public:

  Implicit_fct_delaunay_triangulation_cell_base_3()
    : Cb()
  {
  }

  Implicit_fct_delaunay_triangulation_cell_base_3(Cell_handle c)
    : Cb(c)
  {
  }

  Implicit_fct_delaunay_triangulation_cell_base_3(Vertex_handle v1,
                                                  Vertex_handle v2,
                                                  Vertex_handle v3,
                                                  Vertex_handle v4)
    : Cb(v1,v2,v3,v4)
  {
  }

// Private methods
private:

    /// Copy constructor and operator =() are not implemented.
    Implicit_fct_delaunay_triangulation_cell_base_3(const Implicit_fct_delaunay_triangulation_cell_base_3& toCopy);
    Implicit_fct_delaunay_triangulation_cell_base_3& operator =(const Implicit_fct_delaunay_triangulation_cell_base_3& toCopy);

}; // end of Implicit_fct_delaunay_triangulation_cell_base_3


/// The Implicit_fct_delaunay_triangulation_vertex_base_3 class is the default
/// vertex class of the Implicit_fct_delaunay_triangulation_3 class.
/// It provides the interface requested by the Poisson_implicit_function class.
///
/// @heading Is Model for the Concepts:
/// Model of the ImplicitFctDelaunayTriangulationVertexBase_3 concept.
///
/// @heading Parameters:
/// @param Gt   Geometric traits class / Point_3 is a model of PointWithNormal_3.
/// @param Cb   Vertex base class, model of TriangulationVertexBase_3.

template < typename Gt,
           typename Vb = Triangulation_vertex_base_3<Gt> >
class Implicit_fct_delaunay_triangulation_vertex_base_3 : public Vb
{
// Public types
public:

  // Repeat Triangulation_vertex_base_3 public types
  /// @cond SKIP_IN_MANUAL
  typedef Gt Geom_traits; ///< Geometric traits class / Point_3 is a model of PointWithNormal_3.
  typedef typename Vb::Cell_handle Cell_handle;
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other                       Vb2;
    typedef Implicit_fct_delaunay_triangulation_vertex_base_3<Geom_traits, Vb2> Other;
  };
  /// @endcond

  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;             ///< Model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point_with_normal; ///< Model of PointWithNormal_3
  typedef typename Point_with_normal::Normal Normal; ///< Model of Kernel::Vector_3 concept.

// data members
private:

  FT m_f;               // value of the implicit function
                        // PA: should we make a separate type instead?
                        // (so that the user can decide to run in float or double mode)
  bool m_constrained;   // is vertex constrained?
  unsigned char m_type; // INPUT or STEINER. Default type is INPUT.
  unsigned int m_index; // index in matrix
  double m_average_spacing;   // average spacing
  int m_tag;

// Public methods
public:

  Implicit_fct_delaunay_triangulation_vertex_base_3()
    : Vb()
  {
    m_f = (FT)0.0;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
    m_average_spacing = 0.0;
    m_tag = -1;
  }

  Implicit_fct_delaunay_triangulation_vertex_base_3(const Point& p)
    : Vb(p)
  {
    m_f = 0.0f;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
    m_average_spacing = 0.0;
    m_tag = -1;

  }

  Implicit_fct_delaunay_triangulation_vertex_base_3(const Point& p, Cell_handle c)
    : Vb(p,c)
  {
    m_f = (FT)0.0;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
    m_average_spacing = 0.0;
    m_tag = -1;
  }

  Implicit_fct_delaunay_triangulation_vertex_base_3(Cell_handle c)
    : Vb(c)
  {
    m_f = (FT)0.0;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
    m_average_spacing = 0.0;
    tag = -1;
  }

  /// Is vertex constrained?
  bool  constrained() const { return m_constrained; }
  bool& constrained()       { return m_constrained; }

  /// Get/set the value of the implicit function.
  FT  f() const { return m_f; }
  FT& f()       { return m_f; }

  /// Get/set average spacing at each input point.
  double  average_spacing() const { return m_average_spacing; }
  double& average_spacing()       { return m_average_spacing; }

  /// Get/set the type = INPUT or STEINER.
  unsigned char  type() const { return m_type; }
  unsigned char& type()       { return m_type; }

  /// Get/set the index in matrix.
  unsigned int  index() const { return m_index; }
  unsigned int& index()       { return m_index; }

  /// Get/set normal (vector + orientation).
  const Normal& normal() const { return this->point().normal(); }
  Normal&       normal()       { return this->point().normal(); }

  /// General purpose tag
  int tag() const { return m_tag; }
  int& tag()       { return m_tag; }

// Private methods
private:

    /// Copy constructor and operator =() are not implemented.
    Implicit_fct_delaunay_triangulation_vertex_base_3(const Implicit_fct_delaunay_triangulation_vertex_base_3& toCopy);
    Implicit_fct_delaunay_triangulation_vertex_base_3& operator =(const Implicit_fct_delaunay_triangulation_vertex_base_3& toCopy);

}; // end of Implicit_fct_delaunay_triangulation_vertex_base_3


/// Helper class:
/// Implicit_fct_delaunay_triangulation_default_geom_traits_3
/// changes in a geometric traits class the Point_3 type to
/// a lightweight model of PointWithNormal_3.
///
/// @heading Parameters:
/// @param BaseGt   Kernel's regular geometric traits.
template <class BaseGt>
struct Implicit_fct_delaunay_triangulation_default_geom_traits_3 : public BaseGt
{
  typedef Point_with_normal_3<BaseGt, Lightweight_vector_3<BaseGt> > Point_3;
};


/// The Implicit_fct_delaunay_triangulation_3 class is the default implementation
/// of the ImplicitFctDelaunayTriangulation_3 concept.
/// It provides the interface requested by the Poisson_implicit_function class.
///
/// CAUTION: invalidate_bounding_box() must be called
/// after modifying the points.
///
/// @heading Is Model for the Concepts:
/// Model of the ImplicitFctDelaunayTriangulation_3 concept.
///
/// @heading Parameters:
/// @param BaseGt   Kernel's regular geometric traits.
/// @param Gt       Geometric traits class / Point_3 is a model of PointWithNormal_3.
/// @param Tds      Model of TriangulationDataStructure_3. The cell base class must be
/// a model of ImplicitFctDelaunayTriangulationCellBase_3 and the vertex base class
/// must be a model of ImplicitFctDelaunayTriangulationVertexBase_3.

template <class BaseGt,
          class Gt = Implicit_fct_delaunay_triangulation_default_geom_traits_3<BaseGt>,
          class Tds = Triangulation_data_structure_3<Implicit_fct_delaunay_triangulation_vertex_base_3<Gt>,
                                                     Implicit_fct_delaunay_triangulation_cell_base_3<Gt> > >
class Implicit_fct_delaunay_triangulation_3 : public Delaunay_triangulation_3<Gt,Tds>
{
// Private types
private:

  // Base class
  typedef Delaunay_triangulation_3<Gt,Tds>  Base;

  // Auxiliary class to build an iterator over normals.
  template <class Node>
  struct Project_normal {
    typedef Node                  argument_type;
    typedef typename Node::Normal Normal;
    typedef Normal                result_type;
    Normal&       operator()(Node& x)       const { return x.normal(); }
    const Normal& operator()(const Node& x) const { return x.normal(); }
  };

  // Auxiliary class to build an iterator over input points.
  class Is_steiner_point
  {
  public:
      typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;

      bool operator()(const Finite_vertices_iterator& v) const
      {
        return (v->type() == Implicit_fct_delaunay_triangulation_3::STEINER);
      }
  };

// Public types
public:

  // Repeat base class' types
  /// @cond SKIP_IN_MANUAL
  typedef Tds Triangulation_data_structure;
  typedef Gt  Geom_traits; ///< Geometric traits class / Point_3 is a model of PointWithNormal_3.
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
  typedef typename Geom_traits::Vector_3 Vector;
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;
  typedef typename Geom_traits::Sphere_3 Sphere;

  /// The geometric traits class's Point_3 type is a model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point;             ///< Model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point_with_normal; ///< Model of PointWithNormal_3
  typedef typename Point_with_normal::Normal Normal; ///< Model of Kernel::Vector_3 concept.

  /// Iterator over all normals.
  typedef Iterator_project<Finite_vertices_iterator,
                           Project_normal<Vertex> >  Normal_iterator;

  /// Point type
  static const unsigned char INPUT = 0;
  static const unsigned char STEINER = 1;

  /// Iterator over INPUT vertices.
  typedef Filter_iterator<Finite_vertices_iterator, Is_steiner_point>
                                                    Input_vertices_iterator;

  /// Iterator over INPUT points.
  typedef Iterator_project<Input_vertices_iterator,
                           Project_point<Vertex> >  Input_point_iterator;

// Public methods
public:

  /// Default constructor.
  Implicit_fct_delaunay_triangulation_3()
  {
    m_bounding_box_is_valid = false;
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

  /// Get first iterator over finite vertices normals.
  Normal_iterator normals_begin()
  {
      return Normal_iterator(finite_vertices_begin());
  }
  /// Get past-the-end iterator over finite vertices normals.
  Normal_iterator normals_end()
  {
      return Normal_iterator(finite_vertices_end());
  }

  /// Get first iterator over INPUT vertices.
  Input_vertices_iterator input_vertices_begin() const
  {
      return Input_vertices_iterator(finite_vertices_end(), Is_steiner_point(),
                                     finite_vertices_begin());
  }
  /// Get past-the-end iterator over INPUT vertices.
  Input_vertices_iterator input_vertices_end() const
  {
      return Input_vertices_iterator(finite_vertices_end(), Is_steiner_point());
  }

  /// Get first iterator over INPUT points.
  Input_point_iterator input_points_begin() const
  {
      return Input_point_iterator(input_vertices_begin());
  }
  /// Get past-the-end iterator over INPUT points.
  Input_point_iterator input_points_end() const
  {
      return Input_point_iterator(input_vertices_end());
  }

  /// Get the bounding box of all points.
  Iso_cuboid bounding_box() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_bounding_box;
  }

  /// Get the bounding box of INPUT points.
  Iso_cuboid input_points_bounding_box() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_input_points_bounding_box;
  }

  /// Get the bounding sphere of all points.
  Sphere bounding_sphere() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_bounding_sphere;
  }

  /// Get the bounding sphere of INPUT points.
  Sphere input_points_bounding_sphere() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_input_points_bounding_sphere;
  }

  /// Get the barycenter of all points.
  Point barycenter() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_barycenter;
  }

  /// Get the standard deviation of the distance to barycenter (for all points).
  FT diameter_standard_deviation() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_diameter_standard_deviation;
  }

  /// Update barycenter, bounding box, bounding sphere and standard deviation.
  /// Owner is responsible to call this function after modifying the triangulation.
  void invalidate_bounding_box()
  {
    m_bounding_box_is_valid = false;
  }

  /// Insert point in the triangulation.
  /// Default type is INPUT.
  Vertex_handle insert(const Point& p,
                       unsigned char type = INPUT /* INPUT or STEINER */,
                       Cell_handle start = Cell_handle())
  {
    Vertex_handle v = Base::insert(p, start);

    v->type() = type;
    invalidate_bounding_box();

    return v;
  }

  /// Insert points in the triangulation using a spatial sort.
  /// Default type is INPUT.
  ///
  /// Precondition: the value type of InputIterator must be 'Point'.
  ///
  /// @param first First point to add to pdt.
  /// @param beyond Past-the-end point to add to pdt.
  /// @return the number of inserted points.
  template < class InputIterator >
  int insert(InputIterator first, InputIterator beyond,
             unsigned char type = INPUT /* INPUT or STEINER */)
  {
    int n = number_of_vertices();

    // spatial sorting
    std::vector<Point> points (first, beyond);
    std::random_shuffle (points.begin(), points.end());
    spatial_sort (points.begin(), points.end(), geom_traits());

    Cell_handle hint;
    for (typename std::vector<Point>::const_iterator p = points.begin();
         p != points.end(); ++p)
    {
      Vertex_handle v = insert(*p, type, hint);
      hint = v->cell();
    }

    invalidate_bounding_box();

    return number_of_vertices() - n;
  }

  /// Delaunay refinement callback:
  /// insert STEINER point in the triangulation.
  template <class CellIt>
  Vertex_handle
  insert_in_hole(const Point & p, CellIt cell_begin, CellIt cell_end,
	         Cell_handle begin, int i,
                 unsigned char type = STEINER /* INPUT or STEINER */)
  {
      Vertex_handle v = Base::insert_in_hole(p, cell_begin, cell_end, begin, i);

      v->type() = type;
      invalidate_bounding_box();

      return v;
  }

  /// Index all (finite) vertices following the order of Finite_vertices_iterator.
  /// @return the number of (finite) vertices.
  unsigned int index_vertices()
  {
    unsigned int index = 0;
    for (Finite_vertices_iterator v = finite_vertices_begin();
         v != finite_vertices_end();
         v++)
    {
      v->index() = index++;
    }
    return index;
  }

  /// Index unconstraint vertices following the order of Finite_vertices_iterator.
  /// @return the number of unconstraint vertices.
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

// Private methods:
private:

  /// Compute barycenter, bounding box, bounding sphere and standard deviation.
  void update_bounding_box() const
  {
    typedef CGAL::Min_sphere_d< CGAL::Optimisation_d_traits_3<Gt> > Min_sphere_d;

    if (points_begin() == points_end())
      return;

    // Compute barycenter and bounding boxes
    //
    // LS 06/2008: We should use the functions in PCA component instead.
    //             Unfortunately, the next lines not compile...
    //m_bounding_box = CGAL::bounding_box(points_begin(), points_end());
    //m_barycenter   = CGAL::centroid(points_begin(), points_end());
    //m_input_points_bounding_box = CGAL::bounding_box(input_points_begin(), input_points_end());
    //
    FT xmin,xmax,ymin,ymax,zmin,zmax; // for all points
    xmin = ymin = zmin =  1e38;
    xmax = ymax = zmax = -1e38;
    Vector sum = CGAL::NULL_VECTOR;
    FT nb_points = 0;
    FT input_xmin,input_xmax,input_ymin,input_ymax,input_zmin,input_zmax; // for INPUT points
    input_xmin = input_ymin = input_zmin =  1e38;
    input_xmax = input_ymax = input_zmax = -1e38;
    for (Finite_vertices_iterator it = finite_vertices_begin(); it != finite_vertices_end(); it++)
    {
      const Point& p = it->point();

      // update bounding box of all points
      xmin = (std::min)(p.x(),xmin);
      ymin = (std::min)(p.y(),ymin);
      zmin = (std::min)(p.z(),zmin);
      xmax = (std::max)(p.x(),xmax);
      ymax = (std::max)(p.y(),ymax);
      zmax = (std::max)(p.z(),zmax);

      if (it->type() == INPUT)
      {
        // update bounding box of INPUT points
        input_xmin = (std::min)(p.x(),input_xmin);
        input_ymin = (std::min)(p.y(),input_ymin);
        input_zmin = (std::min)(p.z(),input_zmin);
        input_xmax = (std::max)(p.x(),input_xmax);
        input_ymax = (std::max)(p.y(),input_ymax);
        input_zmax = (std::max)(p.z(),input_zmax);
      }

      // update barycenter of all points
      sum = sum + (p - CGAL::ORIGIN);
      nb_points += 1;
    }
    //
    Point p(xmin,ymin,zmin); // for all points
    Point q(xmax,ymax,zmax);
    m_bounding_box = Iso_cuboid(p,q);
    m_barycenter = CGAL::ORIGIN + sum / nb_points;
    //
    Point input_p(input_xmin,input_ymin,input_zmin); // for INPUT points
    Point input_q(input_xmax,input_ymax,input_zmax);
    m_input_points_bounding_box = Iso_cuboid(input_p,input_q);

    // Compute bounding spheres
    Min_sphere_d ms3(points_begin(), points_end()); // for all points
    m_bounding_sphere = Sphere(ms3.center(), ms3.squared_radius());
    //
    Min_sphere_d input_points_ms3(input_points_begin(), input_points_end()); // for INPUT points
    m_input_points_bounding_sphere = Sphere(input_points_ms3.center(), input_points_ms3.squared_radius());

    // Compute standard deviation of the distance to barycenter (for all points)
    typename Geom_traits::Compute_squared_distance_3 sqd;
    FT sq_radius = 0;
    for (Point_iterator it = points_begin(); it != points_end(); it++)
    {
        sq_radius += sqd(*it, m_barycenter);
    }
    sq_radius /= number_of_vertices();
    m_diameter_standard_deviation = CGAL::sqrt(sq_radius);

    m_bounding_box_is_valid = true;
  }

// Data members
private:

  // Indicate if m_*bounding_box, m_*bounding_sphere, m_barycenter and
  // m_diameter_standard_deviation below are valid.
  mutable bool m_bounding_box_is_valid;

  mutable Iso_cuboid m_bounding_box; // bounding box of all points.
  mutable Iso_cuboid m_input_points_bounding_box; // bounding box of INPUT points.
  mutable Sphere m_bounding_sphere; // bounding sphere of all points.
  mutable Sphere m_input_points_bounding_sphere; // bounding sphere of INPUT points.
  mutable Point m_barycenter; // barycenter of all points.
  mutable FT m_diameter_standard_deviation; // standard deviation of the distance
                                            // to barycenter (for all points).

}; // end of Implicit_fct_delaunay_triangulation_3


CGAL_END_NAMESPACE

#endif // CGAL_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H
