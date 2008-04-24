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
#include <CGAL/surface_reconstruction_assertions.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/boost/graph/properties.h>

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
  typedef typename Point_with_normal::Normal Normal; ///< Model of OrientedNormal_3 concept.

// data members
private:
  FT m_f;               // value of the implicit function
                        // PA: should we make a separate type instead?
                        // (so that the user can decide to run in float or double mode)
  bool m_constrained;   // is vertex constrained?
  unsigned char m_type; // INPUT or STEINER
  unsigned int m_index; // index in matrix

// Public methods
public:

  Implicit_fct_delaunay_triangulation_vertex_base_3()
    : Vb()
  {
    m_f = (FT)0.0;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
  }

  Implicit_fct_delaunay_triangulation_vertex_base_3(const Point& p)
    : Vb(p)
  {
    m_f = 0.0f;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
  }

  Implicit_fct_delaunay_triangulation_vertex_base_3(const Point& p, Cell_handle c)
    : Vb(p,c)
  {
    m_f = (FT)0.0;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
  }

  Implicit_fct_delaunay_triangulation_vertex_base_3(Cell_handle c)
    : Vb(c)
  {
    m_f = (FT)0.0;
    m_type = 0;
    m_constrained = false;
    m_index = 0;
  }

  /// is vertex constrained?
  bool  constrained() const { return m_constrained; }
  bool& constrained()       { return m_constrained; }

  /// Get/set the value of the implicit function.
  FT  f() const { return m_f; }
  FT& f()       { return m_f; }

  /// Get/set the type = INPUT or STEINER.
  unsigned char  type() const { return m_type; }
  unsigned char& type()       { return m_type; }

  /// Get/set the index in matrix.
  unsigned int  index() const { return m_index; }
  unsigned int& index()       { return m_index; }

  /// Get/set normal (vector + orientation).
  const Normal& normal() const { return this->point().normal(); }
  Normal&       normal()       { return this->point().normal(); }

// Private methods
private:

    /// Copy constructor and operator =() are not implemented.
    Implicit_fct_delaunay_triangulation_vertex_base_3(const Implicit_fct_delaunay_triangulation_vertex_base_3& toCopy);
    Implicit_fct_delaunay_triangulation_vertex_base_3& operator =(const Implicit_fct_delaunay_triangulation_vertex_base_3& toCopy);

}; // end of Implicit_fct_delaunay_triangulation_vertex_base_3


/// Helper class: Implicit_fct_delaunay_triangulation_default_geom_traits_3
/// changes in a geometric traits class the Point_3 type to Point_with_normal_3.
///
/// @heading Parameters:
/// @param BaseGt   Kernel's regular geometric traits.
template <class BaseGt>
struct Implicit_fct_delaunay_triangulation_default_geom_traits_3 : public BaseGt
{
  typedef Point_with_normal_3<BaseGt> Point_3;
};


/// The Implicit_fct_delaunay_triangulation_3 class is the default implementation
/// of the ImplicitFctDelaunayTriangulation_3 concept.
/// It provides the interface requested by the Poisson_implicit_function class.
///
/// TODO: Test speed if using Triangulation_hierarchy_3
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

  typedef Delaunay_triangulation_3<Gt,Tds>  Base;

  // Auxiliary class to build a normals iterator
  template <class Node>
  struct Project_normal {
    typedef Node                  argument_type;
    typedef typename Node::Normal Normal;
    typedef Normal                result_type;
    typedef Arity_tag<1> Arity;
    Normal&       operator()(Node& x)       const { return x.normal(); }
    const Normal& operator()(const Node& x) const { return x.normal(); }
  };

// Public types
public:

  // Repeat Delaunay_triangulation_3 public types
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
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
  typedef typename Geom_traits::Sphere_3 Sphere;

  /// The geometric traits class's Point_3 type is a model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point;             ///< Model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point_with_normal; ///< Model of PointWithNormal_3
  typedef typename Point_with_normal::Normal Normal; ///< Model of OrientedNormal_3 concept.

  /// Iterator over normals
  typedef Iterator_project<Finite_vertices_iterator,
                           Project_normal<Vertex> >  Normal_iterator;

  /// Point type
  static const unsigned char INPUT = 0;
  static const unsigned char STEINER = 1;

// Data members
private:

  // Indicate if m_barycenter, m_bounding_box and m_diameter_standard_deviation below are valid
  mutable bool m_bounding_box_is_valid;
  mutable Iso_cuboid_3 m_bounding_box; // Triangulation's bounding box
  mutable Point m_barycenter; // Triangulation's barycenter
  mutable FT m_diameter_standard_deviation; // Triangulation's standard deviation

// Public methods
public:

  /// Default constructor.
  Implicit_fct_delaunay_triangulation_3()
  {
    m_bounding_box_is_valid = false;
  }

  // Default copy constructor and operator =() are fine.

  // Repeat Delaunay_triangulation_3 public methods used below
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

  /// Get the bounding box.
  Iso_cuboid_3 bounding_box() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_bounding_box;
  }

  /// Get bounding sphere.
  Sphere bounding_sphere() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    // Center point
    FT mx = 0.5 * (m_bounding_box.xmax() + m_bounding_box.xmin());
    FT my = 0.5 * (m_bounding_box.ymax() + m_bounding_box.ymin());
    FT mz = 0.5 * (m_bounding_box.zmax() + m_bounding_box.zmin());
    Point center(mx,my,mz);

    // Squared radius
    FT dx = m_bounding_box.xmax() - m_bounding_box.xmin();
    FT dy = m_bounding_box.ymax() - m_bounding_box.ymin();
    FT dz = m_bounding_box.zmax() - m_bounding_box.zmin();
    FT squared_radius = dx*dx + dy*dy + dz*dz;

    return Sphere(center, squared_radius);
  }

  /// Get points barycenter.
  Point barycenter() const
  {
    if (!m_bounding_box_is_valid)
      update_bounding_box();

    return m_barycenter;
  }

  /// Get the standard deviation of the distance to barycenter.
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

  /// Insert point to the triangulation.
  Vertex_handle insert(const Point& p,
                       unsigned char type /* INPUT or STEINER */,
                       Cell_handle start = Cell_handle())
  {
    Vertex_handle v = Base::insert(p, start);
    v->type() = type;

    invalidate_bounding_box();

    return v;
  }

  /// Insert points to the triangulation using a spatial sort.
  ///
  /// Precondition: the value type of InputIterator must 'Point'.
  ///
  /// @param first First point to add to pdt.
  /// @param beyond Past-the-end point to add to pdt.
  /// @return the number of inserted points.
  template < class InputIterator >
  int insert(InputIterator first, InputIterator beyond,
             unsigned char type /* INPUT or STEINER */)
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

  /// Index all (finite) vertices following the order of Finite_vertices_iterator.
  /// @return the number (finite) of vertices.
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

  /// Recompute barycenter, bounding box, bounding sphere and standard deviation.
  void update_bounding_box() const
  {
    // Update bounding box and barycenter.
    // TODO: we should use the functions in PCA component instead.
    FT xmin,xmax,ymin,ymax,zmin,zmax;
    xmin = ymin = zmin =  1e38;
    xmax = ymax = zmax = -1e38;
    Vector v = NULL_VECTOR;
    FT norm = 0;
    CGAL_surface_reconstruction_assertion(points_begin() != points_end());
    for (Point_iterator it = points_begin(); it != points_end(); it++)
    {
        const Point& p = *it;

        // update bbox
        xmin = (std::min)(p.x(),xmin);
        ymin = (std::min)(p.y(),ymin);
        zmin = (std::min)(p.z(),zmin);
        xmax = (std::max)(p.x(),xmax);
        ymax = (std::max)(p.y(),ymax);
        zmax = (std::max)(p.z(),zmax);

        // update barycenter
        v = v + (p - ORIGIN);
        norm += 1;
    }
    //
    Point p(xmin,ymin,zmin);
    Point q(xmax,ymax,zmax);
    m_bounding_box = Iso_cuboid_3(p,q);
    //
    m_barycenter = ORIGIN + v / norm;

    /// Compute standard deviation of the distance to barycenter
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

}; // end of Implicit_fct_delaunay_triangulation_3


/// Helper class: type of the "vertex_point" property map
/// of an Implicit_fct_delaunay_triangulation_3 object.
template <class BaseGt, class Gt, class Tds>
class Implicit_fct_delaunay_triangulation_vertex_point_const_map 
{
public:
    typedef Implicit_fct_delaunay_triangulation_3<BaseGt,Gt,Tds> Triangulation;
    typedef typename Gt::Point_3 Point_3;  

    // Property maps required types
    typedef boost::readable_property_map_tag                    category;
    typedef Point_3                                             value_type;
    typedef value_type                                          reference;
    typedef typename Triangulation::Finite_vertices_iterator    key_type;

    Implicit_fct_delaunay_triangulation_vertex_point_const_map(const Triangulation&) {}

    /// Free function to access the map elements.
    friend inline 
    reference 
    get(const Implicit_fct_delaunay_triangulation_vertex_point_const_map&, key_type v)
    {
      return v->point();
    }
};

/// Free function to get the "vertex_point" property map
/// of an Implicit_fct_delaunay_triangulation_3 object.
template <class BaseGt, class Gt, class Tds>
inline
Implicit_fct_delaunay_triangulation_vertex_point_const_map<BaseGt,Gt,Tds> 
get(vertex_point_t, const Implicit_fct_delaunay_triangulation_3<BaseGt,Gt,Tds>& tr) 
{
  Implicit_fct_delaunay_triangulation_vertex_point_const_map<BaseGt,Gt,Tds> aMap(tr);
  return aMap;
}


/// Helper class: type of the "vertex_normal" property map
/// of an Implicit_fct_delaunay_triangulation_3 object.
template <class BaseGt, class Gt, class Tds>
class Implicit_fct_delaunay_triangulation_vertex_normal_map 
  : public boost::put_get_helper< typename Gt::Point_3::Normal&, 
                                  Implicit_fct_delaunay_triangulation_vertex_normal_map<BaseGt,Gt,Tds> >
{
public:
    typedef Implicit_fct_delaunay_triangulation_3<BaseGt,Gt,Tds> Triangulation;
    typedef typename Gt::Point_3::Normal Normal;  

    // Property maps required types
    typedef boost::lvalue_property_map_tag                      category;
    typedef Normal                                              value_type;
    typedef Normal&                                             reference;
    typedef typename Triangulation::Finite_vertices_iterator    key_type;

    Implicit_fct_delaunay_triangulation_vertex_normal_map(const Triangulation&) {}

    /// Access the map elements.
    reference operator[](key_type v) const { return v->normal(); }
};

/// Free function to get the "vertex_normal" property map
/// of an Implicit_fct_delaunay_triangulation_3 object.
template <class BaseGt, class Gt, class Tds>
inline
Implicit_fct_delaunay_triangulation_vertex_normal_map<BaseGt,Gt,Tds> 
get(boost::vertex_normal_t, const Implicit_fct_delaunay_triangulation_3<BaseGt,Gt,Tds>& tr) 
{
  Implicit_fct_delaunay_triangulation_vertex_normal_map<BaseGt,Gt,Tds> aMap(tr);
  return aMap;
}


/// Helper class: type of the "vertex_index" property map
/// of an Implicit_fct_delaunay_triangulation_3 object.
template <class BaseGt, class Gt, class Tds>
class Implicit_fct_delaunay_triangulation_vertex_index_map 
  : public boost::put_get_helper< unsigned int&, 
                                  Implicit_fct_delaunay_triangulation_vertex_index_map<BaseGt,Gt,Tds> >
{
public:
    typedef Implicit_fct_delaunay_triangulation_3<BaseGt,Gt,Tds> Triangulation;

    // Property maps required types
    typedef boost::lvalue_property_map_tag                      category;
    typedef unsigned int                                        value_type;
    typedef unsigned int&                                       reference;
    typedef typename Triangulation::Finite_vertices_iterator    key_type;

    Implicit_fct_delaunay_triangulation_vertex_index_map(const Triangulation&) {}

    /// Access the map elements.
    reference operator[](key_type v) const { return v->index(); }
};

/// Free function to get the "vertex_index" property map
/// of an Implicit_fct_delaunay_triangulation_3 object.
template <class BaseGt, class Gt, class Tds>
inline
Implicit_fct_delaunay_triangulation_vertex_index_map<BaseGt,Gt,Tds> 
get(boost::vertex_index_t, const Implicit_fct_delaunay_triangulation_3<BaseGt,Gt,Tds>& tr) 
{
  Implicit_fct_delaunay_triangulation_vertex_index_map<BaseGt,Gt,Tds> aMap(tr);
  return aMap;
}


CGAL_END_NAMESPACE

#endif // CGAL_IMPLICIT_FCT_DELAUNAY_TRIANGULATION_H
