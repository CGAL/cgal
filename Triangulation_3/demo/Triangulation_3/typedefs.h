#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#define CGAL_CONCURRENT_TRIANGULATION_3_PROFILING

#include <vector>  //dynamic array
#include <list>  //linked list

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>

// Added for T3 demo
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>

// Use EPEC as Kernel
// Note: the computation of VD requires exact constructions;
//  while computing the triangulation only requires exact predicates.
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Ddefine field type
typedef Kernel::FT	FT;

typedef Kernel::Vector_3		Vector_3;
typedef Kernel::Direction_3		Direction_3;

//typedef Kernel::Point_3		Point_3;
//typedef Kernel::Vector_3	Vector_3;
//typedef Kernel::Segment_3	Segment_3;
//typedef Kernel::Triangle_3	Triangle_3;

// Added for T3 demo

/*
 * The user has several ways to add his own data in the vertex
 * and cell base classes used by the TDS. He can either:
 * 1. use the classes Triangulation vertex base with info
 *  and Triangulation cell base with info, which allow to
 *  add one data member of a user provided type, and give access to it.
 * 2. derive his own classes from the default base classes
 *  Triangulation ds vertex base, and Triangulation ds cell base
 *  (or the geometric versions typically used by the geometric layer,
 *  Triangulation vertex base, and Triangulation cell base).
 * 3. write his own base classes following the requirements given by the concepts
 *  TriangulationCellBase 3 and TriangulationVertexBase 3
 *  (described in page 2494 and page 2495).
 */
/* add index and color to vertex class */
template < class GT, class Vb=CGAL::Triangulation_vertex_base_3<GT> >
class Vertex_base : public Vb
{
public:
  typedef typename Vb::Point	Point;
  typedef typename Vb::Vertex_handle	Vertex_handle;
  typedef typename Vb::Cell_handle		Cell_handle;

  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other	Vb2;
    typedef Vertex_base< GT, Vb2 >	Other;
  };
  Vertex_base()
   : m_isSelected(false) {}
  Vertex_base(const Point& p)
   : Vb(p), m_isSelected(false) {}
  Vertex_base(const Point& p, Cell_handle c)
   : Vb(p, c), m_isSelected(false) {}

  inline bool isSeled() const {  return m_isSelected;  }
  inline void setSeled(bool flag=true) {  m_isSelected = flag;  }

private:
  bool m_isSelected;  // whether it is selected
};

/*
 * Delaunay_triangulation_3<arg1, arg2, arg3>
 * arg1: a model of the DelaunayTriangulationTraits_3 concept
 * arg2: a model of the TriangulationDataStructure_3 concept
 *       default: Triangulation_data_structure_3
 * arg3: Fast_location or Compact_location (default)
 *       Fast_location offers O(logn) time point location, using additional data structure,
 *     good for fast point locations or random point insertions
 *       Compact_location saves memory by avoiding the separate data structure
 *     and point location is then O(n^(1/3)) time
 */
#ifdef CGAL_CONCURRENT_TRIANGULATION_3
typedef CGAL::Spatial_lock_grid_3<
  CGAL::Tag_priority_blocking>                        Lock_ds;
typedef CGAL::Triangulation_data_structure_3< 
  Vertex_base<Kernel>, 
  CGAL::Triangulation_ds_cell_base_3<>, 
  CGAL::Parallel_tag >	                              Tds;
typedef CGAL::Delaunay_triangulation_3<
  Kernel, Tds, CGAL::Default, Lock_ds>	              DT3;

#else
typedef CGAL::Triangulation_data_structure_3< Vertex_base<Kernel> >	Tds;
typedef CGAL::Delaunay_triangulation_3<
  Kernel, Tds/*, CGAL::Fast_location*/>	                            DT3;
#endif

typedef DT3::Object		Object_3;
typedef DT3::Point		Point_3;
typedef DT3::Segment	Segment_3;
typedef DT3::Ray		Ray_3;
typedef DT3::Triangle	Triangle_3;

typedef DT3::Vertex_handle	Vertex_handle;
typedef DT3::Finite_vertices_iterator	vertices_iterator;
typedef DT3::Edge		Edge;
typedef DT3::Finite_edges_iterator	edges_iterator;
typedef DT3::Facet		Facet;
typedef DT3::Finite_facets_iterator	facets_iterator;
typedef DT3::Cell_handle	Cell_handle;
typedef DT3::Finite_cells_iterator	cells_iterator;

#endif
