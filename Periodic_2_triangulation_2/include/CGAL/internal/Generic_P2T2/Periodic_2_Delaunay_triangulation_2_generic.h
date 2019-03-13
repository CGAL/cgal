#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/internal/Generic_P2T2/Periodic_2_triangulation_vertex_base_2_generic.h>
#include <CGAL/internal/Generic_P2T2/Periodic_2_triangulation_face_base_2_generic.h>
#include <CGAL/utility.h>

#include <utility>
#include <iostream>

namespace CGAL {

template < class Gt,
           class TDS = Triangulation_data_structure_2 <
                         Periodic_2_triangulation_vertex_base_2_generic<Gt>,
                         Periodic_2_triangulation_face_base_2_generic<Gt> > >
class Periodic_2_Delaunay_triangulation_2_generic
{
  typedef Periodic_2_Delaunay_triangulation_2_generic<Gt, TDS>          Self;

public:
  typedef TDS                                  Triangulation_data_structure;
  typedef Gt                                   Geom_traits;

  typedef CGAL::Delaunay_triangulation_2<Gt, TDS> Delaunay_triangulation_2;

  typedef typename CGAL::Periodic_2_offset_2   Offset;

  typedef typename Gt::FT                      FT;
  typedef typename Gt::Point_2                 Point;
  typedef typename Gt::Segment_2               Segment;
  typedef typename Gt::Vector_2                Vector;
  typedef typename Gt::Triangle_2              Triangle;

  typedef std::pair<Point, Offset>              Periodic_point;
  typedef array< std::pair<Point, Offset>, 2>   Periodic_segment;
  typedef array< std::pair<Point, Offset>, 3>   Periodic_triangle;
  typedef array< std::pair<Point, Offset>, 4>   Periodic_tetrahedron;

  typedef typename TDS::Face_handle            Face_handle;
  typedef typename TDS::Vertex_handle          Vertex_handle;
  typedef typename TDS::Edge                   Edge;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false                             Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_true                              Periodic_tag;

  typedef std::pair<Vector, Vector>             Basis;
  typedef cpp11::array<Vector, 3>               Voronoi_face_normals;

  template <class InputIterator>
  Periodic_2_Delaunay_triangulation_2_generic(InputIterator first, InputIterator beyond,
                                              const Basis& basis,
                                              const Gt& gt = Gt())
    : gt_(gt)
  {
    reduce_basis(basis);

    Voronoi_face_normals vfns = construct_Voronoi_face_normals(basis);

    std::vector<Point> cpoints = construct_canonical_points(first, beyond, vfns);

    construct_periodic_copies(cpoints, basis, vfns);

    dt2.insert(cpoints.begin(), cpoints.end());

    mark_canonical_simplices();
  }

#ifndef CGAL_CFG_NO_CPP0X_DELETED_AND_DEFAULT_FUNCTIONS
  Periodic_2_Delaunay_triangulation_2_generic& operator=(const Periodic_2_Delaunay_triangulation_2_generic&)=default;
#endif

  // @tmp
  Basis reduce_basis(const Basis& basis) { return basis; }

  // @tmp
  Voronoi_face_normals construct_Voronoi_face_normals(const Basis& basis)
  {
    Vector third = gt_.construct_opposite_vector_2_object()(
                       basis.first + basis.second);
    return CGAL::make_array(basis.first, basis.second, third);
  }

  ///
  Point construct_canonical_point(const Point& p,
                                  const Voronoi_face_normals& vfns)
  {
    Point cp = p;

    int vfn_pos = 0;
    while(vfn_pos < 3)
    {
      const Vector& vfn = vfns[vfn_pos];
      const Vector ptv(CGAL::ORIGIN, cp);

      const FT sp = gt_.compute_scalar_product_2_object()(ptv, vfn) /
                      gt_.compute_scalar_product_2_object()(vfn, vfn);

      if(-0.5 <= sp && sp < 0.5)
      {
        ++vfn_pos;
      }
      else
      {
        Vector tv = vfn;
        tv = gt_.construct_scaled_vector_2_object()(tv, - std::floor(sp + 0.5) );
        cp = gt_.construct_translated_point_2_object()(cp, tv);
        vfn_pos = 0;
      }
    }

    return cp;
  }

  template <class InputIterator>
  std::vector<Point> construct_canonical_points(InputIterator first, InputIterator beyond,
                                                const Voronoi_face_normals& vfns)
  {
    std::vector<Point> canonical_points;

    while(first != beyond)
    {
      const Point& p = *first;
      canonical_points.push_back(construct_canonical_point(p, vfns));
    }
  }

  template <class PointRange>
  void construct_periodic_copies(PointRange& pts,
                                 const Basis& basis,
                                 const Voronoi_face_normals& vfn)
  {

  }

  ///
  bool is_canonical(const Point& p,
                    const Basis& basis,
                    const Voronoi_face_normals& vfn) { }
  bool is_canonical(const Face_handle fh) { }


  void mark_canonical_simplices() { }

private:
  Delaunay_triangulation_2 dt2;
  Geom_traits gt_;
};

} //namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H
