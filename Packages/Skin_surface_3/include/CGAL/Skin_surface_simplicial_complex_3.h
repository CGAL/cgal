#ifndef SKIN_SURFACE_SIMPLICIAL_COMPLEX_3_H
#define SKIN_SURFACE_SIMPLICIAL_COMPLEX_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_incremental_builder_3.h>
// #include <CGAL/Skin_surface_quadratic_surface_3.h>

CGAL_BEGIN_NAMESPACE
template < class SkinSurfaceTraits3, 
	   class Cb = CGAL::Triangulation_vertex_base_3<
             typename SkinSurfaceTraits3::Simplicial_K> >
class Skin_surface_simplicial_vertex_base_3
  : public Cb
{
public:
  typedef typename SkinSurfaceTraits3::Simplicial_K Geom_traits;
  typedef typename SkinSurfaceTraits3::Mesh_K       Mesh_kernel;
  typedef typename Cb::Triangulation_data_structure Tds;
  typedef typename Tds::Vertex_handle               Vertex_handle;
  typedef typename Tds::Cell_handle                 Cell_handle;
  typedef typename Geom_traits::Point_3             Point;
  typedef typename Geom_traits::RT                  RT;
  typedef CGAL::Weighted_point<Point,RT>            Weighted_point;

  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
    typedef Skin_surface_simplicial_vertex_base_3<SkinSurfaceTraits3, Cb2>             Other;
  };

  Skin_surface_simplicial_vertex_base_3() : Cb(), iValue(0) { }
  Skin_surface_simplicial_vertex_base_3(Vertex_handle v0, Vertex_handle v1,
					Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3), iValue(RT(0)) {}

  typename Mesh_kernel::RT iValue;
};


template < class SkinSurfaceTraits3, 
	   class Cb = CGAL::Triangulation_cell_base_3<typename SkinSurfaceTraits3::Simplicial_K> >
class Skin_surface_simplicial_cell_base_3
  : public Cb
{
public:
  typedef typename SkinSurfaceTraits3::Simplicial_K Geom_traits;
  typedef typename Cb::Triangulation_data_structure Tds;
  typedef typename Tds::Vertex_handle               Vertex_handle;
  typedef typename Tds::Cell_handle                 Cell_handle;
  typedef typename Geom_traits::Point_3             Point;
  typedef typename Geom_traits::Tetrahedron_3       Tetrahedron;
  typedef typename Geom_traits::Segment_3           Segment;
  typedef typename Geom_traits::Line_3              Line;
  typedef typename Geom_traits::Plane_3             Plane;
  typedef typename Geom_traits::Triangle_3          Triangle;
  typedef typename Geom_traits::RT                  RT;
  typedef CGAL::Weighted_point<Point,RT>            Weighted_point;

//   typedef Skin_surface_quadratic_surface_3<typename SkinSurfaceTraits3::Mesh_K>               QuadrSurface;
	
  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
    typedef Skin_surface_simplicial_cell_base_3<SkinSurfaceTraits3, Cb2>               Other;
  };

  Skin_surface_simplicial_cell_base_3() : Cb() {
  }
  Skin_surface_simplicial_cell_base_3(Vertex_handle v0, Vertex_handle v1,
				      Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {
  }

//   QuadrSurface *surf;
};
	
template < class GT,
           class Tds = Triangulation_data_structure_3 <
  Skin_surface_simplicial_vertex_base_3<GT>,
  Skin_surface_simplicial_cell_base_3<GT> > >
class Skin_surface_simplicial_complex_3
  : public Triangulation_3<GT,Tds> {
public:
  typedef Skin_surface_simplicial_complex_3<GT, Tds>             Self;
  typedef Triangulation_3<GT, Tds>                Parent;
  typedef typename Parent::Vertex_handle          Vertex_handle;
	
  friend class Triangulation_incremental_builder_3<Self>;
};


CGAL_END_NAMESPACE

#endif // SKIN_SURFACE_SIMPLICIAL_COMPLEX_3_H
