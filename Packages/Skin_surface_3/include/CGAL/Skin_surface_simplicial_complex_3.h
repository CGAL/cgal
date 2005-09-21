#ifndef SKIN_SURFACE_SIMPLICIAL_COMPLEX_3_H
#define SKIN_SURFACE_SIMPLICIAL_COMPLEX_3_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_incremental_builder_3.h>
#include <CGAL/Skin_surface_quadratic_surface_3.h>

#include <CGAL/Simplex_3.h>

CGAL_BEGIN_NAMESPACE

template < class Simplicial_K_,
           class Regular,
           class Vb = CGAL::Triangulation_vertex_base_3<Simplicial_K_> >
class Skin_surface_simplicial_vertex_base_3
  : public Vb
{
public:
  typedef Simplicial_K_                    Geom_traits;
  typedef typename Vb::Triangulation_data_structure Tds;
  typedef typename Tds::Vertex_handle               Vertex_handle;
  typedef typename Tds::Cell_handle                 Cell_handle;
  typedef typename Geom_traits::Point_3             Point;
  typedef typename Geom_traits::RT                  RT;
  typedef CGAL::Weighted_point<Point,RT>            Weighted_point;

  typedef Simplex_3<Regular>                        Rt_Simplex;


  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef Skin_surface_simplicial_vertex_base_3<Simplicial_K_, Regular, Vb2>
                                                           Other;
  };

  Skin_surface_simplicial_vertex_base_3() : Vb() { }
  Skin_surface_simplicial_vertex_base_3(Vertex_handle v0, Vertex_handle v1,
                                        Vertex_handle v2, Vertex_handle v3)
    : Vb(v0, v1, v2, v3) {}

};

template < class SimplicialTraits_3, 
           class MeshTraits_3, 
	   class Regular,
	   class Cb = Triangulation_cell_base_3<SimplicialTraits_3> >
class Skin_surface_simplicial_cell_base_3
  : public Cb
{
public:
  typedef SimplicialTraits_3                        Geom_traits;
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

  typedef MeshTraits_3                                    Mesh_traits_3;
  typedef Skin_surface_quadratic_surface_3<Mesh_traits_3> QuadrSurface;
	
  typedef Simplex_3<Regular>                              Rt_Simplex;
  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other  Cb2;
    typedef Skin_surface_simplicial_cell_base_3<
      SimplicialTraits_3, Mesh_traits_3, Regular, Cb2>              Other;
  };

  Skin_surface_simplicial_cell_base_3() : Cb() {
  }
  Skin_surface_simplicial_cell_base_3(Vertex_handle v0, Vertex_handle v1,
				      Vertex_handle v2, Vertex_handle v3)
    : Cb(v0, v1, v2, v3) {
  }

  QuadrSurface *surf;
};
	
template < class GT,
           class Tds = Triangulation_data_structure_3 <
             Triangulation_vertex_base_3<GT>,
             Triangulation_cell_base_3<GT> > >
class Skin_surface_simplicial_complex_3
  : public Triangulation_3<GT,Tds> {
public:
  typedef Skin_surface_simplicial_complex_3<GT, Tds> Self;
  typedef Triangulation_3<GT, Tds>                   Parent;
  typedef typename Parent::Vertex_handle             Vertex_handle;
	
  friend class Triangulation_incremental_builder_3<Self>;
};


CGAL_END_NAMESPACE

#endif // SKIN_SURFACE_SIMPLICIAL_COMPLEX_3_H
