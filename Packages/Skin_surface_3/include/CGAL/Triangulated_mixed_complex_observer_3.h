#ifndef TRIANGULATED_MIXED_COMPLEX_OBSERVER_3
#define TRIANGULATED_MIXED_COMPLEX_OBSERVER_3

#include <CGAL/Triangulation_simplex_3.h>
#include <CGAL/Skin_surface_quadratic_surface_3.h>

CGAL_BEGIN_NAMESPACE

template <class SkinSurfaceTraits_3,
	  class TriangulatedMixedComplex_3,
	  class RegularTriangulation_3>
class Triangulated_mixed_complex_observer_3 {
public:
  typedef SkinSurfaceTraits_3                        Skin_traits_3;
  typedef typename Skin_traits_3::Regular_traits     Regular_traits;
  typedef typename Skin_traits_3::Triangulated_mixed_complex_kernel
                                              Triangulated_mixed_complex_kernel;
  typedef typename Skin_traits_3::Polyhedron_kernel  Polyhedron_kernel;

  typedef RegularTriangulation_3                     Regular;
  typedef TriangulatedMixedComplex_3                 Triangulated_mixed_complex;

  typedef typename Regular::Vertex_handle            Rt_Vertex_handle;
  typedef typename Regular::Edge                     Rt_Edge;
  typedef typename Regular::Facet                    Rt_Facet;
  typedef typename Regular::Cell_handle              Rt_Cell_handle;
  typedef Triangulation_simplex_3<Regular>           Rt_Simplex;

  typedef typename Regular::Bare_point               Rt_Point;
  typedef typename Regular::Geom_traits              Rt_Geom_traits;
  typedef typename Rt_Geom_traits::RT                Rt_RT;
  typedef typename Regular::Weighted_point           Rt_Weighted_point;

  typedef typename Triangulated_mixed_complex::Vertex_handle Sc_Vertex_handle;
  typedef typename Triangulated_mixed_complex::Edge          Sc_Edge;
  typedef typename Triangulated_mixed_complex::Facet         Sc_Facet;
  typedef typename Triangulated_mixed_complex::Cell_handle   Sc_Cell_handle;
  
  typedef typename Skin_traits_3::Triangulated_mixed_complex_kernel  Triangulated_mixed_complex_kernel;
  typedef typename Triangulated_mixed_complex_kernel::Point_3        Sc_Point;
  typedef typename Triangulated_mixed_complex_kernel::RT             Sc_RT;
  typedef Weighted_point<Sc_Point,Sc_RT>        Sc_Weighted_point;
  typedef Skin_surface_quadratic_surface_3<Polyhedron_kernel> QuadrSurface;
  typedef Skin_surface_sphere_3<Polyhedron_kernel>         Sphere_surface;
  typedef Skin_surface_hyperboloid_3<Polyhedron_kernel>    Hyperboloid_surface;

  typedef typename Polyhedron_kernel::RT                   Mesh_RT;
  typedef typename Polyhedron_kernel::Point_3              Mesh_Point;
  typedef Weighted_point<Mesh_Point,Mesh_RT>    Mesh_Weighted_point;
  
  typedef typename Skin_traits_3::R2P_converter R2P_converter;
  typedef typename Skin_traits_3::T2P_converter T2P_converter;

  Triangulated_mixed_complex_observer_3(Mesh_RT shrink) : 
    shrink(shrink), r2p_converter(SkinSurfaceTraits_3().r2p_converter_object()),
    t2p_converter(SkinSurfaceTraits_3().t2p_converter_object()) {
  }

  void after_vertex_insertion(
    Rt_Simplex &sDel, Rt_Simplex &sVor, Sc_Vertex_handle &vh) 
  {
  }

  void after_cell_insertion(Rt_Simplex &s, Sc_Cell_handle &ch) 
  {
    if (!(s == prev_s)) {
      prev_s = s;
      Rt_Vertex_handle vh;
      Rt_Edge          e;
      Rt_Facet         f;
      Rt_Cell_handle   ch;

      switch (s.dimension()) {
        case 0:
	  vh = s;
	  surf = new Sphere_surface(r2p_converter(vh->point()), shrink, 1);
	  break;
        case 1:
	  e = s;
	  surf = new Hyperboloid_surface(
	    r2p_converter(Rt_Weighted_point(
	      Rt_Geom_traits().construct_weighted_circumcenter_3_object()(
		e.first->vertex(e.second)->point(),
		e.first->vertex(e.third)->point()),
	      Rt_Geom_traits().
		compute_squared_radius_smallest_orthogonal_sphere_3_object()(
		  e.first->vertex(e.second)->point(),
		  e.first->vertex(e.third)->point()))),
	    r2p_converter(
	      e.first->vertex(e.second)->point()-
	      e.first->vertex(e.third)->point()),
	    shrink, 1);
	  break;
        case 2:
	  f = s;
	  surf = new Hyperboloid_surface(
	    r2p_converter(Rt_Weighted_point(
	      Rt_Geom_traits().construct_weighted_circumcenter_3_object()(
		f.first->vertex((f.second+1)&3)->point(),
		f.first->vertex((f.second+2)&3)->point(),
		f.first->vertex((f.second+3)&3)->point()),
	      Rt_Geom_traits().
		compute_squared_radius_smallest_orthogonal_sphere_3_object()(
		  f.first->vertex((f.second+1)&3)->point(),
		  f.first->vertex((f.second+2)&3)->point(),
		  f.first->vertex((f.second+3)&3)->point()))),
	      typename Polyhedron_kernel::Construct_orthogonal_vector_3()(
		r2p_converter(f.first->vertex((f.second+1)&3)->point()),
		r2p_converter(f.first->vertex((f.second+2)&3)->point()),
		r2p_converter(f.first->vertex((f.second+3)&3)->point())),
	    1-shrink, -1);
	  break;
	case 3:
	  ch = s;
	  surf = new Sphere_surface(
	    r2p_converter(Rt_Weighted_point(
	      Rt_Geom_traits().construct_weighted_circumcenter_3_object()(
		ch->vertex(0)->point(),
		ch->vertex(1)->point(),
		ch->vertex(2)->point(),
		ch->vertex(3)->point()),
	      Rt_Geom_traits().
		compute_squared_radius_smallest_orthogonal_sphere_3_object()(
		  ch->vertex(0)->point(),
		  ch->vertex(1)->point(),
		  ch->vertex(2)->point(),
		  ch->vertex(3)->point()))),
	      1-shrink, -1);
      }

    }
    ch->surf = surf;
  }

  Mesh_RT shrink;
  Rt_Simplex prev_s;
  QuadrSurface *surf;
  R2P_converter r2p_converter;
  T2P_converter t2p_converter;
};

CGAL_END_NAMESPACE

#endif // TRIANGULATED_MIXED_COMPLEX_OBSERVER_3
