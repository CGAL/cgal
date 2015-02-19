#define CGAL_NO_DEPRECATION_WARNINGS

#include <CGAL/basic.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Delaunay_d.h>
#include <CGAL/predicates_d.h>
#include <CGAL/constructions_d.h>
#include <CGAL/double.h>
#include <iostream>
#include <CGAL/test_macros.h>
#include <CGAL/use.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_real.h>
typedef leda_integer RT;
typedef leda_real FT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
typedef double FT;
#else
typedef double RT;
typedef double FT;
#endif
#endif

int main()
{
  CGAL::set_pretty_mode ( std::cerr );
  CGAL_KD_SETDTHREAD(193);
  CGAL_TEST_START;
  {
  typedef CGAL::Cartesian_d<FT> Kernel;
  typedef CGAL::Delaunay_d<Kernel> Delaunay_d;
  typedef Delaunay_d::Point_d Point_d;
  typedef Delaunay_d::Lifted_hyperplane_d Hyperplane_d;
  typedef Delaunay_d::Sphere_d Sphere_d;
  typedef Delaunay_d::Simplex_handle Simplex_handle;
  typedef Delaunay_d::Vertex_handle Vertex_handle;
  typedef Delaunay_d::Facet_handle Facet_handle;

  {
    Delaunay_d DT(2);
    CGAL_TEST(DT.empty());
    Point_d p1(0,0,1);
    Point_d p2(1,0,1);
    Point_d p3(0,1,1);
    Point_d p4(1,1,1);
    Point_d p5(1,1,3);
    Point_d p6(3,2,1);
    Point_d p7(3,3,4);
    Vertex_handle v1,v2,v3,v4,v5;
    CGAL_TEST(DT.dimension()==2);
    v1 = DT.insert(p1);
    v2 = DT.insert(p2);
    CGAL_TEST(DT.current_dimension()==1);
    CGAL_TEST(DT.associated_point(v2)==p2);
    Simplex_handle ds = DT.simplex(v2);
    int di = DT.index(v2);
    CGAL_TEST(DT.vertex_of_simplex(ds,di)==v2);
    CGAL_TEST(DT.is_simplex_of_nearest(ds) && DT.is_simplex_of_furthest(ds));
    v3 = DT.insert(p3);
    CGAL_TEST(DT.locate(p4)==Simplex_handle());
    v4 = DT.insert(p4);
    v5 = DT.insert(p5);
    Simplex_handle ds2 = DT.simplex(v5);
    int di2 = DT.index(v5);
    CGAL_TEST(DT.is_simplex_of_nearest(ds2)
              && ! DT.is_simplex_of_furthest(ds2));
    CGAL_TEST(DT.point_of_simplex(ds2,di2) == DT.associated_point(v5));
    CGAL_TEST(DT.opposite_simplex(ds2,1)!=Simplex_handle());
    CGAL_TEST(DT.opposite_simplex(DT.opposite_simplex(ds2,1),
      DT.index_of_vertex_in_opposite_simplex(ds2,1)) == ds2);

    std::list<Simplex_handle> NL = DT.all_simplices(Delaunay_d::NEAREST);
    std::list<Simplex_handle> FL = DT.all_simplices(Delaunay_d::FURTHEST);
    std::list<Facet_handle> Facets = DT.all_facets();
    std::list<Facet_handle>::iterator fit;
    int n_num(0), f_num(0);
    for(fit = Facets.begin(); fit != Facets.end(); ++fit) {
      Hyperplane_d h = DT.base_facet_plane(*fit);
      if ( h.orthogonal_vector().homogeneous(h.dimension()-1)>0 ) ++f_num;
      if ( h.orthogonal_vector().homogeneous(h.dimension()-1)<0 ) ++n_num;
    }
    CGAL_TEST(NL.size()==unsigned(n_num));
    CGAL_TEST(FL.size()==unsigned(f_num));
    Point_d q1 = DT.point_of_simplex(ds,0);
    Point_d q2 = DT.point_of_simplex(ds,1);
    CGAL_TEST(DT.contains(ds,CGAL::midpoint(q1,q2)));
    CGAL_TEST(DT.locate(p6)==Simplex_handle());
    CGAL_TEST(DT.locate(p7)!=Simplex_handle());
    CGAL_TEST(DT.lookup(p6)==Vertex_handle());
    CGAL_TEST(DT.associated_point(DT.lookup(p4))==p4);
    CGAL_TEST(DT.associated_point(DT.nearest_neighbor(p6))==p4);
    std::vector<Point_d> V = make_vector(p1,p2,p5);
    Sphere_d K(2,V.begin(),V.end());
    std::list<Vertex_handle> RL = DT.range_search(K);
    std::list<Vertex_handle>::iterator it = RL.begin();
    CGAL_TEST(*it++==v1 && *it++==v2 && *it++==v5);
    RL = DT.range_search(make_vector(p1,p2,p3));
    CGAL_TEST(RL.size()==4);
    std::list<Point_d> L = DT.all_points();
    CGAL_TEST(*L.begin()==p1);
    // show and graphrep in demos
    DT.clear();
    DT.insert(p5);
    CGAL_TEST(DT.current_dimension()==0);
  }

  }

  {
  typedef CGAL::Homogeneous_d<RT> Kernel;
  typedef CGAL::Delaunay_d<Kernel> Delaunay_d;
  typedef Delaunay_d::Point_d Point_d;
  typedef Delaunay_d::Lifted_hyperplane_d Hyperplane_d;
  typedef Delaunay_d::Sphere_d Sphere_d;
  typedef Delaunay_d::Simplex_handle Simplex_handle;
  typedef Delaunay_d::Vertex_handle Vertex_handle;
  typedef Delaunay_d::Facet_handle Facet_handle;

  {
    Delaunay_d DT(2);
    CGAL_TEST(DT.empty());
    Point_d p1(0,0,1);
    Point_d p2(1,0,1);
    Point_d p3(0,1,1);
    Point_d p4(1,1,1);
    Point_d p5(1,1,3);
    Point_d p6(3,2,1);
    Point_d p7(3,3,4);
    Vertex_handle v1,v2,v3,v4,v5;
    CGAL_TEST(DT.dimension()==2);
    v1 = DT.insert(p1);
    v2 = DT.insert(p2);
    CGAL_TEST(DT.current_dimension()==1);
    CGAL_TEST(DT.associated_point(v2)==p2);
    Simplex_handle ds = DT.simplex(v2);
    int di = DT.index(v2);
    CGAL_TEST(DT.vertex_of_simplex(ds,di)==v2);
    CGAL_TEST(DT.is_simplex_of_nearest(ds) && DT.is_simplex_of_furthest(ds));
    v3 = DT.insert(p3);
    CGAL_TEST(DT.locate(p4)==Simplex_handle());
    v4 = DT.insert(p4);
    v5 = DT.insert(p5);
    Simplex_handle ds2 = DT.simplex(v5);
    int di2 = DT.index(v5);
    CGAL_TEST(DT.is_simplex_of_nearest(ds2)
              && ! DT.is_simplex_of_furthest(ds2));
    CGAL_TEST(DT.point_of_simplex(ds2,di2) == DT.associated_point(v5));
    CGAL_TEST(DT.opposite_simplex(ds2,1)!=Simplex_handle());
    CGAL_TEST(DT.opposite_simplex(DT.opposite_simplex(ds2,1),
      DT.index_of_vertex_in_opposite_simplex(ds2,1)) == ds2);

    std::list<Simplex_handle> NL = DT.all_simplices(Delaunay_d::NEAREST);
    std::list<Simplex_handle> FL = DT.all_simplices(Delaunay_d::FURTHEST);
    std::list<Facet_handle> Facets = DT.all_facets();
    std::list<Facet_handle>::iterator fit;
    int n_num(0), f_num(0);
    for(fit = Facets.begin(); fit != Facets.end(); ++fit) {
      Hyperplane_d h = DT.base_facet_plane(*fit);
      if ( h.orthogonal_vector().homogeneous(h.dimension()-1)>0 ) ++f_num;
      if ( h.orthogonal_vector().homogeneous(h.dimension()-1)<0 ) ++n_num;
    }
    CGAL_TEST(NL.size()==unsigned(n_num));
    CGAL_TEST(FL.size()==unsigned(f_num));
    Point_d q1 = DT.point_of_simplex(ds,0);
    Point_d q2 = DT.point_of_simplex(ds,1);
    CGAL_TEST(DT.contains(ds,CGAL::midpoint(q1,q2)));
    CGAL_TEST(DT.locate(p6)==Simplex_handle());
    CGAL_TEST(DT.locate(p7)!=Simplex_handle());
    CGAL_TEST(DT.lookup(p6)==Vertex_handle());
    CGAL_TEST(DT.associated_point(DT.lookup(p4))==p4);
    CGAL_TEST(DT.associated_point(DT.nearest_neighbor(p6))==p4);
    std::vector<Point_d> V = make_vector(p1,p2,p5);
    Sphere_d K(2,V.begin(),V.end());
    std::list<Vertex_handle> RL = DT.range_search(K);
    std::list<Vertex_handle>::iterator it = RL.begin();
    CGAL_TEST(*it++==v1 && *it++==v2 && *it++==v5);
    RL = DT.range_search(make_vector(p1,p2,p3));
    CGAL_TEST(RL.size()==4);
    std::list<Point_d> L = DT.all_points();
    CGAL_TEST(*L.begin()==p1);
    // show and graphrep in demos
    DT.clear();
    DT.insert(p5);
    CGAL_TEST(DT.current_dimension()==0);
  }

  }

  CGAL_TEST_END;
}

