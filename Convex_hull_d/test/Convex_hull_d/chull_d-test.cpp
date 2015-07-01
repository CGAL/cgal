#define CGAL_NO_DEPRECATION_WARNINGS

#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/test_macros.h>
#include <iostream>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer RT;
#else
#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz RT;
#else
#include <CGAL/double.h>
typedef double RT;
#endif
#endif

int main()
{
  CGAL_KD_SETDTHREAD(11);  
  CGAL::set_pretty_mode ( std::cerr );
  CGAL_TEST_START;
  {
    typedef CGAL::Homogeneous_d<RT> Kernel;
    typedef CGAL::Convex_hull_d<Kernel> Convex_hull_d;
    typedef Convex_hull_d::Point_d Point;
    typedef Convex_hull_d::Hyperplane_d Plane;
    typedef Convex_hull_d::Simplex_handle Simplex_handle;
    typedef Convex_hull_d::Facet_handle Facet_handle;
    typedef Convex_hull_d::Vertex_handle Vertex_handle;

    Convex_hull_d T1(2);  
    const Convex_hull_d* pT1 = &T1;
    Point p1(0,0,1);
    Point p2(1,0,1);
    Point p3(0,1,1);
    Point p4(1,1,1);
    CGAL_TEST(T1.dimension()==2);
    CGAL_TEST(T1.current_dimension()==-1);
    Vertex_handle v1 = T1.insert(p1);
    CGAL_TEST(T1.associated_point(v1)==p1);
    T1.insert(p2);
    CGAL_TEST(T1.current_dimension()==1);
    CGAL_TEST(T1.is_dimension_jump(p3));
    T1.insert(p3);
    Simplex_handle s1 = T1.simplex(v1);
    int i1 = T1.index(v1);
    CGAL_TEST(T1.vertex_of_simplex(s1,i1)==v1);
    CGAL_TEST(T1.associated_point(T1.vertex_of_simplex(s1,1))==
              T1.point_of_simplex(s1,1));
    T1.insert(p3);
    CGAL_TEST((T1.opposite_simplex(T1.opposite_simplex(s1,1),
              T1.index_of_vertex_in_opposite_simplex(s1,1)) == s1));
    std::list<Simplex_handle> L = T1.all_simplices();
    CGAL_TEST(L.size() == 4);
    std::list<Facet_handle> F = T1.all_facets();
    CGAL_TEST(F.size() == 3);
    Facet_handle f1 = *(L.begin());
    CGAL_TEST(T1.associated_point(T1.vertex_of_facet(f1,1))==
              T1.point_of_facet(f1,1));
    CGAL_TEST((T1.opposite_facet(T1.opposite_facet(f1,1),
               T1.index_of_vertex_in_opposite_facet(f1,1)) == f1));
    Point pf0 = T1.point_of_facet(f1,0);
    Point pf1 = T1.point_of_facet(f1,1);
    Plane h = T1.hyperplane_supporting(f1);
    CGAL_TEST(h.has_on(pf0)&&h.has_on(pf1));
    std::list<Facet_handle> G = T1.facets_visible_from(p4);
    CGAL_TEST(G.size()==1);
    CGAL_TEST(T1.bounded_side(p4)==CGAL::ON_UNBOUNDED_SIDE && 
              T1.bounded_side(Point(1,1,10))==CGAL::ON_BOUNDED_SIDE);

    Convex_hull_d::Point_const_iterator   pit;
    Convex_hull_d::Vertex_iterator  vit;
    Convex_hull_d::Simplex_iterator sit;
    for (pit = T1.points_begin(); pit != T1.points_end(); pit++) *pit;
    for (vit = T1.vertices_begin(); vit != T1.vertices_end(); vit++) *vit;
    for (sit = T1.simplices_begin(); sit != T1.simplices_end(); sit++) *sit;
    T1.is_valid();
    T1.clear(2);
    CGAL_TEST(T1.number_of_vertices()==0);
    CGAL_TEST(T1.number_of_facets()==0);
    CGAL_TEST(T1.number_of_simplices()==0);
    std::vector<Point> V = make_vector(p1,p2,p3,p4);
    T1.initialize(V.begin(),V.end());
    Convex_hull_d::Facet_iterator fit;
    int fnum(0);
    for (fit = T1.facets_begin(); fit != T1.facets_end(); ++fnum, ++fit) *fit;
    CGAL_TEST(fnum==4);

#ifndef _MSC_VER // truncation due to name length exceeded
    Convex_hull_d::Hull_vertex_iterator hvit;
    int vnum(0);
    for (hvit = T1.hull_vertices_begin(); 
         hvit != T1.hull_vertices_end(); ++vnum, ++hvit) *hvit; 
    CGAL_TEST(vnum==4);

    Convex_hull_d::Facet_const_iterator fcit; 
    for (fcit = pT1->facets_begin(); fcit != pT1->facets_end(); ++fcit) *fcit;
    Convex_hull_d::Hull_vertex_const_iterator hvcit; 
    for (hvcit = pT1->hull_vertices_begin(); 
         hvcit != pT1->hull_vertices_end(); ++hvcit) *hvcit;
    Convex_hull_d::Hull_point_const_iterator hpcit; 
    for (hpcit = pT1->hull_points_begin(); 
         hpcit != pT1->hull_points_end(); ++hpcit) *hpcit;
#endif
  }


  {
    typedef CGAL::Cartesian<double> Kernel_3;
    typedef CGAL::Convex_hull_d_traits_3<Kernel_3> Kernel;
    typedef CGAL::Convex_hull_d<Kernel> Convex_hull_d;
    typedef Convex_hull_d::Point_d Point;
    typedef Convex_hull_d::Hyperplane_d Plane;
    typedef Convex_hull_d::Simplex_handle Simplex_handle;
    typedef Convex_hull_d::Facet_handle Facet_handle;
    typedef Convex_hull_d::Vertex_handle Vertex_handle;
    Convex_hull_d T2(3);  
    Point p1(0,0,0,1);
    Point p2(1,0,0,1);
    Point p3(0,1,0,1);
    Point p4(0,0,1,1);
    Point p5(1,1,1,1);
    CGAL_TEST(T2.dimension()==3);
    CGAL_TEST(T2.current_dimension()==-1);
    T2.insert(p1);
    T2.insert(p2);
    CGAL_TEST(T2.is_dimension_jump(p3));
    T2.insert(p3);
    CGAL_TEST(T2.current_dimension()==2);
    CGAL_TEST(T2.is_dimension_jump(p4));
    Vertex_handle v1 = T2.insert(p4);
    CGAL_TEST(T2.associated_point(v1)==p4);

    Simplex_handle s1 = T2.simplex(v1);
    int i1 = T2.index(v1);


    CGAL_TEST(T2.vertex_of_simplex(s1,i1)==v1);
    CGAL_TEST(T2.associated_point(T2.vertex_of_simplex(s1,1))==
          T2.point_of_simplex(s1,1));
    CGAL_TEST((T2.opposite_simplex(T2.opposite_simplex(s1,1),
               T2.index_of_vertex_in_opposite_simplex(s1,1)) == s1));

    std::list<Simplex_handle> L = T2.all_simplices();
    CGAL_TEST(L.size() == 5);
    std::list<Facet_handle> F = T2.all_facets();
    CGAL_TEST(F.size() == 4);
    Facet_handle f1 = *(L.begin());
    CGAL_TEST(T2.associated_point(T2.vertex_of_facet(f1,1))==
              T2.point_of_facet(f1,1));
    CGAL_TEST((T2.opposite_facet(T2.opposite_facet(f1,1),
               T2.index_of_vertex_in_opposite_facet(f1,1)) == f1));

    Point pf0 = T2.point_of_facet(f1,0);
    Point pf1 = T2.point_of_facet(f1,1);
    Point pf2 = T2.point_of_facet(f1,2);
    Plane h = T2.hyperplane_supporting(f1);
    CGAL_TEST(h.has_on(pf0)&&h.has_on(pf1)&&h.has_on(pf2));

    std::list<Facet_handle> G = T2.facets_visible_from(p5);
    CGAL_TEST(G.size()==1);
    CGAL_TEST(T2.bounded_side(p5)==CGAL::ON_UNBOUNDED_SIDE && 
              T2.bounded_side(Point(1,1,1,10))==CGAL::ON_BOUNDED_SIDE);
    T2.is_valid();
    T2.clear(3);


  }


  {
    typedef CGAL::Homogeneous<RT> Kernel_3;
    typedef CGAL::Convex_hull_d_traits_3<Kernel_3> Kernel;
    typedef CGAL::Convex_hull_d<Kernel> Convex_hull_d;
    typedef Convex_hull_d::Point_d Point;
    typedef Convex_hull_d::Hyperplane_d Plane;
    typedef Convex_hull_d::Simplex_handle Simplex_handle;
    typedef Convex_hull_d::Facet_handle Facet_handle;
    typedef Convex_hull_d::Vertex_handle Vertex_handle;
    Convex_hull_d T2(3);  
    Point p1(0,0,0,1);
    Point p2(1,0,0,1);
    Point p3(0,1,0,1);
    Point p4(0,0,1,1);
    Point p5(1,1,1,1);
    CGAL_TEST(T2.dimension()==3);
    CGAL_TEST(T2.current_dimension()==-1);
    T2.insert(p1);
    T2.insert(p2);
    CGAL_TEST(T2.is_dimension_jump(p3));
    T2.insert(p3);
    CGAL_TEST(T2.current_dimension()==2);
    CGAL_TEST(T2.is_dimension_jump(p4));
    Vertex_handle v1 = T2.insert(p4);
    CGAL_TEST(T2.associated_point(v1)==p4);

    Simplex_handle s1 = T2.simplex(v1);
    int i1 = T2.index(v1);


    CGAL_TEST(T2.vertex_of_simplex(s1,i1)==v1);
    CGAL_TEST(T2.associated_point(T2.vertex_of_simplex(s1,1))==
          T2.point_of_simplex(s1,1));
    CGAL_TEST((T2.opposite_simplex(T2.opposite_simplex(s1,1),
               T2.index_of_vertex_in_opposite_simplex(s1,1)) == s1));

    std::list<Simplex_handle> L = T2.all_simplices();
    CGAL_TEST(L.size() == 5);
    std::list<Facet_handle> F = T2.all_facets();
    CGAL_TEST(F.size() == 4);
    Facet_handle f1 = *(L.begin());
    CGAL_TEST(T2.associated_point(T2.vertex_of_facet(f1,1))==
              T2.point_of_facet(f1,1));
    CGAL_TEST((T2.opposite_facet(T2.opposite_facet(f1,1),
               T2.index_of_vertex_in_opposite_facet(f1,1)) == f1));

    Point pf0 = T2.point_of_facet(f1,0);
    Point pf1 = T2.point_of_facet(f1,1);
    Point pf2 = T2.point_of_facet(f1,2);
    Plane h = T2.hyperplane_supporting(f1);
    CGAL_TEST(h.has_on(pf0)&&h.has_on(pf1)&&h.has_on(pf2));

    std::list<Facet_handle> G = T2.facets_visible_from(p5);
    CGAL_TEST(G.size()==1);
    CGAL_TEST(T2.bounded_side(p5)==CGAL::ON_UNBOUNDED_SIDE && 
              T2.bounded_side(Point(1,1,1,10))==CGAL::ON_BOUNDED_SIDE);
    T2.is_valid();
    T2.clear(3);


  }
  {
    typedef CGAL::Cartesian_d<double> Kernel;
    typedef CGAL::Convex_hull_d<Kernel> Convex_hull_d;
    typedef Convex_hull_d::Point_d Point;
    typedef Convex_hull_d::Hyperplane_d Plane;
    typedef Convex_hull_d::Simplex_handle Simplex_handle;
    typedef Convex_hull_d::Facet_handle Facet_handle;
    typedef Convex_hull_d::Vertex_handle Vertex_handle;
    Convex_hull_d T2(3);  
    Point p1(0,0,0,1);
    Point p2(1,0,0,1);
    Point p3(0,1,0,1);
    Point p4(0,0,1,1);
    Point p5(1,1,1,1);
    CGAL_TEST(T2.dimension()==3);
    CGAL_TEST(T2.current_dimension()==-1);
    T2.insert(p1);
    T2.insert(p2);
    CGAL_TEST(T2.is_dimension_jump(p3));
    T2.insert(p3);
    CGAL_TEST(T2.current_dimension()==2);
    CGAL_TEST(T2.is_dimension_jump(p4));
    Vertex_handle v1 = T2.insert(p4);
    CGAL_TEST(T2.associated_point(v1)==p4);

    Simplex_handle s1 = T2.simplex(v1);
    int i1 = T2.index(v1);


    CGAL_TEST(T2.vertex_of_simplex(s1,i1)==v1);
    CGAL_TEST(T2.associated_point(T2.vertex_of_simplex(s1,1))==
          T2.point_of_simplex(s1,1));
    CGAL_TEST((T2.opposite_simplex(T2.opposite_simplex(s1,1),
               T2.index_of_vertex_in_opposite_simplex(s1,1)) == s1));

    std::list<Simplex_handle> L = T2.all_simplices();
    CGAL_TEST(L.size() == 5);
    std::list<Facet_handle> F = T2.all_facets();
    CGAL_TEST(F.size() == 4);
    Facet_handle f1 = *(L.begin());
    CGAL_TEST(T2.associated_point(T2.vertex_of_facet(f1,1))==
              T2.point_of_facet(f1,1));
    CGAL_TEST((T2.opposite_facet(T2.opposite_facet(f1,1),
               T2.index_of_vertex_in_opposite_facet(f1,1)) == f1));

    Point pf0 = T2.point_of_facet(f1,0);
    Point pf1 = T2.point_of_facet(f1,1);
    Point pf2 = T2.point_of_facet(f1,2);
    Plane h = T2.hyperplane_supporting(f1);
    CGAL_TEST(h.has_on(pf0)&&h.has_on(pf1)&&h.has_on(pf2));

    std::list<Facet_handle> G = T2.facets_visible_from(p5);
    CGAL_TEST(G.size()==1);
    CGAL_TEST(T2.bounded_side(p5)==CGAL::ON_UNBOUNDED_SIDE && 
              T2.bounded_side(Point(1,1,1,10))==CGAL::ON_BOUNDED_SIDE);
    T2.is_valid();
    T2.clear(3);


  }
  {
    typedef CGAL::Homogeneous_d<RT> Kernel;
    typedef CGAL::Convex_hull_d<Kernel> Convex_hull_d;
    typedef Convex_hull_d::Point_d Point;
    typedef Convex_hull_d::Hyperplane_d Plane;
    typedef Convex_hull_d::Simplex_handle Simplex_handle;
    typedef Convex_hull_d::Facet_handle Facet_handle;
    typedef Convex_hull_d::Vertex_handle Vertex_handle;
    Convex_hull_d T2(3);  
    Point p1(0,0,0,1);
    Point p2(1,0,0,1);
    Point p3(0,1,0,1);
    Point p4(0,0,1,1);
    Point p5(1,1,1,1);
    CGAL_TEST(T2.dimension()==3);
    CGAL_TEST(T2.current_dimension()==-1);
    T2.insert(p1);
    T2.insert(p2);
    CGAL_TEST(T2.is_dimension_jump(p3));
    T2.insert(p3);
    CGAL_TEST(T2.current_dimension()==2);
    CGAL_TEST(T2.is_dimension_jump(p4));
    Vertex_handle v1 = T2.insert(p4);
    CGAL_TEST(T2.associated_point(v1)==p4);

    Simplex_handle s1 = T2.simplex(v1);
    int i1 = T2.index(v1);


    CGAL_TEST(T2.vertex_of_simplex(s1,i1)==v1);
    CGAL_TEST(T2.associated_point(T2.vertex_of_simplex(s1,1))==
          T2.point_of_simplex(s1,1));
    CGAL_TEST((T2.opposite_simplex(T2.opposite_simplex(s1,1),
               T2.index_of_vertex_in_opposite_simplex(s1,1)) == s1));

    std::list<Simplex_handle> L = T2.all_simplices();
    CGAL_TEST(L.size() == 5);
    std::list<Facet_handle> F = T2.all_facets();
    CGAL_TEST(F.size() == 4);
    Facet_handle f1 = *(L.begin());
    CGAL_TEST(T2.associated_point(T2.vertex_of_facet(f1,1))==
              T2.point_of_facet(f1,1));
    CGAL_TEST((T2.opposite_facet(T2.opposite_facet(f1,1),
               T2.index_of_vertex_in_opposite_facet(f1,1)) == f1));

    Point pf0 = T2.point_of_facet(f1,0);
    Point pf1 = T2.point_of_facet(f1,1);
    Point pf2 = T2.point_of_facet(f1,2);
    Plane h = T2.hyperplane_supporting(f1);
    CGAL_TEST(h.has_on(pf0)&&h.has_on(pf1)&&h.has_on(pf2));

    std::list<Facet_handle> G = T2.facets_visible_from(p5);
    CGAL_TEST(G.size()==1);
    CGAL_TEST(T2.bounded_side(p5)==CGAL::ON_UNBOUNDED_SIDE && 
              T2.bounded_side(Point(1,1,1,10))==CGAL::ON_BOUNDED_SIDE);
    T2.is_valid();
    T2.clear(3);


  }
  CGAL_TEST_END;
}

