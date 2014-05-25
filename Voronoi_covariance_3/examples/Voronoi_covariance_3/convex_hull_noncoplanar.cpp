#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Voronoi_covariance_3/Convex_hull_traits_dual_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Voronoi_covariance_3::Convex_hull_traits_dual_3<K>    Hull_traits_dual_3;

typedef CGAL::Polyhedron_3<Hull_traits_dual_3>                      Polyhedron_dual_3;

#include <fstream>
#include <iostream>
#include "include/to_dual.h"

// Write a polyhedron into an OFF file
template <typename R>
void convertToOFF (std::string const& filename, Polyhedron_dual_3& P) {
    typedef Polyhedron_dual_3::Point_iterator                    Point_iterator;
    typedef Polyhedron_dual_3::Facet_iterator                    Facet_iterator;
    typedef Polyhedron_dual_3::Halfedge_handle                   Halfedge_handle;
    typedef Polyhedron_dual_3::Halfedge_around_facet_circulator  Halfedge_facet_circulator;

    std::ofstream file(filename.c_str());

    // 0. number of vertices / number of facets / number of edges
    file << "OFF" << std::endl << P.size_of_vertices() << ' '
        << P.size_of_facets() << " 0" << std::endl;

    // 1. vertices definition
    for (Point_iterator vit = P.points_begin(); vit != P.points_end(); vit++) {
        typename R::Plane_3 p = *vit;
        file << to_dual<R>(p) << std::endl;
    }

    // 2. facets definition
    for (  Facet_iterator i = P.facets_begin(); i != P.facets_end(); ++i) {
        Halfedge_facet_circulator j = i->facet_begin();
        // Facets in polyhedral surfaces are at least triangles.
        CGAL_assertion( CGAL::circulator_size(j) >= 3);
        file << CGAL::circulator_size(j) << ' ';
        do {
            file << ' ' << std::distance(P.vertices_begin(), j->vertex());
        } while ( ++j != i->facet_begin());
        file << std::endl;
    }

    file.close();
}

int main (void) {
    // define polyhedron to hold convex hull
    Polyhedron_dual_3 dual_poly;

    // traits
    Hull_traits_dual_3 dual_traits;

    // dual convex hull
    std::list<K::Plane_3> planes;
    planes.push_back(K::Plane_3(1, 0, 0, -1));
    planes.push_back(K::Plane_3(1, 0, 0, 1));
    planes.push_back(K::Plane_3(0, 1, 0, -1));
    planes.push_back(K::Plane_3(0, 1, 0, 1));
    planes.push_back(K::Plane_3(0, 0, 1, -1));
    planes.push_back(K::Plane_3(0, 0, 1, 1));

    // Data structures
    typedef CGAL::Triangulation_data_structure_2<
        CGAL::Triangulation_vertex_base_with_info_2<int, CGAL::GT3_for_CH3<Hull_traits_dual_3> >,
        CGAL::Convex_hull_face_base_2<int, Hull_traits_dual_3> >                           Tds;
    Tds dual_tds;
    typedef std::list<K::Plane_3>::iterator plane_iter;
    typedef Tds::Vertex_handle Vertex_handle;
    typedef Tds::Face_handle Face_handle;

    // Inspired from code of Convex_hull_3.h:879-916
    // Find three non collinear dual points
    Hull_traits_dual_3::Equal_3 equal = dual_traits.equal_3_object();
    Hull_traits_dual_3::Collinear_3 collinear = dual_traits.collinear_3_object();

    plane_iter plane1_it = planes.begin();
    plane_iter plane2_it = planes.begin();
    plane2_it++;
    while (plane2_it != planes.end() && equal(*plane1_it,*plane2_it))
        plane2_it++;
    plane_iter plane3_it = plane2_it;
    plane3_it++;
    while (plane3_it != planes.end() && collinear(*plane1_it,*plane2_it,*plane3_it))
        ++plane3_it;

    std::cout << "3 Non collinear dual points :" << std::endl;
    std::cout << *plane1_it << std::endl;
    std::cout << *plane2_it << std::endl;
    std::cout << *plane3_it << std::endl;

    // Find a fourth dual point non coplanar to the 3 others
    // ch_quickhull_polyhedron_3:670
    /* CGAL::internal::Convex_hull_3::ch_quickhull_polyhedron_3(planes, plane1_it, plane2_it, plane3_it, dual_poly, dual_traits); */
    typedef Hull_traits_dual_3::Plane_3 Plane_3;
    Hull_traits_dual_3::Construct_plane_3 construct_plane = dual_traits.construct_plane_3_object();
    Plane_3 plane = construct_plane(*plane1_it, *plane2_it, *plane3_it);
    typedef Hull_traits_dual_3::Less_signed_distance_to_plane_3 Dist_compare;
    Dist_compare compare_dist = dual_traits.less_signed_distance_to_plane_3_object();

    Hull_traits_dual_3::Coplanar_3 coplanar = dual_traits.coplanar_3_object();
    std::pair<plane_iter, plane_iter> min_max;
    min_max = CGAL::min_max_element(planes.begin(), planes.end(),
                                    boost::bind(compare_dist, plane, _1, _2),
                                    boost::bind(compare_dist, plane, _1, _2));
    plane_iter plane4_it;
    if (coplanar(*plane1_it, *plane2_it, *plane3_it, *min_max.second)) {
        std::cout << "COPLANAR -> SWAP" << std::endl;
        plane4_it = min_max.first;
        std::swap(*plane1_it, *plane3_it);
    } else {
        std::cout << "COPLANAR -> NON SWAP" << std::endl;
        plane4_it = min_max.second;
        // TODO: remove next line
        std::swap(*plane3_it, *plane2_it);
    }
    std::cout << "Fourth point non coplanar to the others :" << std::endl;
    std::cout << *plane4_it << std::endl;

    if (coplanar(*plane1_it, *plane2_it, *plane3_it, *plane4_it)) {
        std::cout << "ERROR: coplanar" << std::endl;
        return 1;
    }

    std::cout << "Planes in order :" << std::endl;
    std::cout << *plane1_it << std::endl;
    std::cout << *plane2_it << std::endl;
    std::cout << *plane3_it << std::endl;
    std::cout << *plane4_it << std::endl;

    // Create the associated TDS
    Vertex_handle v0 = dual_tds.create_vertex(); v0->set_point(*plane1_it);
    Vertex_handle v1 = dual_tds.create_vertex(); v1->set_point(*plane2_it);
    Vertex_handle v2 = dual_tds.create_vertex(); v2->set_point(*plane3_it);
    Vertex_handle v3 = dual_tds.create_vertex(); v3->set_point(*plane4_it);
    v0->info() = v1->info() = v2->info() = v3->info() = 0;
    Face_handle f0 = dual_tds.create_face(v0,v1,v2);
    Face_handle f1 = dual_tds.create_face(v3,v1,v0);
    Face_handle f2 = dual_tds.create_face(v3,v2,v1);
    Face_handle f3 = dual_tds.create_face(v3,v0,v2);
    dual_tds.set_dimension(2);
    f0->set_neighbors(f2, f3, f1);
    f1->set_neighbors(f0, f3, f2);
    f2->set_neighbors(f0, f1, f3);
    f3->set_neighbors(f0, f2, f1);

    // Remove unuseful planes
    planes.erase(plane1_it);
    planes.erase(plane2_it);
    planes.erase(plane3_it);
    planes.erase(plane4_it);

    // Print planes
    std::cout << "Planes :" << std::endl;
    for (plane_iter pit = planes.begin(); pit != planes.end(); pit++) {
        std::cout << *pit << std::endl;
    }

    // Compute the convex hull
    CGAL::internal::Convex_hull_3::non_coplanar_quickhull_3(planes, dual_tds, dual_traits);

    // Build the associated polyhedron
    CGAL::internal::Convex_hull_3::internal::Build_convex_hull_from_TDS_2<Polyhedron_dual_3::HalfedgeDS,Tds> builder(dual_tds);
    dual_poly.delegate(builder);

    // Print the polyhedron in an OFF file
    convertToOFF<K>("test.off", dual_poly);

    return 0;
}

