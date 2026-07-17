// Copyright (c) 2025
// INRIA Nancy (France), and Université de Lorraine (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Vincent Despré, Camille Lanuel, Marc Pouget, Monique Teillaud

#ifndef CGAL_HYPERBOLIC_DIRICHLET_DOMAIN_2
#define CGAL_HYPERBOLIC_DIRICHLET_DOMAIN_2

#include <CGAL/license/Triangulation_on_hyperbolic_surface_2.h>

#include <CGAL/Delaunay_triangulation_on_hyperbolic_surface_2.h>

namespace CGAL {

// Input: triangulation with a single vertex v
// Output: lift of this triangulation with one lift of v mapped to the origin of the Poincaré disk and all its incident triangles around it.
template<class Traits>
std::vector<std::tuple<typename Delaunay_triangulation_on_hyperbolic_surface_2<Traits>::Dart_const_descriptor, typename Traits::Point_2, typename Traits::Point_2, typename Traits::Point_2>>
unfold(Delaunay_triangulation_on_hyperbolic_surface_2<Traits> & triangulation)
{
    typedef typename Traits::Complex                                   Complex;
    typedef typename Traits::Point_2                                   Point;
    typedef Delaunay_triangulation_on_hyperbolic_surface_2<Traits>     Delaunay_Triangulation;
    typedef typename Delaunay_Triangulation::Anchor                    Anchor;
    typedef typename Delaunay_Triangulation::CMap                      CMap;
    typedef typename Delaunay_Triangulation::Dart_const_descriptor     Dart_const_descriptor;
    typedef CGAL::Hyperbolic_isometry_2<Traits>                        Isometry;

    Anchor & anchor = triangulation.anchor();
    CMap & cmap = triangulation.combinatorial_map();

    std::vector<std::tuple<Dart_const_descriptor,Point,Point,Point>> lifted_triangles;  // vector of lifted triangles (future output)
    std::map<Dart_const_descriptor, Point> positions;  // map that will contain the computed lift of the vertex of each dart

    // create a mark for visited darts
    size_t visited = cmap.get_new_mark();
    cmap.unmark_all(visited);

    // translate the first vertex of the anchor to the origin of the Poincaré disk, and add the positions of the translated vertices of the anchor
    Isometry center_the_drawing = hyperbolic_translation<Traits>(anchor.vertices[0]);
    positions[anchor.dart] = center_the_drawing.evaluate(anchor.vertices[0]);
    positions[triangulation.const_ccw(anchor.dart)] = center_the_drawing.evaluate(anchor.vertices[1]);
    positions[triangulation.const_cw(anchor.dart)] = center_the_drawing.evaluate(anchor.vertices[2]);
    cmap.mark(anchor.dart, visited);

    // add the first triangle (the translated anchor) to the vector of triangles
    std::tuple<Dart_const_descriptor,Point,Point,Point> value = std::make_tuple(anchor.dart, positions[anchor.dart], positions[triangulation.const_ccw(anchor.dart)], positions[triangulation.const_cw(anchor.dart)]);
    lifted_triangles.push_back(value);

    // visit all the darts one by one by turning around the central vertex
    Dart_const_descriptor invader = anchor.dart;
    while( cmap.number_of_unmarked_darts(visited) > 1 ){  // >1 because the first triangle appears twice
        Dart_const_descriptor invaded = triangulation.const_opposite(invader);

        // get the positions of the vertices of the invader's triangle
        const Point& a = positions[triangulation.const_ccw(invader)];
        const Point& b = positions[triangulation.const_cw(invader)];
        const Point& c = positions[invader];
        Complex cross_ratio = triangulation.get_cross_ratio(invader);

        // retieve the positions of the invaded's triangle
        positions[invaded] = a;
        positions[triangulation.const_ccw(invaded)] = c;
        Point d = triangulation.fourth_point_from_cross_ratio(a, b, c, cross_ratio);
        positions[triangulation.const_cw(invaded)] = d;

        // add the three vertices to the vector of lifted triangles
        value = std::make_tuple(invaded, a, c, d);
        lifted_triangles.push_back(value);
        cmap.mark(invaded, visited);

        invader = triangulation.const_ccw(invaded);
    }

    cmap.free_mark(visited);
    return lifted_triangles;
}

// Input: Fundamental domain whose vertices are the same point on the surface
// Output: vertices of a Dirichlet domain centered at the origin of the Poincaré disk
template<class Traits>
std::vector<typename Traits::Hyperbolic_Voronoi_point_2> Dirichlet_vertices(Hyperbolic_fundamental_domain_2<Traits> & domain)
{
    typedef typename Traits::Point_2                                    Point;
    typedef typename Traits::Hyperbolic_Voronoi_point_2                 Voronoi_point;
    typedef Delaunay_triangulation_on_hyperbolic_surface_2<Traits>      Delaunay_Triangulation;
    typedef typename Delaunay_Triangulation::Dart_const_descriptor      Dart_const_descriptor;

    Traits gt;
    Delaunay_Triangulation triangulation = Delaunay_Triangulation(gt,domain);

    std::vector<std::tuple<Dart_const_descriptor,Point,Point,Point>> realized_triangles = unfold<Traits>(triangulation);
    std::vector<Voronoi_point> dirichlet_vertices;
    dirichlet_vertices.reserve(realized_triangles.size());
    typename Traits::Construct_hyperbolic_circumcenter_2 chc = gt.construct_hyperbolic_circumcenter_2_object();
    for (std::tuple<Dart_const_descriptor, Point, Point, Point>& triangle : realized_triangles){
        Voronoi_point circumcenter = chc(std::get<1>(triangle), std::get<2>(triangle), std::get<3>(triangle));
        dirichlet_vertices.push_back(circumcenter);
    }
    return dirichlet_vertices;
}

//////////////////////////////////////////////////////
//       TO JSON OUTPUT
//////////////////////////////////////////////////////
template<class Traits>
std::ostream&
  Dirichlet_to_json(Hyperbolic_fundamental_domain_2<Traits> & domain, std::ostream& s)
{
    typedef typename Traits::Point_2                                    Point;
    typedef typename Traits::Hyperbolic_Voronoi_point_2                 Voronoi_point;

    std::vector<Voronoi_point> vertices = Dirichlet_vertices(domain);

    const std::size_t n = vertices.size();

    s << "{\n";
    s << "  \"type\": " << "Dirichlet" << ",\n";
    s << "  \"size\": " << n << ",\n";

    s << "  \"vertices\": [";
    for (std::size_t k = 0; k < n; ++k)
    {
        if (k > 0) s << "," << std::endl;
        s << "[" << to_double(vertices[k].x()) << ", " << to_double(vertices[k].y()) << "]" ;
    }
    s << "]\n";

    s << "}";

    return s;
}

} // namespace CGAL

#endif  //CGAL_HYPERBOLIC_DIRICHLET_DOMAIN_2
