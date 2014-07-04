#ifndef CGAL_HALFSPACES_INTERSECTION_H
#define CGAL_HALFSPACES_INTERSECTION_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/dual/Convex_hull_traits_dual_3.h>
#include <CGAL/Origin.h>
#include <CGAL/Convex_hull_3.h>
#include <CGAL/intersections.h>
#include <CGAL/assertions.h>
#include <CGAL/Point_inside_polyhedron_3.h>

namespace CGAL
{
    namespace Convex_hull_3
    {
        namespace internal
        {
            // Build the primal polyhedron associated to a dual polyhedron
            // We also need the `origin` that is to say a point inside the primal polyhedron
            template <typename K, class Polyhedron_dual, class Polyhedron>
                class Build_primal_polyhedron :
                    public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> {
                        typedef typename Polyhedron::HalfedgeDS HDS;
                        const Polyhedron_dual & _dual;

                        // Origin
                        typedef typename K::Point_3 Primal_point_3;
                        Primal_point_3 origin;

                        public:
                        Build_primal_polyhedron (const Polyhedron_dual & dual,
                                                 Primal_point_3 o =
                                                 Primal_point_3(0, 0, 0)) : _dual (dual), origin(o)
                        {}

                        void operator () (HDS &hds)
                        {
                            typedef typename K::RT RT;
                            typedef typename K::Point_3 Point_3;

                            // Typedefs for dual
                            typedef typename Polyhedron_dual::Facet Facet;
                            typedef typename Polyhedron_dual::Facet_const_handle
                                Facet_const_handle;
                            typedef typename Polyhedron_dual::Facet_const_iterator
                                Facet_const_iterator;
                            typedef typename Polyhedron_dual::Vertex_const_iterator
                                Vertex_const_iterator;

                            // Typedefs for primal
                            typename CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

                            // Typedefs for intersection
                            typedef typename K::Plane_3 Plane_3;
                            typedef typename K::Line_3 Line_3;
                            typedef boost::optional< boost::variant< Point_3,
                                                                     Line_3,
                                                                     Plane_3 > > result_inter;

                            B.begin_surface(_dual.size_of_facets(),
                                            _dual.size_of_vertices(),
                                            _dual.size_of_vertices());

                            std::map <Facet_const_handle, size_t> primal_vertices;
                            size_t n = 0;

                            // First, computing the primal vertices
                            for (Facet_const_iterator it = _dual.facets_begin();
                                 it != _dual.facets_end(); ++it, ++n) {
                                typename Facet::Halfedge_const_handle h = it->halfedge();
                                // Build the dual plane corresponding to the current facet
                                Plane_3 p1 = h->vertex()->point();
                                Plane_3 p2 = h->next()->vertex()->point();
                                Plane_3 p3 = h->next()->next()->vertex()->point();

                                RT dp1 = p1.d() + origin.x() * p1.a()
                                    + origin.y() * p1.b() + origin.z() * p1.c();
                                RT dp2 = p2.d() + origin.x() * p2.a()
                                    + origin.y() * p2.b() + origin.z() * p2.c();
                                RT dp3 = p3.d() + origin.x() * p3.a()
                                    + origin.y() * p3.b() + origin.z() * p3.c();

                                Plane_3 pp1(p1.a(), p1.b(), p1.c(), dp1);
                                Plane_3 pp2(p2.a(), p2.b(), p2.c(), dp2);
                                Plane_3 pp3(p3.a(), p3.b(), p3.c(), dp3);

                                // Compute the intersection
                                result_inter result = CGAL::intersection(pp1, pp2, pp3);
                                CGAL_assertion_msg(result,
                                                   "halfspaces_intersection: no intersection");
                                CGAL_assertion_msg(boost::get<Point_3>(& *result),
                                                   "halfspaces_intersection: intersection is not a point");

                                const Point_3* pp = boost::get<Point_3>(& *result);

                                // Primal vertex associated to the current dual plane
                                Point_3 ppp(origin.x() + pp->x(),
                                            origin.y() + pp->y(),
                                            origin.z() + pp->z());

                                B.add_vertex(ppp);
                                primal_vertices[it] = n;
                            }

                            // Then, add facets to the primal polyhedron
                            // To do this, for each dual vertex, we circulate around this vertex
                            // and we add an edge between each facet we encounter
                            for (Vertex_const_iterator it = _dual.vertices_begin();
                                 it != _dual.vertices_end(); ++it, ++n) {
                                typename Polyhedron_dual::Halfedge_around_vertex_const_circulator
                                    h0 = it->vertex_begin(), hf = h0;
                                B.begin_facet();
                                do {
                                    B.add_vertex_to_facet(primal_vertices[hf->facet()]);
                                } while (++hf != h0);
                                B.end_facet();
                            }

                            B.end_surface();
                        }
                    };
        } // namespace internal
    } // namespace Convex_hull_3

    // Compute the intersection of halfspaces
    template <class PlaneIterator, class Polyhedron>
        void halfspaces_intersection (PlaneIterator begin, PlaneIterator end,
                                      Polyhedron &P,
                                      typename Polyhedron::Traits::Point_3 const& origin = typename Polyhedron::Traits::Point_3(CGAL::ORIGIN)) {
            typedef typename Polyhedron::Traits::Kernel K;
            typedef Convex_hull_3::Convex_hull_traits_dual_3<K> Hull_traits_dual_3;
            typedef Polyhedron_3<Hull_traits_dual_3> Polyhedron_dual_3;
            typedef Convex_hull_3::internal::Build_primal_polyhedron<K, Polyhedron_dual_3, Polyhedron> Builder;

            Hull_traits_dual_3 dual_traits(origin);

            Polyhedron_dual_3 dual_convex_hull;
            CGAL::convex_hull_3(begin, end, dual_convex_hull, dual_traits);
            Builder build_primal(dual_convex_hull, origin);
            P.delegate(build_primal);

            // Posterior check for the origin inside the cmputed polyhedron
            Point_inside_polyhedron_3<Polyhedron, K> is_inside(P);
            CGAL_assertion_msg(is_inside(origin) == CGAL::ON_BOUNDED_SIDE,
                               "halfspaces_intersection: origin not in the polyhedron");
        }
} // namespace CGAL

#endif // CGAL_HALFSPACES_INTERSECTION_H

