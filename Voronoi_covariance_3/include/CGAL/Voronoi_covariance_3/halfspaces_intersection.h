#ifndef HALFSPACES_INTERSECTION_H
#define HALFSPACES_INTERSECTION_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Voronoi_covariance_3/Convex_hull_traits_dual_3.h>

namespace CGAL
{
    namespace Voronoi_covariance_3
    {
        namespace internal
        {
            template <typename K, class Polyhedron_dual, class Polyhedron>
                class Build_primal_polyhedron :
                    public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS> {
                        typedef typename Polyhedron::HalfedgeDS HDS;
                        const Polyhedron_dual & _dual;

                        public:
                        Build_primal_polyhedron (const Polyhedron_dual & dual) : _dual (dual)
                        {}

                        // Compute the primal point associated to a triple of dual planes

                        void operator () (HDS &hds)
                        {
                            typedef typename K::RT RT;
                            typedef typename K::Point_3 Point_3;

                            // Typedefs for dual
                            typedef typename Polyhedron_dual::Facet Facet;
                            typedef typename Polyhedron_dual::Vertex Vertex;
                            typedef typename Vertex::Point_3 Plane;
                            typedef typename Polyhedron_dual::Facet_const_handle
                                Facet_const_handle;
                            typedef typename Polyhedron_dual::Facet_const_iterator
                                Facet_const_iterator;
                            typedef typename Polyhedron_dual::Vertex_const_iterator
                                Vertex_const_iterator;

                            // Typedefs for primal
                            typename CGAL::Polyhedron_incremental_builder_3<HDS> B(hds, true);

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
                                Plane p1 = h->vertex()->point();
                                Plane p2 = h->next()->vertex()->point();
                                Plane p3 = h->next()->next()->vertex()->point();

                                // Normal to the dual plane
                                RT alpha = (p1.d() * p2.b() - p2.d() * p1.b()) *
                                    (p1.d() * p3.c() - p3.d() * p1.c()) -
                                    (p1.d() * p2.c() - p2.d() * p1.c()) *
                                    (p1.d() * p3.b() - p3.d() * p1.b());

                                RT beta  = (p1.d() * p2.c() - p2.d() * p1.c()) *
                                    (p1.d() * p3.a() - p3.d() * p1.a()) -
                                    (p1.d() * p2.a() - p2.d() * p1.a()) *
                                    (p1.d() * p3.c() - p3.d() * p1.c());

                                RT gamma = (p1.d() * p2.a() - p2.d() * p1.a()) *
                                    (p1.d() * p3.b() - p3.d() * p1.b()) -
                                    (p1.d() * p2.b() - p2.d() * p1.b()) *
                                    (p1.d() * p3.a() - p3.d() * p1.a());

                                // last coefficient of the dual plane equation
                                RT d = (alpha * p1.a() + beta * p1.b() + gamma * p1.c()) / p1.d();

                                // Primal vertex associated to the current dual plane
                                // TODO: add origin
                                // TODO: replace by CGAL::intersection
                                Point_3 p(-alpha / d, -beta / d, -gamma / d);

                                B.add_vertex(p);
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

        template <class PlaneIterator, class Polyhedron, class K>
            void halfspaces_intersection (PlaneIterator begin, PlaneIterator end,
                                          Polyhedron &P, const K &k, typename K::Point_3 const& origin = typename K::Point_3(0, 0, 0)) {
                typedef Convex_hull_traits_dual_3<K> Hull_traits_dual_3;
                typedef Polyhedron_3<Hull_traits_dual_3> Polyhedron_dual_3;
                typedef internal::Build_primal_polyhedron<K, Polyhedron_dual_3, Polyhedron> Builder;

                Hull_traits_dual_3 dual_traits;

                Polyhedron_dual_3 dual_convex_hull;
                /* convex_hull_3(begin, end, dual_convex_hull, dual_traits); */
                Builder build_primal(dual_convex_hull);
                P.delegate(build_primal);
            }

    } // namespace Voronoi_covariance_3
} // namespace CGAL

#endif

