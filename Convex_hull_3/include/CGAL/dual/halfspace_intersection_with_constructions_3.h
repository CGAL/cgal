#ifndef CGAL_HALFSPACES_INTERSECTION_WITH_CONSTRUCTION_3_H
#define CGAL_HALFSPACES_INTERSECTION_WITH_CONSTRUCTION_3_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Origin.h>
#include <CGAL/convex_hull_3.h>

namespace CGAL
{
    namespace internal
    {
        template <class Polyhedron>
        class Build_dual_polyhedron :
                public CGAL::Modifier_base<typename Polyhedron::HalfedgeDS>
        {
            typedef typename Polyhedron::HalfedgeDS HDS;
            typedef typename Polyhedron::Traits::Point_3 Point_3;
            const Polyhedron &_primal;
            Point_3 origin;

            public:
            Build_dual_polyhedron (const Polyhedron & primal,
                                   Point_3 o = Point_3(CGAL::ORIGIN)):
                _primal (primal), origin(o)
            {}

            void operator () (HDS &hds)
            {
                typedef typename Polyhedron::Facet Facet;
                typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
                typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;
                typedef typename Polyhedron::Vertex_const_iterator Vertex_const_iterator;
                typename CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);

                B.begin_surface(_primal.size_of_facets(),
                                _primal.size_of_vertices(),
                                _primal.size_of_vertices());

                // compute coordinates of extreme vertices in the dual polyhedron
                // from primal faces
                std::map<Facet_const_handle, size_t> extreme_points;
                size_t n = 0;

                for (Facet_const_iterator it = _primal.facets_begin();
                     it != _primal.facets_end(); ++it, ++n)
                {
                    typename Facet::Halfedge_const_handle h = it->halfedge();
                    typename Facet::Plane_3 p ( h->vertex()->point(),
                                                h->next()->vertex()->point(),
                                                h->next()->next()->vertex()->point());
                    // translate extreme vertex
                    Point_3 extreme_p = CGAL::ORIGIN + p.orthogonal_vector () / (-p.d());
                    Point_3 translated_extreme_p(extreme_p.x() + origin.x(),
                                                 extreme_p.y() + origin.y(),
                                                 extreme_p.z() + origin.z());
                    B.add_vertex(translated_extreme_p);
                    extreme_points[it] = n;
                }

                // build faces
                for (Vertex_const_iterator it = _primal.vertices_begin();
                     it != _primal.vertices_end(); ++it)
                {
                    assert (it->is_bivalent() == false);

                    typename Polyhedron::Halfedge_around_vertex_const_circulator
                        h0 = it->vertex_begin(), hf = h0;

                    B.begin_facet();
                    do
                    {
                        B.add_vertex_to_facet(extreme_points[hf->facet()]);
                    } while (++hf != h0);
                    B.end_facet();
                }

                B.end_surface();
            }
        };
    } // namespace internal

    template <class PlaneIterator, class Polyhedron, class Traits>
    void halfspace_intersection_with_constructions_3(PlaneIterator pbegin,
                                                     PlaneIterator pend,
                                                     Polyhedron &P,
                                                     const Traits & ch_traits,
                                                     typename Polyhedron::Vertex::Point_3 const& origin)
            {
            typedef typename Kernel_traits<typename Polyhedron::Vertex::Point_3>::Kernel K;
            typedef typename K::Point_3 Point;
            typedef typename K::Plane_3 Plane;
            typedef typename CGAL::internal::Build_dual_polyhedron<Polyhedron> Builder;

            // construct dual points to apply the convex hull
            std::vector<Point> dual_points;
            for (PlaneIterator p = pbegin; p != pend; ++p) {
                // translate plane
                Plane translated_p(p->a(),
                                   p->b(),
                                   p->c(),
                                   p->d() + origin.x() * p->a() + origin.y() * p->b() + origin.z() * p->c());
                dual_points.push_back(CGAL::ORIGIN + translated_p.orthogonal_vector () / (-translated_p.d()));
            }

            Polyhedron ch;
            CGAL::convex_hull_3(dual_points.begin(), dual_points.end(), ch, ch_traits);

            Builder build_dual (ch, origin);
            P.delegate(build_dual);
        }

    template <class PlaneIterator, class Polyhedron>
    void halfspace_intersection_with_constructions_3(PlaneIterator pbegin,
                                                     PlaneIterator pend,
                                                     Polyhedron &P,
                                                     typename Polyhedron::Vertex::Point_3 const& origin)
    {
        typedef typename Kernel_traits<typename Polyhedron::Vertex::Point_3>::Kernel K;
        typedef typename K::Point_3 Point_3;
        typedef typename internal::Convex_hull_3::Default_traits_for_Chull_3<Point_3>::type Traits;

        halfspace_intersection_with_constructions_3(pbegin, pend, P, Traits(), origin);
    }
} // namespace CGAL

#endif // CGAL_HALFSPACES_INTERSECTION_WITH_CONSTRUCTION_3_H

