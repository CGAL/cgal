#ifndef CGAL_HALFSPACES_INTERSECTION_WITH_CONSTRUCTION_3_H
#define CGAL_HALFSPACES_INTERSECTION_WITH_CONSTRUCTION_3_H

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
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
            const Polyhedron &_primal;

            public:
            Build_dual_polyhedron (const Polyhedron & primal):
                _primal (primal)
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
                    B.add_vertex(CGAL::ORIGIN + p.orthogonal_vector () / (-p.d()));
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

    template <class PlaneIterator, class Polyhedron>
    void
    halfspaces_intersection_with_construction_3(PlaneIterator pbegin,
                                                PlaneIterator pend,
                                                Polyhedron &P)
        {
            typedef typename Polyhedron::Traits::Kernel K;
            typedef typename K::Point_3 Point;
            typedef typename CGAL::internal::Build_dual_polyhedron<Polyhedron> Builder;

            // construct dual points to apply the convex hull
            std::vector<Point> dual_points;
            for (PlaneIterator p = pbegin; p != pend; ++p)
                dual_points.push_back(CGAL::ORIGIN + p->orthogonal_vector () / (-p->d()));

            Polyhedron ch;
            CGAL::convex_hull_3(dual_points.begin(), dual_points.end(), ch);

            Builder build_dual (ch);
            P.delegate(build_dual);
        }
} // namespace CGAL

#endif // CGAL_HALFSPACES_INTERSECTION_WITH_CONSTRUCTION_3_H

