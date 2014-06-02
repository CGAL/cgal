#ifndef NAIVECOLLISIONDETECTOR_HEADER
#define NAIVECOLLISIONDETECTOR_HEADER
#include "ICollisionDetector.h"
#include <CGAL/intersections.h>
#include <CGAL/Boolean_set_operations_2.h>

template <class Kernel_, class Container_> class NaiveCollisionDetector : public ICollisionDetector< Kernel_,  Container_> {
public:
    NaiveCollisionDetector() {}
    typedef CGAL::Polygon_2<Kernel_> Polygon_2;
    typedef typename Polygon_2::Edge_const_iterator           Edge_iterator ;
    typedef typename Polygon_2::Traits::Segment_2           Segment_2 ;
    //typedef typename Polygon_2::Vertex_circulator          Vertex_circulator;
    //typedef typename
    virtual bool checkCollision(const Polygon_2 &p, const Polygon_2 &q) {
        bool intersect = false;

        if (CGAL::do_intersect(p, q)) {
            return true;
        } else {
            return false;
        }

        /*
        for (Edge_iterator itr1 = p.edges_begin();itr1!=p.edges_end();++itr1)
        {
            for (Edge_iterator itr2 = q.edges_begin();itr2!=q.edges_end();++itr2)
            {
                //intersect = intersect CGAL::Do_intersect(*itr1,*itr2);
                 CGAL::Object result = CGAL::intersection(*itr1,*itr2);
                 if (const CGAL::Point_2<Kernel> *ipoint = CGAL::object_cast<CGAL::Point_2<Kernel> >(&result)) {
                    // handle the point intersection case with *ipoint.
                    // check if the intersection point is the vertex of one of the edges
                     if (((*ipoint)==itr1->source()) || ((*ipoint)==itr1->target())
                         ||((*ipoint)==itr2->source()) || ((*ipoint)==itr2->target()))
                        intersect = false;
                     else
                        intersect = true;
                 } else
                if (const CGAL::Segment_2<Kernel> *iseg = CGAL::object_cast<CGAL::Segment_2<Kernel> >(&result)) {
                    // handle the segment intersection case with *iseg. this is the case where the edges are touching
                    intersect = false;
                } else {
                // handle the no intersection case.
                }

                if (intersect)
                    return intersect;
            }

        }

        if (p.bounded_side(*q.vertices_begin())== CGAL::ON_BOUNDED_SIDE || q.bounded_side(*p.vertices_begin())== CGAL::ON_BOUNDED_SIDE)
            return true;
        return intersect;
        */

    }
};


#endif
