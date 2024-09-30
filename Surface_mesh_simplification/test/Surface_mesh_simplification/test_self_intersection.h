#include "basics.h"

#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/intersections.h>

typedef CGAL::Bbox_3                                                                Bbox;
typedef Kernel::Triangle_3                                                          Triangle;
typedef Kernel::Segment_3                                                           Segment;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double, 3, Facet_const_handle>  Box;

std::vector<Triangle> triangles;

struct Intersect_facets
{
  void operator()(const Box* b, const Box* c) const
  {
    Halfedge_const_handle h = b->handle()->halfedge();
    // check for shared edge --> no intersection
    if(h->opposite()->facet() == c->handle()
          || h->next()->opposite()->facet() == c->handle()
          || h->next()->next()->opposite()->facet() == c->handle())
        return;

    // check for shared vertex --> maybe intersection, maybe not
    Halfedge_const_handle g = c->handle()->halfedge();
    Halfedge_const_handle v;
    if(h->vertex() == g->vertex())
        v = g;
    if(h->vertex() == g->next()->vertex())
        v = g->next();
    if(h->vertex() == g->next()->next()->vertex())
        v = g->next()->next();
    if(v == Halfedge_const_handle()) {
        h = h->next();
        if(h->vertex() == g->vertex())
            v = g;
        if(h->vertex() == g->next()->vertex())
            v = g->next();
        if(h->vertex() == g->next()->next()->vertex())
            v = g->next()->next();
        if(v == Halfedge_const_handle()) {
            h = h->next();
            if(h->vertex() == g->vertex())
                v = g;
            if(h->vertex() == g->next()->vertex())
                v = g->next();
            if(h->vertex() == g->next()->next()->vertex())
                v = g->next()->next();
        }
    }

    if(v != Halfedge_const_handle()) {
        // found shared vertex:
        assert(h->vertex() == v->vertex());
        // geometric check if the opposite segments intersect the triangles
        Triangle t1(h->vertex()->point(),
                      h->next()->vertex()->point(),
                      h->next()->next()->vertex()->point());
        Triangle t2(v->vertex()->point(),
                      v->next()->vertex()->point(),
                      v->next()->next()->vertex()->point());
        Segment  s1(h->next()->vertex()->point(),
                      h->next()->next()->vertex()->point());
        Segment  s2(v->next()->vertex()->point(),
                      v->next()->next()->vertex()->point());
        if(CGAL::do_intersect(t1, s2)) {
            //cerr << "Triangles intersect (t1,s2):\n    T1: " << t1
            //     << "\n    T2 :" << t2 << endl;
            triangles.push_back(t1);
            triangles.push_back(t2);
        } else if(CGAL::do_intersect(t2, s1)) {
            //cerr << "Triangles intersect (t2,s1):\n    T1: " << t1
            //     << "\n    T2 :" << t2 << endl;
            triangles.push_back(t1);
            triangles.push_back(t2);
        }
        return;
    }

    // check for geometric intersection
    Triangle t1(h->vertex()->point(),
                  h->next()->vertex()->point(),
                  h->next()->next()->vertex()->point());
    Triangle t2(g->vertex()->point(),
                  g->next()->vertex()->point(),
                  g->next()->next()->vertex()->point());

    if(CGAL::do_intersect(t1, t2)) {
        //cerr << "Triangles intersect:\n    T1: " << t1 << "\n    T2 :"
        //     << t2 << endl;
        triangles.push_back(t1);
        triangles.push_back(t2);
    }
  }
};

bool Is_self_intersecting(const Surface& s)
{
  std::vector<Box> boxes;
  boxes.reserve(s.size_of_facets());
  for(Facet_const_iterator i = s.facets_begin(); i != s.facets_end(); ++i)
  {
      boxes.push_back(
          Box(i->halfedge()->vertex()->point().bbox()
              + i->halfedge()->next()->vertex()->point().bbox()
              + i->halfedge()->next()->next()->vertex()->point().bbox(),
                i));
  }
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(s.size_of_facets());
  for(std::vector<Box>::iterator j = boxes.begin(); j != boxes.end(); ++j)
  {
      box_ptr.push_back(&*j);
  }
  CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(),
                                  Intersect_facets(), std::ptrdiff_t(2000));

  return triangles.size() > 0;
}
