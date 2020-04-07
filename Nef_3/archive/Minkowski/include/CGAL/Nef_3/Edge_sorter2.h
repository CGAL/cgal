#ifndef CGAL_NEF3_EDGE_SORTER2_H
#define CGAL_NEF3_EDGE_SORTER2_H

#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_3/SNC_point_locator.h>


namespace CGAL {

template<typename Nef_, typename Container>
class Edge_sorter2 : public Modifier_base<typename Nef_::SNC_and_PL> {

  typedef Nef_                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL    SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure SNC_structure;
  typedef typename Nef_polyhedron::Items         Items;
  typedef CGAL::SNC_decorator<SNC_structure>     Base;
  typedef CGAL::SNC_point_locator<Base>          SNC_point_locator;
  typedef CGAL::SNC_constructor<Items, SNC_structure>
    SNC_constructor;

  typedef typename SNC_structure::Vertex_handle     Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle   SVertex_handle;
  typedef typename SNC_structure::Halfedge_iterator SVertex_iterator;
  typedef typename SNC_structure::Halfedge_handle   Halfedge_handle;
  typedef typename SNC_structure::SHalfedge_handle   SHalfedge_handle;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Vector_3 Vector_3;

  typedef typename Container::iterator Iterator;

  struct sort_edges {
    bool operator()(Halfedge_handle e1, Halfedge_handle e2) {
      return e1->source()->point() > e2->source()->point();
    }
  };

  bool split_at(Segment_3 s1, Segment_3 s2, Point_3& ip2) {

    Point_3 ip1;
    Vector_3 vec1(cross_product(s2.to_vector(),Vector_3(1,0,0)));
    Plane_3 pl1(s2.source(), vec1);
    Object o1 = intersection(pl1,s1);
    if(!assign(ip1,o1))
      return false;
    Vector_3 vec2(cross_product(s1.to_vector(),Vector_3(1,0,0)));
    Plane_3 pl2(s1.source(), vec2);
    Object o2 = intersection(pl2,s2);
    if(!assign(ip2,o2) || ip2 == s2.source() || ip2 == s2.target())
      return false;
    if(ip2 > ip1)
      return true;
    return false;
  }

  Container& c;
  SNC_structure* sncp;
  SNC_point_locator* pl;

 public:
  Edge_sorter2(Container& cin) : c(cin) {}

  void operator()(SNC_and_PL& sncpl) {

    sncp = sncpl.sncp;
    pl = sncpl.pl;
    SNC_constructor C(*sncp);

    //    std::cerr << "edge_sorter2 " << c.size() << std::endl;
    std::sort(c.begin(), c.end(), sort_edges());
    Iterator esi1,esi2,esi3;
    for(esi1 = c.begin(); esi1 != c.end(); ++esi1) {
      //      std::cerr << "1: " << (*esi1)->source()->point() << "->" << (*esi1)->twin()->source()->point() << std::endl;
      esi2 = esi1;
      ++esi2;
      //      if(esi2 != c.end())
      //        std::cerr << "2: " << (*esi2)->source()->point() << "->" << (*esi2)->twin()->source()->point() << std::endl;
      while(esi2!=c.end() &&
            (*esi1)->twin()->source()->point().x() <
            (*esi2)->source()->point().x()) {
        if((*esi1)->source() == (*esi2)->source()) {
          ++esi2;
          continue;
        }
        Point_3 ip;
        bool b = split_at(Segment_3((*esi1)->source()->point(),
                                    (*esi1)->twin()->source()->point()),
                          Segment_3((*esi2)->source()->point(),
                                    (*esi2)->twin()->source()->point()),ip);
        //        std::cerr << "split " << b << std::endl;
        if(b) {
          //          std::cerr << ip << std::endl;
          Vertex_handle v;
          Halfedge_handle e = (*esi2);
          v = C.create_from_edge(e,ip);
          pl->add_vertex(v);

          SVertex_iterator svi = v->svertices_begin();
          SVertex_handle svf, svb;
          if(svi->point() == e->point()) {
            svf = svi;
            svb = ++svi;
          } else {
            svb = svi;
            svf = ++svi;
          }

          svb->twin() = e;
          svf->twin() = e->twin();
          e->twin()->twin() = svf;
          e->twin() = svb;

          pl->add_edge(svf);
          pl->add_edge(svb);

          esi3 = esi2;
          ++esi3;
          while(esi3 != c.end() &&
                (*esi3)->source()->point() >
                svf->source()->point())
            ++esi3;
          c.insert(esi3, svf);
        }

        ++esi2;
        //        std::cerr << "2: " << (*esi2)->source()->point() << "->" << (*esi2)->twin()->source()->point() << std::endl;
      }
    }
  }
};

} //namespace CGAL
#endif //CGAL_NEF3_EDGE_SORTER2_H
