// Copyright (c) 2005-2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_CD3_EDGE_SORTER_H
#define CGAL_CD3_EDGE_SORTER_H

#include <CGAL/license/Convex_decomposition_3.h>


#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_point_locator.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Convex_decomposition_3/Insert_vertex_into_edge.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 547
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename Halfedge_handle, typename Point_comparison>
class Compare_halfedges_in_reflex_edge_sorter
{
public:
        typedef Point_comparison Compare_points;
        Compare_points comp;
        Compare_halfedges_in_reflex_edge_sorter() {}
        bool operator()(Halfedge_handle e0, Halfedge_handle e1) const
        {
                return comp(e0->source()->point(), e1->source()->point());
        }
};

template<typename Point_3,
  typename FTComparison,
  typename SplitTest,
  typename VertexInserter,
  typename Container>
class Generic_edge_sorter
{
        typedef typename Container::iterator Iterator;
        typedef typename Container::value_type Edge_handle;
        typedef typename Container::key_compare EdgeComparison;
        typedef typename EdgeComparison::Compare_points PointComparison;

        PointComparison compare_points;
        FTComparison compare_fts;

public:
        Generic_edge_sorter() {}

        void operator()(Container& c, SplitTest& st, VertexInserter& vi)
        {
//                CGAL_NEF_SETDTHREAD(547);

                CGAL_assertion_code(
                        for(Iterator ti = c.begin(); ti != c.end(); ++ti)
                                CGAL_assertion(compare_points((*ti)->source()->point(),
                                                                                          (*ti)->twin()->source()->point())));
                CGAL_NEF_TRACEN("edge_sorter " << c.size());

//                std::sort(c.begin(), c.end(), Sort_edges());

                Iterator esi1,esi2,esi3;
                for(esi1 = c.begin(); esi1 != c.end(); ++esi1)
                {
                        CGAL_NEF_TRACEN("1: " << (*esi1)->source()->point() <<
                                                        "->" << (*esi1)->twin()->source()->point());
                        esi2 = esi1;
                        ++esi2;
                        CGAL_assertion_code(
                                if(esi2!=c.end())
                                        CGAL_NEF_TRACEN("2: " << (*esi2)->source()->point() <<
                                                                        "->" << (*esi2)->twin()->source()->point()););
                        while(esi2!=c.end() &&
                                  compare_fts((*esi2)->source()->point().x(),
                                                          (*esi1)->twin()->source()->point().x()))
                        {
                                if((*esi1)->source() == (*esi2)->source())
                                {
                                        ++esi2;
                                        continue;
                                }
                                Point_3 ip;
                                bool b = st(*esi1, *esi2, ip);
                                if(b)
                                {
                                        CGAL_NEF_TRACEN("ip " << ip);

                                        Edge_handle e = *esi2;
                                        Edge_handle svf = vi(e, ip);

                                        esi3 = esi2;
                                        ++esi3;
                                        while(esi3 != c.end() &&
                                                  compare_points((*esi3)->source()->point(),
                                                                                 svf->source()->point()))
                                                ++esi3;

                                        c.insert(esi3, svf);
                                }

                                ++esi2;

                                CGAL_assertion_code(
                                        if(esi2!=c.end())
                                                CGAL_NEF_TRACEN("y2: " << (*esi2)->source()->point() <<
                                                                                "->" << (*esi2)->twin()->source()->point()););
                        } // while
                } //while
        } // operator()

};

template<class Kernel, class PointComparison>
class Need_to_split {

  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Plane_3 Plane_3;
  typedef typename Kernel::Vector_3 Vector_3;

  PointComparison compare_points;
 public:
  Need_to_split() {}

  template<typename Edge_handle>
  bool operator()(Edge_handle e1, Edge_handle e2, Point_3& ip2) const {
    Segment_3 s1(e1->source()->point(), e1->twin()->source()->point());
    Segment_3 s2(e2->source()->point(), e2->twin()->source()->point());
    return operator()(s1, s2, ip2);
  }

  bool operator()(Segment_3 s1, Segment_3 s2, Point_3& ip2) const {
    CGAL_NEF_TRACEN("need_to_split " << s1);
    CGAL_NEF_TRACEN("         " << s2);
    Point_3 ip1;
    Vector_3 vec1(cross_product(s2.to_vector(),Vector_3(1,0,0)));
    Plane_3 pl1(s2.source(), vec1);
    CGAL_NEF_TRACEN("pl1 " << pl1);
    CGAL_assertion(pl1.has_on(s2.source()));
    CGAL_assertion(pl1.has_on(s2.target()));
    CGAL_assertion(pl1.has_on(s2.source()+Vector_3(1,0,0)));
    Object o1 = intersection(pl1,s1);
    if(!assign(ip1,o1))
      return false;
    CGAL_NEF_TRACEN("ip1 " << ip1);
    Vector_3 vec2(cross_product(s1.to_vector(),Vector_3(1,0,0)));
    Plane_3 pl2(s1.source(), vec2);
    CGAL_NEF_TRACEN("pl2 " << pl2);
    CGAL_assertion(pl2.has_on(s1.source()));
    CGAL_assertion(pl2.has_on(s1.target()));
    CGAL_assertion(pl2.has_on(s1.source()+Vector_3(1,0,0)));
    Object o2 = intersection(pl2,s2);
    if(!assign(ip2,o2) || ip2 == s2.source() || ip2 == s2.target())
      return false;

    CGAL_NEF_TRACEN("ips " << ip1 << ", " << ip2);
    if(compare_points(ip2, ip1))
      return true;
    return false;
  }
};

template<typename Nef_, typename FTComparison, typename Container>
  class Edge_sorter : public CGAL::Modifier_base<typename Nef_::SNC_and_PL> {

  typedef Nef_                                   Nef_polyhedron;
  typedef typename Nef_polyhedron::SNC_and_PL    SNC_and_PL;
  typedef typename Nef_polyhedron::SNC_structure SNC_structure;
  typedef typename Nef_polyhedron::Items         Items;
  typedef typename Nef_polyhedron::Kernel        Kernel;
  typedef CGAL::SNC_decorator<SNC_structure>     Base;
  typedef CGAL::SNC_point_locator<Base>          SNC_point_locator;
  typedef CGAL::SNC_constructor<Items, SNC_structure>
    SNC_constructor;

  typedef typename SNC_structure::Vertex_handle     Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle   SVertex_handle;
  typedef typename SNC_structure::Halfedge_iterator SVertex_iterator;
  typedef typename SNC_structure::Halfedge_handle   Halfedge_handle;
  typedef typename SNC_structure::SHalfedge_handle  SHalfedge_handle;
  typedef typename SNC_structure::Segment_3 Segment_3;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Vector_3 Vector_3;

  typedef typename Container::iterator Iterator;
  typedef typename Container::const_iterator Const_iterator;
  typedef typename Container::key_compare::Compare_points PointComparison;

  typedef typename CGAL::Insert_vertex_into_edge<SNC_structure, SNC_point_locator> InsertVertex;

  typedef typename CGAL::Generic_edge_sorter<
    Point_3,
    FTComparison,
    Need_to_split<Kernel, PointComparison>,
    InsertVertex,
    Container> GenericEdgeSorter;

  Container& c;
  SNC_structure* sncp;
  SNC_point_locator* pl;

 public:
  Edge_sorter(Container& cin) : c(cin) {}

  void operator()(SNC_and_PL& sncpl) {
    sncp = sncpl.sncp;
    pl = sncpl.pl;
    InsertVertex iv(*sncp, *pl);
    GenericEdgeSorter ges;
    Need_to_split<Kernel, PointComparison> nts;
    ges(c, nts, iv);
  }

  bool check() const {
    Const_iterator ci;
    for(ci = c.begin(); ci != c.end(); ++ci) {
      if(compare_points((*ci)->twin()->source()->point(),
                              (*ci)->source()->point())) {
        CGAL_error_msg("wrong orientation");
        return false;
      }
    }

    for(ci = c.begin(); ci != c.end(); ++ci) {
      Const_iterator ci2(ci);
      ++ci2;
      if(ci2 == c.end()) break;
      if(compare_points((*ci2)->source()->point(),
                        (*ci)->source()->point())) {
        CGAL_error_msg("wrong sorting");
        return false;
      }
    }

    for(ci = c.begin(); ci != c.end(); ++ci) {
      Const_iterator ci2(ci);
      ++ci2;
      while(ci2 != c.end() &&
            compare_fts((*ci)->twin()->source()->point().x(),
                        (*ci2)->source()->point().x())) {

        Point_3 ip;
        bool b = need_to_split(Segment_3((*ci)->source()->point(),
                                         (*ci)->twin()->source()->point()),
                               Segment_3((*ci2)->source()->point(),
                                         (*ci2)->twin()->source()->point()),ip);
        if(b) {
          CGAL_error_msg("vertical dependence");
          return false;
        }
        ++ci2;
      }
    }

    return true;
  }
};

} //namespace CGAL
#endif //CGAL_CD3_EDGE_SORTER_H
