// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :     Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_3_ID_SUPPORT_HANDLER
#define CGAL_NEF_3_ID_SUPPORT_HANDLER

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_S2/ID_support_handler.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Unique_hash_map.h>
#include <map>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 131
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename Decorator>
class ID_support_handler<SNC_indexed_items, Decorator> {

  typedef typename Decorator::SVertex_handle SVertex_handle;
  typedef typename Decorator::SHalfedge_handle SHalfedge_handle;

  typedef typename Decorator::SVertex_const_handle SVertex_const_handle;
  typedef typename Decorator::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename Decorator::SHalfloop_const_handle SHalfloop_const_handle;
  typedef typename Decorator::Halffacet_const_handle Halffacet_const_handle;

  typedef CGAL::Unique_hash_map<Halffacet_const_handle, int> F2E;
  CGAL::Unique_hash_map<Halffacet_const_handle, F2E> f2m;
  std::map<int, int> hash;

 public:
  ID_support_handler() {}

  int get_hash(int i) {
    int root(i);
    while(hash[root] != root)
      root = hash[root];
    while(hash[i] != i) {
      int tmp = hash[i];
      hash[i] = root;
      i = tmp;
    }
    return root;
  }

  void set_hash(int i, int parent) {
    CGAL_assertion(parent <= i);
    hash[get_hash(i)] = parent;
  }

  template<typename Handle>
  void initialize_hash(Handle h) {
    hash[h->get_index()] = h->get_index();
  }
  void initialize_hash(int i) {
        hash[i] = i;
  }

  void hash_facet_pair(SVertex_handle sv,
                       Halffacet_const_handle f1,
                       Halffacet_const_handle f2) {
    CGAL_NEF_TRACEN("hash_facet_pair " << sv->point() << std::endl
                    << "  " << f1->plane() << &f1 << std::endl
                    << " "  << f2->plane() << &f2);

    if(f2m[f1][f2]==0) {
      sv->set_index();
      f2m[f1][f2] = sv->get_index();
      CGAL_NEF_TRACEN("insert " << sv->point() << &*sv
                      << ": " << f2m[f1][f2]);
      CGAL_NEF_TRACEN("not defined, yet");
    }
    else {
      CGAL_NEF_TRACEN("access " << sv->point() << &*sv);
      sv->set_index(f2m[f1][f2]);
    }
  }

  void handle_support(SVertex_handle sv,
                      SHalfedge_const_handle se1,
                      SHalfedge_const_handle se2) {
    CGAL_NEF_TRACEN("handle support ee " << sv->point());
    Halffacet_const_handle f1 = se1->get_index_facet();
    if(f1->is_twin()) f1 = f1->twin();
    Halffacet_const_handle f2 = se2->get_index_facet();
    if(f2->is_twin()) f2 = f2->twin();
    hash_facet_pair(sv, f1,f2);
  }

  void handle_support(SVertex_handle sv,
                      SHalfloop_const_handle sl1,
                      SHalfloop_const_handle sl2) {
    Halffacet_const_handle f1 = sl1->get_index_facet();
    if(f1->is_twin()) f1 = f1->twin();
    Halffacet_const_handle f2 = sl2->get_index_facet();
    if(f2->is_twin()) f2 = f2->twin();
    hash_facet_pair(sv, f1,f2);
  }

  void handle_support(SVertex_handle sv,
                      SHalfloop_const_handle sl1,
                      SHalfedge_const_handle se2,
                      bool inverse_order = false) {
    CGAL_NEF_TRACEN("handle support el " << sv->point());
    Halffacet_const_handle f1 = sl1->get_index_facet();
    if(f1->is_twin()) f1 = f1->twin();
    Halffacet_const_handle f2 = se2->get_index_facet();
    if(f2->is_twin()) f2 = f2->twin();
    if(inverse_order)
      hash_facet_pair(sv, f2,f1);
    else
      hash_facet_pair(sv, f1,f2);
  }

  void handle_support(SVertex_handle sv,
                      SHalfedge_const_handle se1,
                      SHalfloop_const_handle sl2) {
    handle_support(sv,sl2,se1,true);
  }

  void handle_support(SVertex_handle sv,
                      SHalfedge_const_handle,
                      SVertex_const_handle sv2) {
    CGAL_NEF_TRACEN("handle support ev " << sv->point());
    sv->set_index(sv2->get_index());
  }

  void handle_support(SVertex_handle sv,
                      SVertex_const_handle sv1,
                      SHalfedge_const_handle se2) {
    CGAL_NEF_TRACEN("handle support ve " << sv->point());
    handle_support(sv,se2,sv1);
  }

  void handle_support(SVertex_handle sv,
                      SVertex_const_handle sv1,
                      SVertex_const_handle sv2) {
    CGAL_NEF_TRACEN("handle support vv " << sv->point());
    CGAL_NEF_TRACEN("  supported by " << sv1->point());
    CGAL_NEF_TRACEN("  supported by " << sv2->point());
    if(sv1->get_index() < sv2->get_index())
      sv->set_index(sv1->get_index());
    else
      sv->set_index(sv2->get_index());
  }

  void handle_support(SVertex_handle sv,
                      SVertex_const_handle sv0) {
    CGAL_NEF_TRACEN("handle support v " << sv->point());
    sv->set_index(sv0->get_index());
  }

  void handle_support(SVertex_handle sv,
                      SVertex_const_handle sv1,
                      SHalfloop_const_handle) {
    CGAL_NEF_TRACEN("handle support vl " << sv->point());
    sv->set_index(sv1->get_index());
  }

  void handle_support(SVertex_handle sv,
                      SHalfloop_const_handle sl1,
                      SVertex_const_handle sv2) {
    CGAL_NEF_TRACEN("handle support lv " << sv->point());
    handle_support(sv, sv2, sl1);
  }

  void handle_support(SHalfedge_handle se,
                      SHalfedge_const_handle se1) {
    if(!equal_not_opposite(se->circle(), se1->circle()))
       se1 = se1->twin();
    se->set_index(se1->get_index());
    se->twin()->set_index(se1->twin()->get_index());
    CGAL_NEF_TRACEN("se " << se->source()->point()
                    << "->" << se->twin()->source()->point()
                    << "|" << se->circle());
    CGAL_NEF_TRACEN("se1 " << se1->get_index());
    CGAL_NEF_TRACEN("se1->twin() " << se1->twin()->get_index());
    CGAL_NEF_TRACEN("result " << se->get_index()
                    << ", " << se->twin()->get_index());
  }

  void handle_support(SHalfedge_handle se,
                      SHalfloop_const_handle sl1) {
    if(!equal_not_opposite(se->circle(), sl1->circle()))
       sl1 = sl1->twin();
    CGAL_assertion(se->circle() == sl1->circle());
    se->set_index(sl1->get_index());
    se->twin()->set_index(sl1->twin()->get_index());
    CGAL_NEF_TRACEN("se " << se->source()->point()
                    << "->" << se->twin()->source()->point()
                    << "|" << se->circle());
    CGAL_NEF_TRACEN("sl1 " << sl1->get_index());
    CGAL_NEF_TRACEN("sl1->twin() " << sl1->twin()->get_index());
    CGAL_NEF_TRACEN("result " << se->get_index()
                    << ", " << se->twin()->get_index());
  }

  void handle_support(SHalfedge_handle se,
                      SHalfedge_const_handle se1,
                      SHalfedge_const_handle se2) {
    if(!equal_not_opposite(se->circle(), se1->circle()))
      se1 = se1->twin();
    if(!equal_not_opposite(se->circle(), se2->circle()))
      se2 = se2->twin();
    CGAL_assertion(se->circle() == se1->circle());
    CGAL_assertion(se->circle() == se2->circle());
    CGAL_NEF_TRACEN("se " << se->source()->point()
                    << "->" << se->twin()->source()->point()
                    << "|" << se->circle());
    CGAL_NEF_TRACEN("se1 " << se1->get_index());
    CGAL_NEF_TRACEN("se2 " << se2->get_index());
    CGAL_NEF_TRACEN("se1->twin() " << se1->twin()->get_index());
    CGAL_NEF_TRACEN("se2->twin() " << se2->twin()->get_index());

    int index1 = get_hash(se1->get_index());
    int index2 = get_hash(se2->get_index());
    if(index1 < index2) {
      se->set_index(index1);
      set_hash(se2->get_index(), index1);
    } else {
      se->set_index(index2);
      set_hash(se1->get_index(), index2);
    }

    index1 = get_hash(se1->twin()->get_index());
    index2 = get_hash(se2->twin()->get_index());
    if(index1 < index2) {
      se->twin()->set_index(index1);
      set_hash(se1->twin()->get_index(), index1);
    } else {
      se->twin()->set_index(index2);
      set_hash(se1->twin()->get_index(), index2);
    }

    CGAL_NEF_TRACEN("result " << se->get_index()
                    << ", " << se->twin()->get_index());
  }

  void handle_support(SHalfedge_handle se,
                      SHalfedge_const_handle se1,
                      SHalfloop_const_handle sl2) {

    if(!equal_not_opposite(se->circle(), se1->circle()))
      se1 = se1->twin();
    if(!equal_not_opposite(se->circle(), sl2->circle()))
      sl2 = sl2->twin();
    CGAL_assertion(se->circle() == se1->circle());
    CGAL_assertion(se->circle() == sl2->circle());
    CGAL_NEF_TRACEN("se " << se->source()->point()
                    << "->" << se->twin()->source()->point()
                    << "|" << se->circle());
    CGAL_NEF_TRACEN("se1 " << se1->get_index());
    CGAL_NEF_TRACEN("sl2 " << sl2->get_index());
    CGAL_NEF_TRACEN("se1->twin() " << se1->twin()->get_index());
    CGAL_NEF_TRACEN("sl2->twin() " << sl2->twin()->get_index());

    int index1 = get_hash(se1->get_index());
    int index2 = get_hash(sl2->get_index());
    if(index1 < index2) {
      se->set_index(index1);
      set_hash(se1->get_index(), index1);
      set_hash(sl2->get_index(), index1);
    } else {
      se->set_index(index2);
      set_hash(se1->get_index(), index2);
      set_hash(sl2->get_index(), index2);
    }

    index1 = get_hash(se1->twin()->get_index());
    index2 = get_hash(sl2->twin()->get_index());
    if(index1 < index2) {
      se->twin()->set_index(index1);
      set_hash(sl2->twin()->get_index(), index1);
    } else {
      se->twin()->set_index(index2);
      set_hash(se1->twin()->get_index(), index2);
    }

    CGAL_NEF_TRACEN("result " << se->get_index()
                    << ", " << se->twin()->get_index());
  }

  void handle_support(SHalfedge_handle se,
                      SHalfloop_const_handle sl1,
                      SHalfedge_const_handle se2) {
    handle_support(se,se2,sl1);
  }

  void handle_support(SHalfedge_handle se,
                      SHalfloop_const_handle sl1,
                      SHalfloop_const_handle sl2) {
    if(!equal_not_opposite(se->circle(), sl1->circle()))
      sl1 = sl1->twin();
    if(!equal_not_opposite(se->circle(), sl2->circle()))
      sl2 = sl2->twin();
    CGAL_assertion(se->circle() == sl1->circle());
    CGAL_assertion(se->circle() == sl2->circle());
    CGAL_NEF_TRACEN("se " << se->source()->point()
                    << "->" << se->twin()->source()->point()
                    << "|" << se->circle());
    CGAL_NEF_TRACEN("sl1 " << sl1->get_index());
    CGAL_NEF_TRACEN("sl2 " << sl2->get_index());
    CGAL_NEF_TRACEN("sl1->twin() " << sl1->twin()->get_index());
    CGAL_NEF_TRACEN("sl2->twin() " << sl2->twin()->get_index());

    int index1 = get_hash(sl1->get_index());
    int index2 = get_hash(sl2->get_index());
    if(index1 < index2) {
      se->set_index(index1);
      set_hash(sl2->get_index(), index1);
    } else {
      se->set_index(index2);
      set_hash(sl1->get_index(), index2);
    }

    index1 = get_hash(sl1->twin()->get_index());
    index2 = get_hash(sl2->twin()->get_index());
    if(index1 < index2) {
      se->twin()->set_index(index1);
      set_hash(sl2->twin()->get_index(), index1);
    } else {
      se->twin()->set_index(index2);
      set_hash(sl1->twin()->get_index(), index2);
    }

    CGAL_NEF_TRACEN("result " << se->get_index()
                    << ", " << se->twin()->get_index());
  }

};

} //namespace CGAL
#endif // CGAL_NEF_3_ID_SUPPORT_HANDLER
