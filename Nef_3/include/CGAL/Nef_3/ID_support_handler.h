// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: 
// $Id: 
// 
//
// Author(s)     :     Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_ID_SUPPORT_HANDLER
#define CGAL_ID_SUPPORT_HANDLER

#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Unique_hash_map.h>

CGAL_BEGIN_NAMESPACE

template<typename Items, typename Decorator>
class ID_support_handler {

  typedef typename Decorator::SVertex_handle SVertex_handle;
  typedef typename Decorator::SHalfedge_handle SHalfedge_handle;
  
  typedef typename Decorator::SVertex_const_handle SVertex_const_handle;
  typedef typename Decorator::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename Decorator::SHalfloop_const_handle SHalfloop_const_handle;
  
 public:
  ID_support_handler() {}
  
  void handle_support(SVertex_handle , 
		      SHalfedge_const_handle ,
		      SHalfedge_const_handle ) {} 
  
  void handle_support(SVertex_handle ,
		      SHalfloop_const_handle ,
		      SHalfloop_const_handle ) {}
  
  void handle_support(SVertex_handle,
		      SHalfloop_const_handle,
		      SHalfedge_const_handle) {} 
  
  void handle_support(SVertex_handle,
		      SHalfedge_const_handle,
		      SHalfloop_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SHalfedge_const_handle,
		      SVertex_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SVertex_const_handle,
		      SHalfedge_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SVertex_const_handle,
		      SVertex_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SVertex_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SVertex_const_handle,
		      SHalfloop_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SHalfloop_const_handle,
		      SVertex_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfedge_const_handle,
		      SHalfedge_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfedge_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfloop_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfedge_const_handle,
		      SHalfloop_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfloop_const_handle,
		      SHalfedge_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfloop_const_handle,
		      SHalfloop_const_handle) {}
};

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
  
 public:
  ID_support_handler() {}
  
  void hash_facet_pair(SVertex_handle sv, 
		       Halffacet_const_handle f1, 
		       Halffacet_const_handle f2) {
    //      std::cerr << "hash_facet_pair " << sv->point() << std::endl
    //		<< "  " << f1->plane() << &f1 << std::endl
    //		<< " "  << f2->plane() << &f2 << std::endl;
    
    if(f2m[f1][f2]==0) {
      sv->set_index();
      f2m[f1][f2] = sv->get_index();
      //	std::cerr << "insert " << sv->point() << &*sv 
      //		  << ": " << f2m[f1][f2] << std::endl;
      //	std::cerr << "not defined, yet" << std::endl;
    }
    else {
      //	std::cerr << "access " << sv->point() << &*sv << std::endl;
      sv->set_index(f2m[f1][f2]);
      //	std::cerr << "combine " << sv->point() << "+" << f2e[f2]->point() << std::endl;
    }
  }
  
  void handle_support(SVertex_handle sv, 
		      SHalfedge_const_handle se1,
		      SHalfedge_const_handle se2) {
    //      std::cerr << "handle support ee " << sv->point() << std::endl;
    Halffacet_const_handle f1 = se1->facet();
    if(f1->is_twin()) f1 = f1->twin();
    Halffacet_const_handle f2 = se2->facet();
    if(f2->is_twin()) f2 = f2->twin();
    hash_facet_pair(sv, f1,f2);
  }
  
  void handle_support(SVertex_handle sv,
		      SHalfloop_const_handle sl1,
		      SHalfloop_const_handle sl2) {
    Halffacet_const_handle f1 = sl1->facet();
    if(f1->is_twin()) f1 = f1->twin();
    Halffacet_const_handle f2 = sl2->facet();
    if(f2->is_twin()) f2 = f2->twin();
    hash_facet_pair(sv, f1,f2);
  }
  
  void handle_support(SVertex_handle sv,
		      SHalfloop_const_handle sl1,
		      SHalfedge_const_handle se2, 
		      bool inverse_order = false) {
    //      std::cerr << "handle support el " << sv->point() << std::endl;
    Halffacet_const_handle f1 = sl1->facet();
    if(f1->is_twin()) f1 = f1->twin();
    Halffacet_const_handle f2 = se2->facet();
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
    //      std::cerr << "handle support ev " << sv->point() << std::endl;
    sv->set_index(sv2->get_index());
  }
  
  void handle_support(SVertex_handle sv,
		      SVertex_const_handle sv1,
		      SHalfedge_const_handle se2) {
    //      std::cerr << "handle support ve " << sv->point() << std::endl;
    handle_support(sv,se2,sv1);
  }
  
  void handle_support(SVertex_handle sv,
		      SVertex_const_handle sv1,
		      SVertex_const_handle sv2) {
    //      std::cerr << "handle support vv " << sv->point() << std::endl;
    //      std::cerr << "  supported by " << sv1->point() << std::endl;
    //      std::cerr << "  supported by " << sv2->point() << std::endl;
    if(sv1->get_index() < sv2->get_index())
      sv->set_index(sv1->get_index());
    else
      sv->set_index(sv2->get_index());
  }
  
  void handle_support(SVertex_handle sv,
		      SVertex_const_handle sv0) {
    //      std::cerr << "handle support v " << sv->point() << std::endl;
    sv->set_index(sv0->get_index());
  }
  
  void handle_support(SVertex_handle sv,
		      SVertex_const_handle sv1,
		      SHalfloop_const_handle) {
    //      std::cerr << "handle support vl " << sv->point() << std::endl;
    sv->set_index(sv1->get_index());
  }
  
  void handle_support(SVertex_handle sv,
		      SHalfloop_const_handle sl1,
		      SVertex_const_handle sv2) {
    //      std::cerr << "handle support lv " << sv->point() << std::endl;
    handle_support(sv, sv2, sl1);
  }
  
  void handle_support(SHalfedge_handle se,
		      SHalfedge_const_handle se1) {
    CGAL_assertion(se->circle() == se1->circle());
    se->set_index(se1->get_index());
    se->twin()->set_index(se1->twin()->get_index());
    //      std::cerr << "se " << se->source()->point()
    //		<< "->" << se->twin()->source()->point() 
    //		<< "|" << se->circle() << std::endl;
    //      std::cerr << "se1 " << se1->get_index() << std::endl;
    //      std::cerr << "se1->twin() " << se1->twin()->get_index() << std::endl;
    //      std::cerr << "result " << se->get_index() 
    //		<< ", " << se->twin()->get_index() << std::endl;
  }
  
  void handle_support(SHalfedge_handle se,
		      SHalfloop_const_handle sl1) {
    CGAL_assertion(se->circle() == sl1->circle());
    se->set_index(sl1->get_index());
    se->twin()->set_index(sl1->twin()->get_index());
    //      std::cerr << "se " << se->source()->point()
    //		<< "->" << se->twin()->source()->point() 
    //		<< "|" << se->circle() << std::endl;
    //      std::cerr << "sl1 " << sl1->get_index() << std::endl;
    //      std::cerr << "sl1->twin() " << sl1->twin()->get_index() << std::endl;
    //      std::cerr << "result " << se->get_index() 
    //		<< ", " << se->twin()->get_index() << std::endl;
  }
  
  int set_se_index(int index1, int index2,
		   bool v1, bool v2) {
    if(v1) {
      if(v2) {
	return index1<index2 ? index1 : index2;
      } else {
	return index1;
      }
    }
    if(v2)
      return index2;
    return index1<index2 ? index1 : index2;
  }
  
  void handle_support(SHalfedge_handle se,
		      SHalfedge_const_handle se1,
		      SHalfedge_const_handle se2) {
    CGAL_assertion(se->circle() == se1->circle());
    CGAL_assertion(se->circle() == se2->circle());
    //      std::cerr << "se " << se->source()->point()
    //		<< "->" << se->twin()->source()->point() 
    //		<< "|" << se->circle() << std::endl;
    //      std::cerr << vs1 << vs2 << vt1 << vt2 << std::endl;
    //      std::cerr << "se1 " << se1->get_index() << std::endl;
    //      std::cerr << "se2 " << se2->get_index() << std::endl;
    //      std::cerr << "se1->twin() " << se1->twin()->get_index() << std::endl;
    //      std::cerr << "se2->twin() " << se2->twin()->get_index() << std::endl;
    /*
      se->set_index(set_se_index(se1->get_index(),
      se2->get_index(),
      vs1, vs2));
      se->twin()->set_index(set_se_index(se1->twin()->get_index(),
      se2->twin()->get_index(),
      vt1, vt2));
    */
    if(se1->get_index()<se2->get_index()) {
      //	if(se1->get_index() < hash[se2->get_index()]) 
      //	  hash[se2->get_index()] = se1->get_index();
      se->set_index(se1->get_index());
    } else {
      //	if(se2->get_index() < hash[se1->get_index()])
      //	  hash[se1->get_index()] = se2->get_index();
      se->set_index(se2->get_index());
    }
    if(se1->twin()->get_index()<se2->twin()->get_index()) {
      //	if(se1->twin()->get_index() < hash[se2->twin()->get_index()])
      //	  hash[se2->twin()->get_index()] = se1->twin()->get_index();
      se->twin()->set_index(se1->twin()->get_index());
    } else {
      //	if(se2->twin()->get_index() < hash[se1->twin()->get_index()])
      //	  hash[se1->twin()->get_index()] = se2->twin()->get_index();
      se->twin()->set_index(se2->twin()->get_index());
    }
    //      std::cerr << "result " << se->get_index() 
    //		<< ", " << se->twin()->get_index() << std::endl;
  }    
  
  void handle_support(SHalfedge_handle se,
		      SHalfedge_const_handle se1,
		      SHalfloop_const_handle sl2) {
    CGAL_assertion(se->circle() == se1->circle());
    CGAL_assertion(se->circle() == sl2->circle());
    if(se1->get_index()<sl2->get_index()) {
      //	if(se1->get_index() < hash[sl2->get_index()]) 
      //	  hash[sl2->get_index()] = se1->get_index();
      se->set_index(se1->get_index());
    } else {
      //	if(sl2->get_index() < hash[se1->get_index()])
      //	  hash[se1->get_index()] = sl2->get_index();
      se->set_index(sl2->get_index());
    }
    if(se1->twin()->get_index()<sl2->twin()->get_index()) {
      //	if(se1->twin()->get_index() < hash[sl2->twin()->get_index()])
      //	  hash[sl2->twin()->get_index()] = se1->twin()->get_index();
      se->twin()->set_index(se1->twin()->get_index());
    } else {
      //	if(sl2->twin()->get_index() < hash[se1->twin()->get_index()])
      //	  hash[se1->twin()->get_index()] = sl2->twin()->get_index();
      se->twin()->set_index(sl2->twin()->get_index());
    }
  }
  
  void handle_support(SHalfedge_handle se,
		      SHalfloop_const_handle sl1,
		      SHalfedge_const_handle se2) {
    handle_support(se,se2,sl1);
  }    
  
  void handle_support(SHalfedge_handle se,
		      SHalfloop_const_handle sl1,
		      SHalfloop_const_handle sl2) {
    CGAL_assertion(se->circle() == sl1->circle());
    CGAL_assertion(se->circle() == sl2->circle());
    if(sl1->get_index()<sl2->get_index()) {
      //	if(sl1->get_index() < hash[sl2->get_index()]) 
      //	  hash[sl2->get_index()] = sl1->get_index();
      se->set_index(sl1->get_index());
    } else {
      //	if(sl2->get_index() < hash[sl1->get_index()])
      //	  hash[sl1->get_index()] = sl2->get_index();
      se->set_index(sl2->get_index());
    }
    if(sl1->twin()->get_index()< sl2->twin()->get_index()) {
      //	if(sl1->twin()->get_index() < hash[sl2->twin()->get_index()])
      //	  hash[sl2->twin()->get_index()] = sl1->twin()->get_index();
      se->twin()->set_index(sl1->twin()->get_index());
    } else {
      //	if(sl2->twin()->get_index() < hash[sl1->twin()->get_index()])
      //	  hash[sl1->twin()->get_index()] = sl2->twin()->get_index();
      se->twin()->set_index(sl2->twin()->get_index());
    }
  } 
  
  //    int& hash_index(const int i) { return hash[i]; }    
};

CGAL_END_NAMESPACE
#endif // CGAL_ID_SUPPORT_HANDLER
