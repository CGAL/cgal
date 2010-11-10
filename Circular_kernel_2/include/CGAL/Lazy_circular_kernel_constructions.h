// Copyright (c) 2003-2008  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s)     : Monique Teillaud, Sylvain Pion, Andreas Fabri, Constantinos Tsirogiannis

// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_LAZY_CIRCULAR_KERNEL_CONSTRUCTIONS_H
#define CGAL_LAZY_CIRCULAR_KERNEL_CONSTRUCTIONS_H

#include <CGAL/Lazy.h>
#include <CGAL/Profile_counter.h>

//TODO : More if-else's with object cast for all the possible 
//       types that could be returned by object-returning functions

namespace CGAL {

template <typename LK>
Object
make_lazy_CK(const Object& eto)
{
  typedef typename LK::AK AK;
  typedef typename LK::EK EK;
  typedef typename LK::E2A E2A;

  const std::pair<typename EK::Circular_arc_point_2, unsigned> *ptr;

  if((     ptr = 
           object_cast<std::pair<typename EK::Circular_arc_point_2,
	                         unsigned > >(&eto))){
           return make_object(std::make_pair(typename LK::Circular_arc_point_2(
           typename LK::Circular_arc_point_2::Rep(
	   new Lazy_rep_0<typename AK::Circular_arc_point_2, 
	   typename EK::Circular_arc_point_2, E2A>(ptr->first))),ptr->second));
  } else if(const typename EK::Circular_arc_2* ptr = 
            object_cast<typename EK::Circular_arc_2>(&eto)){
            return make_object(typename LK::Circular_arc_2(typename LK::Circular_arc_2::Rep(
	    new Lazy_rep_0<typename AK::Circular_arc_2, 
	    typename EK::Circular_arc_2, E2A>(*ptr))));
  } else if(const typename EK::Line_arc_2* ptr = 
            object_cast<typename EK::Line_arc_2>(&eto)){
            return make_object(typename LK::Line_arc_2(typename LK::Line_arc_2::Rep(
	    new Lazy_rep_0<typename AK::Line_arc_2, 
	    typename EK::Line_arc_2, E2A>(*ptr))));
  } else{
    std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's" << std::endl;
  }            
  return Object();
}


// This is the magic functor for functors that write their result as Objects into an output iterator

template <typename LK, typename AC, typename EC>
struct Lazy_construct_intersections_2 {

  static const bool Protection = true;

  typedef typename LK::AK AK;
  typedef typename LK::EK EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;
  typedef void result_type;
  typedef Lazy<Object, Object, EFT, E2A> Lazy_object;
  typedef Lazy<std::vector<Object>, std::vector<Object>, EFT, E2A> Lazy_vector;
  AC ac;
  EC ec;

public:

  // In the example we intersect two Lazy<Segment>s
  // and write into a back_inserter(list<Object([Lazy<Point>,Lazy<Segment>]) >)
  template <typename L1, typename L2, typename OutputIterator>
  OutputIterator
  operator()(const L1& l1, const L2& l2, OutputIterator it) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    try {
      Protect_FPU_rounding<Protection> P;
      Lazy_vector lv(new Lazy_rep_with_vector_2<AC, EC, E2A, L1, L2>(ac, ec, l1, l2));

      for(unsigned int i = 0; i < lv.approx().size(); i++){
        
	const std::pair<typename AK::Circular_arc_point_2, unsigned> *temp_p;
      
	if((temp_p=object_cast<std::pair<typename AK::Circular_arc_point_2,
	                         unsigned > >(& (lv.approx()[i])))){
	  *it = make_object(std::make_pair(typename LK::Circular_arc_point_2(
                             typename LK::Circular_arc_point_2::Rep(
	  new Lazy_rep_1<Ith<typename AK::Circular_arc_point_2>, 
	  Ith<typename EK::Circular_arc_point_2>, E2A, 
	  Lazy_vector>(Ith<typename AK::Circular_arc_point_2>(i,false), 
	  Ith<typename EK::Circular_arc_point_2>(i,false), lv))),temp_p->second)); 
	  ++it;
	} else if(object_cast<typename AK::Circular_arc_2>(& (lv.approx()[i]))){
	  *it = make_object(typename LK::Circular_arc_2( typename LK::Circular_arc_2::Rep(
	  new Lazy_rep_1<Ith<typename AK::Circular_arc_2>, 
	  Ith<typename EK::Circular_arc_2>,E2A, 
	  Lazy_vector>(Ith<typename AK::Circular_arc_2>(i), 
	  Ith<typename EK::Circular_arc_2>(i), lv))));
	  ++it;
	} else if(object_cast<typename AK::Line_arc_2>(& (lv.approx()[i]))){
	  *it = make_object(typename LK::Line_arc_2( typename LK::Line_arc_2::Rep(
	  new Lazy_rep_1<Ith<typename AK::Line_arc_2>, 
	  Ith<typename EK::Line_arc_2>,E2A, 
	  Lazy_vector>(Ith<typename AK::Line_arc_2>(i), 
	  Ith<typename EK::Line_arc_2>(i), lv))));
	  ++it;
	} else{
	       std::cout << "UNEXPECTED CONSTRUCT_INTERSECTIONS_2 PRODUCT" << std::endl;
	       std::cout << lv.approx()[i].type().name() << std::endl;
	       CGAL_error();
	}
      }
      
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
      // TODO: Instead of using a vector, write an iterator adapter
      std::vector<Object> exact_objects;
      ec(CGAL::exact(l1), CGAL::exact(l2), std::back_inserter(exact_objects));
      for (std::vector<Object>::iterator oit = exact_objects.begin();
	   oit != exact_objects.end();
	   oit++){
	*it = make_lazy_CK<LK>(*oit);
	++it;
      }
    }
    return it;
  }
};

template <typename LK, typename AC, typename EC>
struct Lazy_make_x_monotone_2 {

  static const bool Protection = true;

  typedef typename LK::AK AK;
  typedef typename LK::EK EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;
  typedef void result_type;
  typedef Lazy<Object, Object, EFT, E2A> Lazy_object;
  typedef Lazy<std::vector<Object>, std::vector<Object>, EFT, E2A> Lazy_vector;
  AC ac;
  EC ec;

public:

  // In the example we intersect two Lazy<Segment>s
  // and write into a back_inserter(list<Object([Lazy<Point>,Lazy<Segment>]) >)
  template <typename L1, typename OutputIterator>
  OutputIterator
  operator()(const L1& l1,OutputIterator it) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    try {
      Protect_FPU_rounding<Protection> P;
      Lazy_vector lv(new Lazy_rep_with_vector_1<AC, EC, E2A, L1>(ac, ec, l1));
      // lv.approx() is a std::vector<Object([AK::Point_2,AK::Segment_2])>
      // that is, when we get here we have constructed all approximate results
      for(unsigned int i = 0; i < lv.approx().size(); i++){
          if(object_cast<typename AK::Circular_arc_2>(& (lv.approx()[i]))){
	  *it = make_object(typename LK::Circular_arc_2(
          typename LK::Circular_arc_2::Rep(  
	  new Lazy_rep_1<Ith<typename AK::Circular_arc_2>, 
	  Ith<typename EK::Circular_arc_2>,E2A, 
	  Lazy_vector>(Ith<typename AK::Circular_arc_2>(i), 
	  Ith<typename EK::Circular_arc_2>(i), lv))));
	  ++it;
	} else if(object_cast<typename AK::Line_arc_2>(& (lv.approx()[i]))){
	  *it = make_object(typename LK::Line_arc_2(typename LK::Line_arc_2::Rep(
	  new Lazy_rep_1<Ith<typename AK::Line_arc_2>, 
	  Ith<typename EK::Line_arc_2>,E2A, 
	  Lazy_vector>(Ith<typename AK::Line_arc_2>(i), 
	  Ith<typename EK::Line_arc_2>(i), lv))));
	  ++it;
	} 
	else {
	  CGAL_error_msg( "UNEXPECTED MAKE_X_MONOTONE PRODUCT");
	}
      }
      
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
      // TODO: Instead of using a vector, write an iterator adapter
      std::vector<Object> exact_objects;
      ec(CGAL::exact(l1), std::back_inserter(exact_objects));
      for (std::vector<Object>::iterator oit = exact_objects.begin();
	   oit != exact_objects.end();
	   oit++){
	*it = make_lazy_CK<LK>(*oit);
	++it;
      }
    }
    return it;
  }
};


template <typename LK, typename AC, typename EC>
struct Lazy_advanced_make_x_monotone_2 {

  static const bool Protection = true;

  typedef typename LK::AK AK;
  typedef typename LK::EK EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;
  typedef void result_type;
  typedef Lazy<Object, Object, EFT, E2A> Lazy_object;
  typedef Lazy<std::vector<Object>, std::vector<Object>, EFT, E2A> Lazy_vector;
  AC ac;
  EC ec;

public:

  // In the example we intersect two Lazy<Segment>s
  // and write into a back_inserter(list<Object([Lazy<Point>,Lazy<Segment>]) >)
  template <typename L1, typename OutputIterator>
  OutputIterator
  operator()(const L1& l1,OutputIterator it) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    try {
      Protect_FPU_rounding<Protection> P;
      Lazy_vector lv(new Lazy_rep_with_vector_1<AC, EC, E2A, L1>(ac, ec, l1));
      // lv.approx() is a std::vector<Object([AK::Point_2,AK::Segment_2])>
      // that is, when we get here we have constructed all approximate results
      for(unsigned int i = 0; i < lv.approx().size(); i++){
          if(object_cast<typename AK::Circular_arc_2>(& (lv.approx()[i].first))){
	  *it = std::make_pair(make_object(typename LK::Circular_arc_2(typename LK::Circular_arc_2::Rep(
	  new Lazy_rep_1<Ith<typename AK::Circular_arc_2>, 
	  Ith<typename EK::Circular_arc_2>,E2A, 
	  Lazy_vector>(Ith<typename AK::Circular_arc_2>(i), 
	  Ith<typename EK::Circular_arc_2>(i), lv)))),lv.approx()[i].second);
	  ++it;
	} else if(object_cast<typename AK::Line_arc_2>(& (lv.approx()[i].first))){
	  *it = std::make_pair(make_object(typename LK::Line_arc_2(typename LK::Line_arc_2::Rep(
	  new Lazy_rep_1<Ith<typename AK::Line_arc_2>, 
	  Ith<typename EK::Line_arc_2>,E2A, 
	  Lazy_vector>(Ith<typename AK::Line_arc_2>(i), 
	  Ith<typename EK::Line_arc_2>(i), lv)))),lv.approx()[i].second);
	  ++it;
	} 
	else{
	  CGAL_error_msg( "UNEXPECTED ADVANCED_MAKE_X_MONOTONE PRODUCT");
	}
      }
      
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
      // TODO: Instead of using a vector, write an iterator adapter
      std::vector<Object> exact_objects;
      ec(CGAL::exact(l1), std::back_inserter(exact_objects));
      for (std::vector<Object>::iterator oit = exact_objects.begin();
	   oit != exact_objects.end();
	   ++oit){
        CGAL_error_msg( "Unfinished code !!!");
	//*it = std::make_pair(make_lazy_CK<LK>((*oit).first),(*oit).second);
	++it;
      }
    }
    return it;
  }
};

// The following functor returns an Object with a Lazy<Something> inside
// As the nested kernels return Objects of AK::Something and EK::Something
// we have to unwrap them from the Object, and wrap them in a Lazy<Something>
//
// TODO: write operators for more than two arguments. For the current kernel we only need two for       
//    Construct_intersections_2 and one for Make_x_monotone_2

template <typename LK, typename AC, typename EC>
struct Lazy_construction_object_CK {

  static const bool Protection = true;
  typedef typename LK::AK AK;
  typedef typename LK::EK EK;
  typedef typename EK::FT EFT;
  typedef typename LK::E2A E2A;
  typedef typename AC::result_type AT;
  typedef typename EC::result_type ET;
  typedef Object result_type;

  typedef Lazy<Object, Object, EFT, E2A> Lazy_object;
  AC ac;
  EC ec;

public:


template <typename L1>
  result_type
  operator()(const L1& l1) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    try {
      Protect_FPU_rounding<Protection> P;
      Lazy_object lo(new Lazy_rep_1<AC, EC, E2A, L1>(ac, ec, l1));
      
       std::pair<typename AK::Circular_arc_point_2, unsigned> *temp_p;
 
         if((temp_p=object_cast<std::pair<typename AK::Circular_arc_point_2, unsigned > >(& (lo.approx())))){
	typedef Lazy_rep_1<Object_cast<typename AK::Circular_arc_point_2>, 
	        Object_cast<typename EK::Circular_arc_point_2>, E2A, Lazy_object> Lcr;
	        Lcr * lcr = new Lcr(Object_cast<typename AK::Circular_arc_point_2>(), 
	        Object_cast<typename EK::Circular_arc_point_2>(), lo); 
	        return make_object(std::make_pair(typename LK::Circular_arc_point_2(
                  typename LK::Circular_arc_point_2::Rep(lcr)),temp_p->second));
      } else if(object_cast<typename AK::Circular_arc_point_2>(& (lo.approx()))){
	typedef Lazy_rep_1<Object_cast<typename AK::Circular_arc_point_2>, 
	        Object_cast<typename EK::Circular_arc_point_2>, E2A, Lazy_object> Lcr;
	        Lcr * lcr = new Lcr(Object_cast<typename AK::Circular_arc_point_2>(), 
	        Object_cast<typename EK::Circular_arc_point_2>(), lo); 
	        return make_object(typename LK::Circular_arc_point_2(typename LK::Circular_arc_point_2::Rep(lcr)));
      } else if(object_cast<typename AK::Circular_arc_2>(& (lo.approx()))){
	typedef Lazy_rep_1<Object_cast<typename AK::Circular_arc_2>, 
	        Object_cast<typename EK::Circular_arc_2>, E2A, Lazy_object> Lcr;
	        Lcr * lcr = new Lcr(Object_cast<typename AK::Circular_arc_2>(), 
		Object_cast<typename EK::Circular_arc_2>(), lo); 
	        return make_object(typename LK::Circular_arc_2(typename LK::Circular_arc_2::Rep(lcr)));}
	else if(object_cast<typename AK::Line_arc_2>(& (lo.approx()))){
	typedef Lazy_rep_1<Object_cast<typename AK::Line_arc_2>, 
	        Object_cast<typename EK::Line_arc_2>, E2A, Lazy_object> Lcr;
	        Lcr * lcr = new Lcr(Object_cast<typename AK::Line_arc_2>(), 
		Object_cast<typename EK::Line_arc_2>(), lo); 
	        return make_object(typename LK::Line_arc_2(typename LK::Line_arc_2::Rep(lcr)));}
        else {
	std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's" << std::endl;
      }
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
      ET eto = ec(CGAL::exact(l1));
      return make_lazy_CK<LK>(eto);
    }
    return Object();
  }

  
  
  template <typename L1, typename L2>
  result_type
  operator()(const L1& l1, const L2& l2) const
  {
    CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);
    try {
      Protect_FPU_rounding<Protection> P;
      Lazy_object lo(new Lazy_rep_2<AC, EC, E2A, L1, L2>(ac, ec, l1, l2));
 
      std::pair<typename AK::Circular_arc_point_2, unsigned> *temp_p;    
      
      
    if((temp_p=object_cast<std::pair<typename AK::Circular_arc_point_2, unsigned > >(& (lo.approx())))){
	typedef Lazy_rep_1<Object_cast<typename AK::Circular_arc_point_2>, 
	        Object_cast<typename EK::Circular_arc_point_2>, E2A, Lazy_object> Lcr;
	        Lcr * lcr = new Lcr(Object_cast<typename AK::Circular_arc_point_2>(), 
	        Object_cast<typename EK::Circular_arc_point_2>(), lo); 
	        return make_object(std::make_pair(typename LK::Circular_arc_point_2(
                  typename LK::Circular_arc_point_2::Rep(lcr)),temp_p->second));
      } else if(object_cast<typename AK::Circular_arc_point_2>(& (lo.approx()))){
	typedef Lazy_rep_1<Object_cast<typename AK::Circular_arc_point_2>, 
	        Object_cast<typename EK::Circular_arc_point_2>, E2A, Lazy_object> Lcr;
	        Lcr * lcr = new Lcr(Object_cast<typename AK::Circular_arc_point_2>(), 
	        Object_cast<typename EK::Circular_arc_point_2>(), lo); 
	        return make_object(typename LK::Circular_arc_point_2(typename LK::Circular_arc_point_2::Rep(lcr)));
      } else if(object_cast<typename AK::Circular_arc_2>(& (lo.approx()))){
	typedef Lazy_rep_1<Object_cast<typename AK::Circular_arc_2>, 
	        Object_cast<typename EK::Circular_arc_2>, E2A, Lazy_object> Lcr;
	        Lcr * lcr = new Lcr(Object_cast<typename AK::Circular_arc_2>(), 
		Object_cast<typename EK::Circular_arc_2>(), lo); 
	        return make_object(typename LK::Circular_arc_2(typename LK::Circular_arc_2::Rep(lcr)));}
	else if(object_cast<typename AK::Line_arc_2>(& (lo.approx()))){
	typedef Lazy_rep_1<Object_cast<typename AK::Line_arc_2>, 
	        Object_cast<typename EK::Line_arc_2>, E2A, Lazy_object> Lcr;
	        Lcr * lcr = new Lcr(Object_cast<typename AK::Line_arc_2>(), 
		Object_cast<typename EK::Line_arc_2>(), lo); 
	        return make_object(typename LK::Line_arc_2(typename LK::Line_arc_2::Rep(lcr)));}
        else {
	std::cerr << "object_cast inside Lazy_construction_rep::operator() failed. It needs more else if's" << std::endl;
      }
    } catch (Uncertain_conversion_exception) {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);
      Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
      ET eto = ec(CGAL::exact(l1), CGAL::exact(l2));
      return make_lazy_CK<LK>(eto);
    }
    return Object();
  }
  
};

} //namespace CGAL

#endif // CGAL_LAZY_CIRCULAR_KERNEL_CONSTRUCTIONS_H
