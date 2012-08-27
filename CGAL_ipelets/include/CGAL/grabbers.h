// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//
// Author(s)     : Sebastien Loriot, Sylvain Pion


#ifndef CGAL_GRABBER_H
#define CGAL_GRABBER_H

#include <boost/function_output_iterator.hpp>

namespace CGAL{
  
namespace internal{

template <class Kernel, class output_iterator>
class Point_grabber{
  output_iterator out;
public:  
  Point_grabber(output_iterator it):out(it){}
  
  void operator()(const typename Kernel::Point_2& p){
    *out++=p;
  }
    
  void operator()(const typename Kernel::Segment_2& s){
    *out++=s[0];
    *out++=s[1];
  }
  
  template<class Container>
  void operator()(const CGAL::Polygon_2<Kernel,Container>& p){
    for(typename CGAL::Polygon_2<Kernel,Container>::Vertex_iterator it=
        p.vertices_begin();it!=p.vertices_end();++it)
      *out++= *it;
  }
};

template<class Kernel,class output_iterator>
boost::function_output_iterator<Point_grabber<Kernel,output_iterator> >
point_grabber(output_iterator it){
  return boost::make_function_output_iterator(Point_grabber<Kernel,output_iterator>(it));
}


//Segments
template <class Kernel, class output_iterator>
class Segment_grabber{
  output_iterator out;
public:  
  Segment_grabber(output_iterator it):out(it){}
  
  void operator()(const typename Kernel::Segment_2& s){
    *out++=s;
  }
  
  template<class Container>
  void operator()(const CGAL::Polygon_2<Kernel,Container>& p){
    for(typename CGAL::Polygon_2<Kernel,Container>::Edge_const_iterator
        it=p.edges_begin();it!=p.edges_end();++it)
      *out++= *it;
  }
};


template<class Kernel,class output_iterator>
boost::function_output_iterator<Segment_grabber<Kernel,output_iterator> >
segment_grabber(output_iterator it){
  return boost::make_function_output_iterator(Segment_grabber<Kernel,output_iterator>(it));
}


//Weighted points
template <class Kernel,class output_iterator>
class Wpoint_grabber{
  output_iterator out;
  typedef CGAL::Weighted_point<typename Kernel::Point_2,typename Kernel::FT> Self;
public:  
  Wpoint_grabber(output_iterator it):out(it){}

  void operator()(const Self& p){
    *out++=p;
  }
    
  void operator()(const typename Kernel::Point_2& p){
    *out++=Self(p,0);
  }
    
  void operator()(const typename Kernel::Circle_2& c){
    *out++=Self(c.center(),c.squared_radius());
  }
  
  void operator()(const typename Kernel::Segment_2& s){
    *out++=Self(s[0],0);
    *out++=Self(s[1],0);
  }  
  
  template<class Container>
  void operator()(const CGAL::Polygon_2<Kernel,Container>& p){
    for(typename CGAL::Polygon_2<Kernel,Container>::Vertex_iterator 
        it=p.vertices_begin();it!=p.vertices_end();++it)
      *out++= Self(*it,0);
  }
};

template<class Kernel,class output_iterator>
boost::function_output_iterator<Wpoint_grabber<Kernel,output_iterator> >
wpoint_grabber(output_iterator it){
  return boost::make_function_output_iterator(Wpoint_grabber<Kernel,output_iterator>(it));
}

}//internal
}//CGAL
#endif //CGAL_GRABBER_H
