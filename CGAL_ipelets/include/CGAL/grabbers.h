// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot, Sylvain Pion


#ifndef CGAL_GRABBER_H
#define CGAL_GRABBER_H

#include <boost/function_output_iterator.hpp>

namespace CGAL{

template <class Kernel, class Container>
class Polygon_2;


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
  typedef typename Kernel::Weighted_point_2 Self;
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
