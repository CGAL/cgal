// Copyright (c) 2005  Utrecht University
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
// $Id$ $Date$
// 
//
// Author(s)     : Ovidiu Grigore, Remco Veltkamp  <ovidiu.grigore@cs.uu.nl>



#ifndef CGAL_SIMPLIFY_POLYLINE_H
#define CGAL_SIMPLIFY_POLYLINE_H

#include <CGAL/polyap_fct.h>
#include <CGAL/polyap_traits.h>

CGAL_BEGIN_NAMESPACE

template <class Method, class ErrorAssessment, class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_number_of_points {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t n_pt_bound, 
			    typename DistTraits::FT &error, 
			    OutputIterator result)
  {
   std::cout << "generic" << std::endl; 
   return result;
  }
};

template <class Method, class ErrorAssessment, class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_error {
  
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t& n_pt_bound, 
			    typename DistTraits::FT error, 
			    OutputIterator result)
  {
   std::cout << "generic" << std::endl; 
   return result;
  }
};


struct Local{};
struct Global{};

struct Dynamic_programmming {};
struct Graph_search {};
struct Convex_hull_graph_search {};
struct Iterative_graph_search {};
struct Recursive_split {};

////////////// Recursive Split 
template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_number_of_points<Recursive_split, Local, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t n_pt_bound, 
			    typename DistTraits::FT &error, 
			    OutputIterator result)
  {
   return polygonal_approximation_RS_bnp_lea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_number_of_points<Recursive_split, Global, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t n_pt_bound, 
			    typename DistTraits::FT &error, 
			    OutputIterator result)
  {
   return polygonal_approximation_RS_bnp_gea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_error<Recursive_split, Local, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t& n_pt_bound, 
			    typename DistTraits::FT error, 
			    OutputIterator result)
  {
   return polygonal_approximation_RS_be_lea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_error<Recursive_split, Global, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t& n_pt_bound, 
			    typename DistTraits::FT error, 
			    OutputIterator result)
  {
   return polygonal_approximation_RS_be_gea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};



//////////////  Dynamic Programmming
template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_number_of_points<Dynamic_programmming, Local, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t n_pt_bound, 
			    typename DistTraits::FT &error, 
			    OutputIterator result)
  {
   return polygonal_approximation_DP_bnp_lea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_number_of_points<Dynamic_programmming, Global, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t n_pt_bound, 
			    typename DistTraits::FT &error, 
			    OutputIterator result)
  {
   return polygonal_approximation_DP_bnp_gea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_error<Dynamic_programmming, Local, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t& n_pt_bound, 
			    typename DistTraits::FT error, 
			    OutputIterator result)
  {
   return polygonal_approximation_DP_be_lea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_error<Dynamic_programmming, Global, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t& n_pt_bound, 
			    typename DistTraits::FT error, 
			    OutputIterator result)
  {
   return polygonal_approximation_DP_be_gea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};


//////////////  Graph Search
template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_number_of_points<Graph_search, Local, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t n_pt_bound, 
			    typename DistTraits::FT &error, 
			    OutputIterator result)
  {
   return polygonal_approximation_GS_bnp_lea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_number_of_points<Graph_search, Global, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t n_pt_bound, 
			    typename DistTraits::FT &error, 
			    OutputIterator result)
  {
   return polygonal_approximation_GS_bnp_gea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_error<Graph_search, Local, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t& n_pt_bound, 
			    typename DistTraits::FT error, 
			    OutputIterator result)
  {
   return polygonal_approximation_GS_be_lea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_error<Graph_search, Global, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t& n_pt_bound, 
			    typename DistTraits::FT error, 
			    OutputIterator result)
  {
   return polygonal_approximation_GS_be_gea(begin, 
					     beyond, 
					     n_pt_bound, 
					     error, 
					     result,
					     DistTraits());
  }
};

//////////////  Convex Hull Graph Search (exists only for Local)
template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_number_of_points<Convex_hull_graph_search, Local, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t n_pt_bound, 
			    typename DistTraits::FT &error, 
			    OutputIterator result)
  {
   return polygonal_approximation_CHGS_bnp_lea(begin, 
					       beyond, 
					       n_pt_bound, 
					       error, 
					       result);
  }
};

template <class DistTraits, class InputIterator,class OutputIterator>
struct Simplify_polyline_bound_error<Convex_hull_graph_search, Local, DistTraits, InputIterator, OutputIterator> {
  static OutputIterator fct(InputIterator begin, 
			    InputIterator beyond, 
			    std::size_t& n_pt_bound, 
			    typename DistTraits::FT error, 
			    OutputIterator result)
  {
   return polygonal_approximation_CHGS_be_lea(begin, 
					      beyond, 
					      n_pt_bound, 
					      error, 
					      result);
  }
};
//////////////////////////////////////////////////////////////////////////////////////////////////////
template<class Method, class ErrorAssessment, class DistTraits, class InputIterator,class OutputIterator>
OutputIterator simplify_polyline_bound_number_of_points(InputIterator begin, 
                                                        InputIterator beyond, 
                                                        std::size_t n_pt_bound, 
                                                        typename DistTraits::FT &error, 
                                                        OutputIterator result)
{
  return Simplify_polyline_bound_number_of_points<Method,ErrorAssessment,DistTraits,InputIterator,OutputIterator>::fct
    (begin, 
     beyond, 
     n_pt_bound, 
     error, 
     result);
}

template<class Method, class ErrorAssessment, class DistTraits, class InputIterator,class OutputIterator>
OutputIterator simplify_polyline_bound_error(InputIterator begin, 
					     InputIterator beyond, 
					     std::size_t n_pt_bound, 
					     typename DistTraits::FT &error, 
					     OutputIterator result)
{
  return Simplify_polyline_bound_error<Method,ErrorAssessment,DistTraits,InputIterator,OutputIterator>::fct
    (begin, 
     beyond, 
     n_pt_bound, 
     error, 
     result);
}




CGAL_END_NAMESPACE

#endif // CGAL_SIMPLIFY_POLYLINE_H
