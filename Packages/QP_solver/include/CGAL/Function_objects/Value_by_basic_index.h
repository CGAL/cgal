#ifndef CGAL_FUNCTION_OBJECTS_VALUE_BY_BASIC_INDEX_H
#define CGAL_FUNCTION_OBJECTS_VALUE_BY_BASIC_INDEX_H

#include <functional>
#include <iterator>

CGAL_BEGIN_NAMESPACE

template < class RndAccIt >
class Value_by_basic_index : public std::unary_function<
    int, typename std::iterator_traits<RndAccIt>::value_type > {

  public:
    typedef typename
    std::unary_function<
      int, typename std::iterator_traits
         <RndAccIt>::value_type >::result_type
    result_type;

    Value_by_basic_index( RndAccIt x_B_O_it,            int n_original,
			  RndAccIt x_B_S_it = x_B_O_it, int n_slack = 0)
	: o( x_B_O_it), s( x_B_S_it),
	  l( n_original), u( n_original+n_slack),
	  z( 0)
	{ }

    result_type  operator () ( int i) const
        {
	    if ( i < 0) return z;
	    if ( i >= l && i < u) return s[ i];
	    return o[ i];
	}

  private:
    RndAccIt     o, s;
    int          l, u;
    result_type  z;
};

CGAL_END_NAMESPACE
  
#endif
