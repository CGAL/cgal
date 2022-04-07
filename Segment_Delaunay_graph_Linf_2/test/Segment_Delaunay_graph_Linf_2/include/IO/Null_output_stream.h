#ifndef CGAL_NULL_OUTPUT_STREAM_H
#define CGAL_NULL_OUTPUT_STREAM_H

#include <CGAL/basic.h>


namespace CGAL {


struct Null_output_stream {};

#if defined(__INTEL_COMPILER)
template<class T>
Null_output_stream&
operator<<(Null_output_stream& nos, const T&);
#endif

template<class T>
Null_output_stream&
operator<<(Null_output_stream& nos, const T&)
{
  return nos;
}


} //namespace CGAL

#endif // CGAL_NULL_OUTPUT_STREAM_H
