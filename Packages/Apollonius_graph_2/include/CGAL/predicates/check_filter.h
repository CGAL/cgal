#ifndef CGAL_CHECK_FILTER_H
#define CGAL_CHECK_FILTER_H

#include <CGAL/Filtered_exact.h>

#undef CGAL_IA_NEW_FILTERS

namespace CGAL {

template < class T>
void must_be_filtered(const T&)
{}

template < class CT, class ET, class Type, bool Protection, class Cache>
void must_be_filtered(const Filtered_exact<CT, ET, Type, Protection,
		      Cache> &)
{ dont_compile(); }

}

#endif
