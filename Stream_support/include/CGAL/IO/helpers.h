// Copyright (c) 2015 GeometryFactory
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STREAM_SUPPORT_HELPERS_H
#define CGAL_STREAM_SUPPORT_HELPERS_H

#include <CGAL/array.h>
#include <CGAL/assertions.h>
#include <CGAL/Container_helper.h>
#include <CGAL/Has_member.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/type_traits/is_iterator.h>

#include <boost/mpl/logical.hpp>
#include <boost/mpl/has_xxx.hpp>
#include <boost/range/has_range_iterator.hpp>

#include <algorithm>
#include <string>

namespace CGAL {
namespace IO {
namespace internal {

template <typename Kernel>
void fill_point(const double x, const double y, const double /*z*/, const double w,
                CGAL::Point_2<Kernel>& pt)
{
  typedef typename Kernel::FT FT;
  pt = CGAL::Point_2<Kernel>(FT(x)/FT(w), FT(y)/FT(w));
}

template <typename Kernel>
void fill_point(const double x, const double y, const double z, const double w,
                CGAL::Point_3<Kernel>& pt)
{
  typedef typename Kernel::FT FT;
  pt = CGAL::Point_3<Kernel>(FT(x)/FT(w), FT(y)/FT(w), FT(z)/FT(w));
}

// something else (arrays...)
template <typename Point>
void fill_point(const double x, const double y, const double z, const double w,
                Point& pt)
{
  CGAL::internal::resize(pt, 3);
  pt[0] = x / w; pt[1] = y / w; pt[2] = z / w;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

static inline std::string get_file_extension(const std::string fname)
{
  std::string::size_type dot(fname.rfind("."));
  if(dot == std::string::npos)
    return std::string();

  std::string ext = fname.substr(dot+1, fname.length() - dot - 1);
  std::transform(ext.begin(), ext.end(), ext.begin(),
                 [](char c) {
                   return static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
                 });

  return ext;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Functions like 'write_OFF' can take :
// - write_OFF(stream, point_set)
// - write_OFF(stream, pointrange)
// - write_OFF(stream, polygon_mesh)
// so a way to distinguish is needed
BOOST_MPL_HAS_XXX_TRAIT_DEF(Point_set)

template <typename T>
struct is_Point_set_3 : has_Point_set<T> { };

// Point_set_3 and strings also functions as ranges, but we want to match polygon soups here
template <typename T>
struct is_Range
  : public std::bool_constant<
             boost::has_range_const_iterator<T>::value && // should be a range
             !is_Point_set_3<T>::value && // but not a Point_set_3
             !std::is_convertible_v<T, std::string> > // or a std::string / char [x]
{ };

template <class T>
inline constexpr bool is_Range_v = is_Range<T>::value;

// For polygon meshes
template <typename T>
struct is_Point_set_or_Range_or_Iterator
  : public std::bool_constant<is_Point_set_3<T>::value || is_Range<T>::value || is_iterator<T>::value >
{ };

template <class T>
inline constexpr bool is_Point_set_or_Range_or_Iterator_v = is_Point_set_or_Range_or_Iterator<T>::value;

} // end namespace internal
} // end namespace IO
} // namespace CGAL

#endif // CGAL_STREAM_SUPPORT_HELPERS_H
