// Copyright (c) 2017  GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_IO_PLY_READ_PLY_POINTS_H
#define CGAL_IO_PLY_READ_PLY_POINTS_H


#include <CGAL/config.h>

#include <CGAL/IO/PLY.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/IO/io.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <boost/version.hpp>
#include <boost/cstdint.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>

namespace CGAL {

namespace IO {

// documentation in ../PLY.h
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename ... PropertyHandler>
bool read_PLY_with_properties(std::istream& is,
                              PointOutputIterator output,
                              PropertyHandler&& ... properties)
{
  if(!is)
    return false;

  internal::PLY_reader reader(true);

  if(!(reader.init(is)))
  {
    is.setstate(std::ios::failbit);
    return false;
  }

  for(std::size_t i = 0; i < reader.number_of_elements(); ++ i)
  {
    internal::PLY_element& element = reader.element(i);

    for(std::size_t j = 0; j < element.number_of_items(); ++ j)
    {
      for(std::size_t k = 0; k < element.number_of_properties(); ++ k)
      {
        internal::PLY_read_number* property = element.property(k);
        property->get(is);

        if(is.fail())
          return false;
      }

      if(element.name() == "vertex" || element.name() == "vertices")
      {
        OutputIteratorValueType new_element;
        internal::process_properties(element, new_element, std::forward<PropertyHandler>(properties)...);
        *(output ++) = new_element;
      }
    }
  }

  return true;
}

/// \cond SKIP_IN_MANUAL

template <typename OutputIterator,
          typename ... PropertyHandler>
bool read_PLY_with_properties(std::istream& is,
                              OutputIterator output,
                              PropertyHandler&& ... properties)
{
  typedef typename value_type_traits<OutputIterator>::type OutputValueType;

  return read_PLY_with_properties<OutputValueType>(is, output, std::forward<PropertyHandler>(properties)...);
}

/// \endcond


// documented in ../PLY.h
template <typename OutputIteratorValueType,
typename PointOutputIterator,
 typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool read_PLY(std::istream& is,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np ,
              std::enable_if_t<CGAL::is_iterator<PointOutputIterator>::value>*
              )
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;

  // basic geometric types
  typedef Point_set_processing_3_np_helper<PointRange, CGAL_NP_CLASS> NP_helper;
  typedef typename NP_helper::Point_map PointMap;
  typedef typename NP_helper::Normal_map NormalMap;

  PointMap point_map = NP_helper::get_point_map(np);
  NormalMap normal_map = NP_helper::get_normal_map(np);

  return read_PLY_with_properties(is, output,
                                  make_ply_point_reader(point_map),
                                  make_ply_normal_reader(normal_map));
}

// documentation in ../PLY.h
template <typename OutputIteratorValueType,
          typename PointOutputIterator,
          typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT>
bool read_PLY(const std::string& fname,
              PointOutputIterator output,
              const CGAL_NP_CLASS& np,
              std::enable_if_t<CGAL::is_iterator<PointOutputIterator>::value>*
              )
{
  const bool binary = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, internal_np::use_binary_mode), true);
  if(binary)
  {
    std::ifstream is(fname, std::ios::binary);
    CGAL::IO::set_mode(is, CGAL::IO::BINARY);
    return read_PLY<OutputIteratorValueType>(is, output, np);
  }
  else
  {
    std::ifstream is(fname);
    CGAL::IO::set_mode(is, CGAL::IO::ASCII);
    return read_PLY<OutputIteratorValueType>(is, output, np);
  }
}

/// \cond SKIP_IN_MANUAL

// variants with default output iterator value type
template <typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(std::istream& is, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values(),
              std::enable_if_t<CGAL::is_iterator<OutputIterator>::value>* = nullptr)
{
  return read_PLY<typename value_type_traits<OutputIterator>::type>(is, output, np);
}

template <typename OutputIterator,typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_PLY(const std::string& fname, OutputIterator output, const CGAL_NP_CLASS& np = parameters::default_values(),
              std::enable_if_t<CGAL::is_iterator<OutputIterator>::value>* = nullptr)
{
  return read_PLY<typename value_type_traits<OutputIterator>::type>(fname, output, np);
}

/// \endcond

} // namespace IO

} // namespace CGAL

#undef TRY_TO_GENERATE_POINT_PROPERTY
#undef TRY_TO_GENERATE_SIZED_FACE_PROPERTY
#undef TRY_TO_GENERATE_FACE_PROPERTY

#endif // CGAL_IO_PLY_READ_PLY_POINTS_H
