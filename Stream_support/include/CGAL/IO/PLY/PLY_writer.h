// Copyright (c) 2017 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_IO_PLY_PLY_WRITER_H
#define CGAL_IO_PLY_PLY_WRITER_H

#include <CGAL/IO/io.h>

namespace CGAL {
namespace IO {
namespace internal {

template <typename T> inline void property_header_type (std::ostream& stream)
{
  CGAL_assertion_msg (false, "Unknown PLY type");
  stream << "undefined_type";
}

template <> inline void property_header_type<char> (std::ostream& stream) { stream << "char"; }
template <> inline void property_header_type<signed char> (std::ostream& stream) { stream << "char"; }
template <> inline void property_header_type<unsigned char> (std::ostream& stream) { stream << "uchar"; }
template <> inline void property_header_type<short> (std::ostream& stream) { stream << "short"; }
template <> inline void property_header_type<unsigned short> (std::ostream& stream) { stream << "ushort"; }
template <> inline void property_header_type<int> (std::ostream& stream) { stream << "int"; }
template <> inline void property_header_type<unsigned int> (std::ostream& stream) { stream << "uint"; }
template <> inline void property_header_type<float> (std::ostream& stream) { stream << "float"; }
template <> inline void property_header_type<double> (std::ostream& stream) { stream << "double"; }

template <typename T>
void property_header (std::ostream& stream, const PLY_property<T>& prop)
{
  stream << "property ";
  property_header_type<T>(stream);
  stream << " " << prop.name << std::endl;
}

template <typename T>
void property_header (std::ostream& stream, const PLY_property<std::vector<T> >& prop)
{
  stream << "property list uchar ";
  property_header_type<T>(stream);
  stream << " " << prop.name << std::endl;
}

template <std::size_t N>
struct Properties_header
{
  template <class PLY_property_tuple>
  static void write(std::ostream& stream, PLY_property_tuple& wrappers)
  {
    Properties_header<N-1>::write(stream, wrappers);
    property_header (stream, std::get<N+1>(wrappers));
  }
};
template <>
struct Properties_header<0>
{
  template <class PLY_property_tuple>
  static void write(std::ostream& stream, PLY_property_tuple& wrappers)
  {
    property_header (stream, std::get<1>(wrappers));
  }
};

template <typename PropertyMap,
          typename ... T>
void output_property_header(std::ostream& stream,
                            std::tuple<PropertyMap, PLY_property<T>... >&& current)
{
  Properties_header<sizeof...(T)-1>::write(stream, current);
}

template <typename PropertyMap,
          typename T>
void output_property_header(std::ostream& stream,
                            std::pair<PropertyMap, PLY_property<T> >&& current)
{
  property_header(stream, current.second);
}

template <typename PropertyMap,
          typename T,
          typename NextPropertyHandler,
          typename ... PropertyHandler>
void output_property_header(std::ostream& stream,
                            std::pair<PropertyMap, PLY_property<T> >&& current,
                            NextPropertyHandler&& next,
                            PropertyHandler&& ... properties)
{
  property_header(stream, current.second);
  output_property_header(stream, std::forward<NextPropertyHandler>(next),
                         std::forward<PropertyHandler>(properties)...);
}
template <typename PropertyMap,
          typename ... T,
          typename NextPropertyHandler,
          typename ... PropertyHandler>
void output_property_header(std::ostream& stream,
                            std::tuple<PropertyMap, PLY_property<T>... >&& current,
                            NextPropertyHandler&& next,
                            PropertyHandler&& ... properties)
{
  Properties_header<sizeof...(T)-1>::write(stream, current);
  output_property_header(stream, std::forward<NextPropertyHandler>(next),
                         std::forward<PropertyHandler>(properties)...);
}

template <typename ForwardIterator,
          typename PropertyMap>
void property_write(std::ostream& stream, ForwardIterator it, PropertyMap map)
{
  stream << CGAL::IO::oformat(get(map, *it));
}

template <typename T>
inline T no_char_character(const T& t) { return t; }
inline int no_char_character(const char& t) { return int(t); }
inline int no_char_character(const signed char& t) { return int(t); }
inline int no_char_character(const unsigned char& t) { return int(t); }

template <typename ForwardIterator,
          typename PropertyMap,
          typename T>
void simple_property_write(std::ostream& stream,
                           ForwardIterator it,
                           std::pair<PropertyMap, PLY_property<T> > map)
{
  if(CGAL::IO::get_mode(stream) == CGAL::IO::ASCII)
    stream << no_char_character(get(map.first, *it));
  else
  {
    typename boost::property_traits<PropertyMap>::value_type value = get(map.first, *it);
    stream.write(reinterpret_cast<char*>(&value), sizeof(value));
  }
}

template <typename ForwardIterator,
          typename PropertyMap,
          typename T>
void simple_property_write(std::ostream& stream,
                           ForwardIterator it,
                           std::pair<PropertyMap, PLY_property<std::vector<T> > > map)
{
  const typename PropertyMap::reference value = get(map.first, *it);

  if(CGAL::IO::get_mode(stream) == CGAL::IO::ASCII)
  {
    stream << value.size();
    for(std::size_t i = 0; i < value.size(); ++ i)
      stream << " " << no_char_character(value[i]);
  }
  else
  {
    unsigned char size = static_cast<unsigned char>(value.size());
    stream.write(reinterpret_cast<char*>(&size), sizeof(size));
    for(std::size_t i = 0; i < value.size(); ++ i)
    {
      T t = T(value[i]);
      stream.write(reinterpret_cast<char*>(&t), sizeof(t));
    }
  }
}

template <typename ForwardIterator,
          typename PropertyMap,
          typename ... T>
void output_properties(std::ostream& stream,
                       ForwardIterator it,
                       std::tuple<PropertyMap, PLY_property<T>... >&& current)
{
  property_write(stream, it, std::get<0>(current));
  if(CGAL::IO::get_mode(stream) == CGAL::IO::ASCII)
    stream << std::endl;
}

template <typename ForwardIterator,
          typename PropertyMap,
          typename T>
void output_properties(std::ostream& stream,
                       ForwardIterator it,
                       std::pair<PropertyMap, PLY_property<T> >&& current)
{
  simple_property_write(stream, it, std::forward<std::pair<PropertyMap, PLY_property<T> > >(current));
  if(get_mode(stream) == CGAL::IO::ASCII)
    stream << std::endl;
}

template <typename ForwardIterator,
          typename PropertyMap,
          typename T,
          typename NextPropertyHandler,
          typename ... PropertyHandler>
void output_properties(std::ostream& stream,
                       ForwardIterator it,
                       std::pair<PropertyMap, PLY_property<T> >&& current,
                       NextPropertyHandler&& next,
                       PropertyHandler&& ... properties)
{
  simple_property_write(stream, it, current);
  if(get_mode(stream) == CGAL::IO::ASCII)
    stream << " ";
  output_properties(stream, it, std::forward<NextPropertyHandler>(next),
                    std::forward<PropertyHandler>(properties)...);
}

template <typename ForwardIterator,
          typename PropertyMap,
          typename ... T,
          typename NextPropertyHandler,
          typename ... PropertyHandler>
void output_properties(std::ostream& stream,
                       ForwardIterator it,
                       std::tuple<PropertyMap, PLY_property<T>... >&& current,
                       NextPropertyHandler&& next,
                       PropertyHandler&& ... properties)
{
  property_write(stream, it, std::get<0>(current));
  if(get_mode(stream) == CGAL::IO::ASCII)
    stream << " ";
  output_properties(stream, it, std::forward<NextPropertyHandler>(next),
                    std::forward<PropertyHandler>(properties)...);
}

// Printer classes used by Point_set_3 and Surface_mesh(translate a
// property map to a PLY property)

template <typename Index>
class Abstract_property_printer
{
public:
  virtual ~Abstract_property_printer() { }
  virtual void print(std::ostream& stream, const Index& index) = 0;
};

template <typename Index, typename PropertyMap>
class Property_printer
  : public Abstract_property_printer<Index>
{
  PropertyMap m_pmap;
public:
  Property_printer(const PropertyMap& pmap) : m_pmap(pmap) { }

  virtual void print(std::ostream& stream, const Index& index)
  {
    stream << get(m_pmap, index);
  }
};

template <typename Index,
          typename PropertyMap,
          typename Type = typename boost::property_traits<PropertyMap>::value_type>
class Simple_property_printer
  : public Abstract_property_printer<Index>
{
  PropertyMap m_pmap;
public:
  Simple_property_printer(const PropertyMap& pmap) : m_pmap(pmap) { }

  virtual void print(std::ostream& stream, const Index& index)
  {
    if(get_mode(stream) == CGAL::IO::ASCII)
      stream << get(m_pmap, index);
    else
    {
      Type t = Type(get(m_pmap, index));
      stream.write(reinterpret_cast<char*>(&t), sizeof(t));
    }
  }
};

template <typename Index, typename PropertyMap>
class Char_property_printer
  : public Abstract_property_printer<Index>
{
  typedef typename boost::property_traits<PropertyMap>::value_type Type;

  PropertyMap m_pmap;

public:
  Char_property_printer(const PropertyMap& pmap) : m_pmap(pmap) { }

  virtual void print(std::ostream& stream, const Index& index)
  {
    if(get_mode(stream) == CGAL::IO::ASCII)
      stream << int(get(m_pmap, index));
    else
    {
      Type t = get(m_pmap, index);
      stream.write(reinterpret_cast<char*>(&t), sizeof(t));
    }
  }
};

} // namespace internal
} // namespace IO
} // namespace CGAL

#endif // CGAL_IO_PLY_PLY_WRITER_H
