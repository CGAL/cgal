// Copyright (c) 2017 GeometryFactory
//
// This file is part of CGAL (www.cgal.org);
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Simon Giraudot

#ifndef CGAL_IO_PLY_PLY_READER_H
#define CGAL_IO_PLY_PLY_READER_H

#include <CGAL/Container_helper.h>
#include <CGAL/IO/io.h>
#include <CGAL/type_traits/is_iterator.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/property_map.h>

#include <cstdint>
#include <boost/range/value_type.hpp>

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#define TRY_TO_GENERATE_PROPERTY(STD_TYPE, T_TYPE, TYPE)                \
  if(type == STD_TYPE  || type == T_TYPE)                              \
    m_elements.back().add_property(new PLY_read_typed_number< TYPE >(name, format))

#define TRY_TO_GENERATE_SIZED_LIST_PROPERTY(STD_SIZE_TYPE, T_SIZE_TYPE, SIZE_TYPE, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE) \
  if((size_type == STD_SIZE_TYPE  || size_type == T_SIZE_TYPE) &&      \
  (index_type == STD_INDEX_TYPE || index_type == T_INDEX_TYPE))     \
    m_elements.back().add_property(new PLY_read_typed_list_with_typed_size< SIZE_TYPE , INDEX_TYPE >(name, format))

#define TRY_TO_GENERATE_LIST_PROPERTY(STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE) \
  TRY_TO_GENERATE_SIZED_LIST_PROPERTY("uchar", "uint8", std::uint8_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE); \
  else TRY_TO_GENERATE_SIZED_LIST_PROPERTY("char", "int8", std::int8_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE); \
  else TRY_TO_GENERATE_SIZED_LIST_PROPERTY("ushort", "uint16", std::uint16_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE); \
  else TRY_TO_GENERATE_SIZED_LIST_PROPERTY("short", "int16", std::int16_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE); \
  else TRY_TO_GENERATE_SIZED_LIST_PROPERTY("uint", "uint32", std::uint32_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE); \
  else TRY_TO_GENERATE_SIZED_LIST_PROPERTY("int", "int32", std::int32_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE)

namespace CGAL {

namespace IO {

/// \cond SKIP_IN_MANUAL

// PLY types:
// name        type        number of bytes
// ---------------------------------------
// char       character                 1
// uchar      unsigned character        1
// short      short integer             2
// ushort     unsigned short integer    2
// int        integer                   4
// uint       unsigned integer          4
// float      single-precision float    4
// double     double-precision float    8

template <typename T>
struct PLY_property
{
  typedef T type;
  const char* name;
  PLY_property(const char* name) : name(name) { }
};

// Use a double property for all kernels...
template <typename FT> struct Convert_FT        { typedef double type; };
// ...except if kernel uses type float
template <>            struct Convert_FT<float> { typedef float type;  };

template <typename PointOrVectorMap>
struct Get_FT_from_map
{
  typedef typename Convert_FT<typename Kernel_traits<
                                typename boost::property_traits<
                                  PointOrVectorMap>::value_type>::Kernel::FT>::type type;
};

template <typename PointMap>
std::tuple<PointMap,
           typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3,
           PLY_property<typename Get_FT_from_map<PointMap>::type>,
           PLY_property<typename Get_FT_from_map<PointMap>::type>,
           PLY_property<typename Get_FT_from_map<PointMap>::type> >
make_ply_point_reader(PointMap point_map)
{
  return std::make_tuple(point_map, typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3(),
                         PLY_property<typename Get_FT_from_map<PointMap>::type>("x"),
                         PLY_property<typename Get_FT_from_map<PointMap>::type>("y"),
                         PLY_property<typename Get_FT_from_map<PointMap>::type>("z"));
}

template <typename VectorMap>
std::tuple<VectorMap,
           typename Kernel_traits<typename VectorMap::value_type>::Kernel::Construct_vector_3,
           PLY_property<typename Get_FT_from_map<VectorMap>::type>,
           PLY_property<typename Get_FT_from_map<VectorMap>::type>,
           PLY_property<typename Get_FT_from_map<VectorMap>::type> >
make_ply_normal_reader(VectorMap normal_map)
{
  return std::make_tuple(normal_map, typename Kernel_traits<typename VectorMap::value_type>::Kernel::Construct_vector_3(),
                         PLY_property<typename Get_FT_from_map<VectorMap>::type>("nx"),
                         PLY_property<typename Get_FT_from_map<VectorMap>::type>("ny"),
                         PLY_property<typename Get_FT_from_map<VectorMap>::type>("nz"));
}

template <typename PointMap>
std::tuple<PointMap,
           PLY_property<typename Get_FT_from_map<PointMap>::type>,
           PLY_property<typename Get_FT_from_map<PointMap>::type>,
           PLY_property<typename Get_FT_from_map<PointMap>::type> >
make_ply_point_writer(PointMap point_map)
{
  return std::make_tuple(point_map,
                         PLY_property<typename Get_FT_from_map<PointMap>::type>("x"),
                         PLY_property<typename Get_FT_from_map<PointMap>::type>("y"),
                         PLY_property<typename Get_FT_from_map<PointMap>::type>("z"));
}

template <typename VectorMap>
std::tuple<VectorMap,
           PLY_property<typename Get_FT_from_map<VectorMap>::type>,
           PLY_property<typename Get_FT_from_map<VectorMap>::type>,
           PLY_property<typename Get_FT_from_map<VectorMap>::type> >
make_ply_normal_writer(VectorMap normal_map)
{
  return std::make_tuple(normal_map,
                         PLY_property<typename Get_FT_from_map<VectorMap>::type>("nx"),
                         PLY_property<typename Get_FT_from_map<VectorMap>::type>("ny"),
                         PLY_property<typename Get_FT_from_map<VectorMap>::type>("nz"));
}

namespace internal {

class PLY_read_number
{
protected:
  std::string m_name;
  std::size_t m_format;

public:
  PLY_read_number(std::string name, std::size_t format)
    : m_name(name), m_format(format) { }
  virtual ~PLY_read_number() { }

  const std::string& name() const { return m_name; }

  virtual void get(std::istream& stream) const = 0;

  // The two following functions prevent the stream to only extract
  // ONE character (= what the types char imply) by requiring
  // explicitly an integer object when reading the stream
  void read_ascii(std::istream& stream, char& c) const
  {
    short s;
    if(stream >> s)
      c = static_cast<char>(s);
    else
    {
      c = 0;
      stream.clear(std::ios::badbit);
    }
  }

  void read_ascii(std::istream& stream, signed char& c) const
  {
    short s;
    if(stream >> s)
      c = static_cast<signed char>(s);
    else
    {
      c = 0;
      stream.clear(std::ios::badbit);
    }
  }

  void read_ascii(std::istream& stream, unsigned char& c) const
  {
    unsigned short s;
    if(stream >> s)
      c = static_cast<unsigned char>(s);
    else
    {
      c = 0;
      stream.clear(std::ios::badbit);
    }
  }

  void read_ascii(std::istream& stream, float& t) const
  {
    if(!(stream >> IO::iformat(t)))
      stream.clear(std::ios::badbit);
  }

  void read_ascii(std::istream& stream, double& t) const
  {
    if(!(stream >> IO::iformat(t)))
      stream.clear(std::ios::badbit);
  }

  // Default template when Type is not a char type
  template <typename Type>
  void read_ascii(std::istream& stream, Type& t) const
  {
    if(!(stream >> t))
      stream.clear(std::ios::badbit);
  }

  template <typename Type>
  Type read(std::istream& stream) const
  {
    if(m_format == 0) // ASCII
    {
      Type t;
      read_ascii(stream, t);
      return t;
    }
    else // Binary (2 = little endian)
    {
      union
      {
        char uChar[sizeof(Type)];
        Type type;
      } buffer;

      std::size_t size = sizeof(Type);

      stream.read(buffer.uChar, size);

      if(m_format == 2) // Big endian
      {
        for(std::size_t i = 0; i < size / 2; ++ i)
        {
          unsigned char tmp = buffer.uChar[i];
          buffer.uChar[i] = buffer.uChar[size - 1 - i];
          buffer.uChar[size - 1 - i] = tmp;
        }
      }
      return buffer.type;
    }
    return Type();
  }
};

template <typename Type>
class PLY_read_typed_number : public PLY_read_number
{
  mutable Type m_buffer;
public:
  PLY_read_typed_number(std::string name, std::size_t format)
    : PLY_read_number(name, format)
  { }

  void get(std::istream& stream) const { m_buffer =(this->read<Type>(stream)); }

  const Type& buffer() const { return m_buffer; }
};

template <typename Type>
class PLY_read_typed_list
  : public PLY_read_number
{
protected:
  mutable std::vector<Type> m_buffer;

public:
  PLY_read_typed_list(std::string name, std::size_t format)
    : PLY_read_number(name, format)
  { }

  virtual void get(std::istream& stream) const = 0;

  const std::vector<Type>& buffer() const { return m_buffer; }
};

template <typename SizeType, typename IndexType>
class PLY_read_typed_list_with_typed_size
  : public PLY_read_typed_list<IndexType>
{
public:
  PLY_read_typed_list_with_typed_size(std::string name, std::size_t format)
    : PLY_read_typed_list<IndexType>(name, format)
  { }

  void get(std::istream& stream) const
  {
    std::size_t size = static_cast<std::size_t>(this->template read<SizeType>(stream));
    this->m_buffer.resize(size);
    for(std::size_t i = 0; i < size; ++ i)
      this->m_buffer[i] = this->template read<IndexType>(stream);
  }
};

class PLY_element
{
  std::string m_name;
  std::size_t m_number;
  std::vector<PLY_read_number*> m_properties;

public:
  PLY_element(const std::string& name, std::size_t number)
    : m_name(name), m_number(number)
  { }

  PLY_element(const PLY_element& other)
    : m_name(other.m_name), m_number(other.m_number), m_properties(other.m_properties)
  {
    const_cast<PLY_element&>(other).m_properties.clear();
  }

  PLY_element& operator=(const PLY_element& other)
  {
    m_name = other.m_name;
    m_number = other.m_number;
    m_properties = other.m_properties;
    const_cast<PLY_element&>(other).m_properties.clear();
    return *this;
  }

  ~PLY_element()
  {
    for(std::size_t i = 0; i < m_properties.size(); ++ i)
      delete m_properties[i];
  }

  const std::string& name() const { return m_name; }
  std::size_t number_of_items() const { return m_number; }
  std::size_t number_of_properties() const { return m_properties.size(); }

  PLY_read_number* property(std::size_t idx) { return m_properties[idx]; }

  void add_property(PLY_read_number* read_number)
  {
    m_properties.push_back(read_number);
  }

  template <typename Type>
  bool has_property(const char* tag)
  {
    return has_property(tag, Type());
  }

  template <typename Type>
  bool has_property(const char* tag, const std::vector<Type>&)
  {
    for(std::size_t i = 0; i < number_of_properties(); ++ i)
      if(m_properties[i]->name() == tag)
        return (dynamic_cast<PLY_read_typed_list<Type>*>(m_properties[i]) != nullptr);
    return false;
  }

  template <typename Type>
  bool has_property(const char* tag, Type)
  {
    for(std::size_t i = 0; i < number_of_properties(); ++ i)
      if(m_properties[i]->name() == tag)
        return (dynamic_cast<PLY_read_typed_number<Type>*>(m_properties[i]) != nullptr);
    return false;
  }

  bool has_property(const char* tag, double)
  {
    for(std::size_t i = 0; i < number_of_properties(); ++ i)
      if(m_properties[i]->name() == tag)
        return (dynamic_cast<PLY_read_typed_number<double>*>(m_properties[i]) != nullptr
                                                                                 || dynamic_cast<PLY_read_typed_number<float>*>(m_properties[i]) != nullptr);

    return false;
  }

  template <typename Type>
  void assign(Type& t, const char* tag)
  {
    for(std::size_t i = 0; i < number_of_properties(); ++ i)
      if(m_properties[i]->name() == tag)
      {
        PLY_read_typed_number<Type>*
            property = dynamic_cast<PLY_read_typed_number<Type>*>(m_properties[i]);
        CGAL_assertion(property != nullptr);
        t = property->buffer();
        return;
      }
    t = {};
  }

  template <typename Type>
  void assign(std::vector<Type>& t, const char* tag)
  {
    for(std::size_t i = 0; i < number_of_properties(); ++ i)
      if(m_properties[i]->name() == tag)
      {
        PLY_read_typed_list<Type>*
            property = dynamic_cast<PLY_read_typed_list<Type>*>(m_properties[i]);
        CGAL_assertion(property != nullptr);
        t = property->buffer();
        return;
      }
    t = {};
  }

  void assign(double& t, const char* tag)
  {
    for(std::size_t i = 0; i < number_of_properties(); ++ i)
      if(m_properties[i]->name() == tag)
      {
        PLY_read_typed_number<double>*
            property_double = dynamic_cast<PLY_read_typed_number<double>*>(m_properties[i]);
        if(property_double == nullptr)
        {
          PLY_read_typed_number<float>*
              property_float = dynamic_cast<PLY_read_typed_number<float>*>(m_properties[i]);
          CGAL_assertion(property_float != nullptr);
          t = property_float->buffer();
        }
        else
          t = property_double->buffer();

        return;
      }
    t = {};
  }
};

class PLY_reader
{
  std::vector<PLY_element> m_elements;
  std::string m_comments;
  bool m_verbose;

public:
  PLY_reader(bool verbose) : m_verbose(verbose) { }

  std::size_t number_of_elements() const { return m_elements.size(); }
  PLY_element& element(std::size_t idx)
  {
    return m_elements[idx];
  }

  const std::string& comments() const { return m_comments; }

  template <typename Stream>
  bool init(Stream& stream)
  {
    std::size_t lineNumber = 0; // current line number
    enum Format { ASCII = 0, BINARY_LITTLE_ENDIAN = 1, BINARY_BIG_ENDIAN = 2};
    Format format = ASCII;

    std::string line;
    std::istringstream iss;

    while(getline(stream,line))
    {
      iss.clear();
      iss.str(line);
      ++ lineNumber;

      // Reads file signature on first line
      if(lineNumber == 1)
      {
        std::string signature;
        if(!(iss >> signature) || (signature != "ply"))
        {
          // if wrong file format
          if(m_verbose)
            std::cerr << "Error: incorrect file format line " << lineNumber << " of file" << std::endl;
          return false;
        }
      }

      // Reads format on 2nd line
      else if(lineNumber == 2)
      {
        std::string tag, format_string, version;
        if( !(iss >> tag >> format_string >> version) )
        {
          if(m_verbose)
            std::cerr << "Error line " << lineNumber << " of file" << std::endl;
          return false;
        }
        if(format_string == "ascii") format = ASCII;
        else if(format_string == "binary_little_endian") format = BINARY_LITTLE_ENDIAN;
        else if(format_string == "binary_big_endian") format = BINARY_BIG_ENDIAN;
        else
        {
          if(m_verbose)
            std::cerr << "Error: unknown file format \"" << format_string << "\" line " << lineNumber << std::endl;
          return false;
        }
      }

      // Comments and vertex properties
      else
      {
        std::string keyword;
        if(!(iss >> keyword))
        {
          if(m_verbose)
            std::cerr << "Error line " << lineNumber << " of file" << std::endl;
          return false;
        }

        if(keyword == "property")
        {
          std::string type, name;
          if(!(iss >> type >> name))
          {
            if(m_verbose)
              std::cerr << "Error line " << lineNumber << " of file" << std::endl;
            return false;
          }


          if(type == "list") // Special case
          {
            std::string size_type = name;
            std::string index_type;
            name.clear();
            if(!(iss >> index_type >> name))
            {
              if(m_verbose)
                std::cerr << "Error line " << lineNumber << " of file" << std::endl;
              return false;
            }

            TRY_TO_GENERATE_LIST_PROPERTY("char", "int8", std::int8_t);
            else TRY_TO_GENERATE_LIST_PROPERTY("uchar", "uint8", std::uint8_t);
            else TRY_TO_GENERATE_LIST_PROPERTY("short", "int16", std::int16_t);
            else TRY_TO_GENERATE_LIST_PROPERTY("ushort", "uint16", std::uint16_t);
            else TRY_TO_GENERATE_LIST_PROPERTY("int", "int32", std::int32_t);
            else TRY_TO_GENERATE_LIST_PROPERTY("uint", "uint32", std::uint32_t);
            else TRY_TO_GENERATE_LIST_PROPERTY("float", "float32", float);
            else TRY_TO_GENERATE_LIST_PROPERTY("double", "float64", double);
          }
          else
          {
            TRY_TO_GENERATE_PROPERTY("char", "int8", std::int8_t);
            else TRY_TO_GENERATE_PROPERTY("uchar", "uint8", std::uint8_t);
            else TRY_TO_GENERATE_PROPERTY("short", "int16", std::int16_t);
            else TRY_TO_GENERATE_PROPERTY("ushort", "uint16", std::uint16_t);
            else TRY_TO_GENERATE_PROPERTY("int", "int32", std::int32_t);
            else TRY_TO_GENERATE_PROPERTY("uint", "uint32", std::uint32_t);
            else TRY_TO_GENERATE_PROPERTY("float", "float32", float);
            else TRY_TO_GENERATE_PROPERTY("double", "float64", double);
          }

          continue;
        }
        else if(keyword == "comment")
        {
          std::string str = iss.str();
          if(str.size() > 8)
          {
            std::copy(str.begin() + 8, str.end(), std::back_inserter(m_comments));
            m_comments += "\n";
          }
        }
        else if(keyword == "element")
        {
          std::string type;
          std::size_t number;
          if(!(iss >> type >> number))
          {
            if(m_verbose)
              std::cerr << "Error line " << lineNumber << " of file" << std::endl;
            return false;
          }

          m_elements.push_back(PLY_element(type, number));
        }
        // When end_header is reached, stop loop and begin reading points
        else if(keyword == "end_header")
          break;
      }
    }
    return true;
  }

  ~PLY_reader()
  {
  }

};

template <class Reader, class T>
void get_value(Reader& r, T& v, PLY_property<T>& wrapper)
{
  return r.assign(v, wrapper.name);
}

template <std::size_t N>
struct Filler
{
  template <class Reader, class Value_tuple, class PLY_property_tuple>
  static void fill(Reader& r, Value_tuple& values, PLY_property_tuple wrappers)
  {
    get_value(r, std::get<N>(values), std::get<N+2>(wrappers));
    Filler<N-1>::fill(r, values, wrappers);
  }
};

template<int ...>
struct seq { };

template<int N, int ...S>
struct gens : gens<N-1, N-1, S...> { };

template<int ...S>
struct gens<0, S...> {
  typedef seq<S...> type;
};

template<class ValueType, class Functor, class Tuple, int ...S>
ValueType call_functor(Functor f, Tuple t, seq<S...>) {
  return f(std::get<S>(t) ...);
}

template <class ValueType, class Functor, typename ... T>
ValueType call_functor(Functor f, std::tuple<T...>& t)
{
  return call_functor<ValueType>(f, t, typename gens<sizeof...(T)>::type());
}

template<>
struct Filler<0>
{
  template <class Reader, class Value_tuple, class PLY_property_tuple>
  static void fill(Reader& r, Value_tuple& values, PLY_property_tuple wrappers)
  {
    get_value(r, std::get<0>(values), std::get<2>(wrappers));
  }
};

template <typename OutputValueType,
          typename PropertyMap,
          typename Constructor,
          typename ... T>
void process_properties(PLY_element& element, OutputValueType& new_element,
                        std::tuple<PropertyMap, Constructor, PLY_property<T>...>&& current)
{
  typedef typename boost::property_traits<PropertyMap>::value_type PmapValueType;

  std::tuple<T...> values;
  Filler<sizeof...(T)-1>::fill(element, values, current);
  PmapValueType new_value = call_functor<PmapValueType>(std::get<1>(current), values);
  put(std::get<0>(current), new_element, new_value);
}

template <typename OutputValueType,
          typename PropertyMap,
          typename Constructor,
          typename ... T,
          typename NextPropertyBinder,
          typename ... PropertyMapBinders>
void process_properties(PLY_element& element, OutputValueType& new_element,
                        std::tuple<PropertyMap, Constructor, PLY_property<T>...>&& current,
                        NextPropertyBinder&& next,
                        PropertyMapBinders&& ... properties)
{
  typedef typename boost::property_traits<PropertyMap>::value_type PmapValueType;

  std::tuple<T...> values;
  Filler<sizeof...(T)-1>::fill(element, values, current);
  PmapValueType new_value = call_functor<PmapValueType>(std::get<1>(current), values);
  put(std::get<0>(current), new_element, new_value);

  process_properties(element, new_element, std::forward<NextPropertyBinder>(next),
                     std::forward<PropertyMapBinders>(properties)...);
}


template <typename OutputValueType, typename PropertyMap, typename T>
void process_properties(PLY_element& element, OutputValueType& new_element,
                        std::pair<PropertyMap, PLY_property<T> >&& current)
{
  T new_value = T();
  element.assign(new_value, current.second.name);
  put(current.first, new_element, new_value);
}

template <typename OutputValueType, typename PropertyMap, typename T,
          typename NextPropertyBinder, typename ... PropertyMapBinders>
void process_properties(PLY_element& element, OutputValueType& new_element,
                        std::pair<PropertyMap, PLY_property<T> >&& current,
                        NextPropertyBinder&& next,
                        PropertyMapBinders&& ... properties)
{
  T new_value = T();
  element.assign(new_value, current.second.name);
  put(current.first, new_element, new_value);
  process_properties(element, new_element, std::forward<NextPropertyBinder>(next),
                     std::forward<PropertyMapBinders>(properties)...);
}

template <typename Integer, class PolygonRange, class ColorOutputIterator>
bool read_PLY_faces(std::istream& in,
                    PLY_element& element,
                    PolygonRange& polygons,
                    ColorOutputIterator fc_out,
                    const char* vertex_indices_tag,
                    std::enable_if_t<CGAL::is_iterator<ColorOutputIterator>::value>* = nullptr)
{
  typedef CGAL::IO::Color                                 Color_rgb;

  bool has_colors = false;
  std::string rtag = "r", gtag = "g", btag = "b";

  if((element.has_property<std::uint8_t>("red") || element.has_property<std::uint8_t>("r")) &&
     (element.has_property<std::uint8_t>("green") || element.has_property<std::uint8_t>("g")) &&
     (element.has_property<std::uint8_t>("blue") || element.has_property<std::uint8_t>("b")))
  {
    has_colors = true;
    if(element.has_property<std::uint8_t>("red"))
    {
      rtag = "red";
      gtag = "green";
      btag = "blue";
    }
  }

  for(std::size_t j = 0; j < element.number_of_items(); ++ j)
  {
    for(std::size_t k = 0; k < element.number_of_properties(); ++ k)
    {
      PLY_read_number* property = element.property(k);
      property->get(in);

      if(in.fail())
        return false;
    }

    std::tuple<std::vector<Integer>, std::uint8_t, std::uint8_t, std::uint8_t> new_face;

    if(has_colors)
    {
      process_properties(element, new_face,
                         std::make_pair(CGAL::make_nth_of_tuple_property_map<0>(new_face),
                                        PLY_property<std::vector<Integer> >(vertex_indices_tag)),
                         std::make_pair(CGAL::make_nth_of_tuple_property_map<1>(new_face),
                                        PLY_property<std::uint8_t>(rtag.c_str())),
                         std::make_pair(CGAL::make_nth_of_tuple_property_map<2>(new_face),
                                        PLY_property<std::uint8_t>(gtag.c_str())),
                         std::make_pair(CGAL::make_nth_of_tuple_property_map<3>(new_face),
                                        PLY_property<std::uint8_t>(btag.c_str())));

      *fc_out++ = Color_rgb(get<1>(new_face), get<2>(new_face), get<3>(new_face));
    }
    else
    {
      process_properties(element, new_face,
                         std::make_pair(CGAL::make_nth_of_tuple_property_map<0>(new_face),
                                        PLY_property<std::vector<Integer> >(vertex_indices_tag)));
    }

    polygons.emplace_back();
    ::CGAL::internal::resize(polygons.back(), get<0>(new_face).size());
    for(std::size_t i = 0; i < get<0>(new_face).size(); ++ i)
      polygons.back()[i] = std::size_t(get<0>(new_face)[i]);
  }

  return true;
}

template <typename Integer, class PolygonRange, class ColorRange>
bool read_PLY_faces(std::istream& in,
                    PLY_element& element,
                    PolygonRange& polygons,
                    ColorRange& fcolors,
                    const char* vertex_indices_tag,
                    std::enable_if_t<
                      boost::has_range_const_iterator<ColorRange>::value
                    >* = nullptr)
{
  return read_PLY_faces<Integer>(in, element, polygons, std::back_inserter(fcolors), vertex_indices_tag);
}

} // namespace PLY
} // namespace internal

#ifndef CGAL_NO_DEPREACTED_CODE
using IO::PLY_property;
using IO::make_ply_normal_reader;
using IO::make_ply_normal_writer;
using IO::make_ply_point_reader;
using IO::make_ply_point_writer;
#endif

/// \endcond

} // namespace CGAL

#endif // CGAL_IO_PLY_PLY_READER_H
