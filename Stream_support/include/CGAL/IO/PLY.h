// Copyright (c) 2015  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Simon Giraudot

#ifndef CGAL_IO_PLY_H
#define CGAL_IO_PLY_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/IO/io.h>
#include <CGAL/property_map.h>

#include <boost/cstdint.hpp>

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#define TRY_TO_GENERATE_PROPERTY(STD_TYPE, T_TYPE, TYPE)                \
  if (type == STD_TYPE  || type == T_TYPE)                              \
    m_elements.back().add_property (new PLY_read_typed_number< TYPE > (name, format))

#define TRY_TO_GENERATE_SIZED_LIST_PROPERTY(STD_SIZE_TYPE, T_SIZE_TYPE, SIZE_TYPE, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE) \
  if ((size_type == STD_SIZE_TYPE  || size_type == T_SIZE_TYPE) &&      \
      (index_type == STD_INDEX_TYPE || index_type == T_INDEX_TYPE))     \
    m_elements.back().add_property (new PLY_read_typed_list_with_typed_size< SIZE_TYPE , INDEX_TYPE > (name, format))

#define TRY_TO_GENERATE_LIST_PROPERTY(STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE) \
  TRY_TO_GENERATE_SIZED_LIST_PROPERTY("uchar", "uint8", boost::uint8_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE); \
  else TRY_TO_GENERATE_SIZED_LIST_PROPERTY("ushort", "uint16", boost::uint16_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE); \
  else TRY_TO_GENERATE_SIZED_LIST_PROPERTY("uint", "uint32", boost::uint32_t, STD_INDEX_TYPE, T_INDEX_TYPE, INDEX_TYPE)


/// \cond SKIP_IN_MANUAL

namespace CGAL {

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
  PLY_property (const char* name) : name (name) { }
};

// Use a double property for all kernels...
template <typename FT> struct Convert_FT        { typedef double type; };
// ...except if kernel uses type float
template <>            struct Convert_FT<float> { typedef float type;  };

template <typename PointOrVectorMap>
struct Get_FT_from_map
{
  typedef typename Convert_FT
  <typename Kernel_traits
   <typename boost::property_traits
    <PointOrVectorMap>::value_type>::Kernel::FT>::type type;
};

template <typename PointMap>
std::tuple<PointMap,
           typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3,
           PLY_property<typename Get_FT_from_map<PointMap>::type>,
           PLY_property<typename Get_FT_from_map<PointMap>::type>,
           PLY_property<typename Get_FT_from_map<PointMap>::type> >
make_ply_point_reader(PointMap point_map)
{
  return std::make_tuple (point_map, typename Kernel_traits<typename PointMap::value_type>::Kernel::Construct_point_3(),
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
  return std::make_tuple (normal_map, typename Kernel_traits<typename VectorMap::value_type>::Kernel::Construct_vector_3(),
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
  return std::make_tuple (point_map,
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
  return std::make_tuple (normal_map,
                          PLY_property<typename Get_FT_from_map<VectorMap>::type>("nx"),
                          PLY_property<typename Get_FT_from_map<VectorMap>::type>("ny"),
                          PLY_property<typename Get_FT_from_map<VectorMap>::type>("nz"));
}

namespace internal {

namespace PLY {

class PLY_read_number
{
protected:
  std::string m_name;
  std::size_t m_format;

public:
  PLY_read_number (std::string name, std::size_t format)
    : m_name (name), m_format (format) { }
  virtual ~PLY_read_number() { }

  const std::string& name () const { return m_name; }

  virtual void get (std::istream& stream) const = 0;

  // The two following functions prevent the stream to only extract
  // ONE character (= what the types char imply) by requiring
  // explicitely an integer object when reading the stream
  void read_ascii (std::istream& stream, char& c) const
  {
    short s;
    stream >> s;
    c = static_cast<char>(s);
  }
  void read_ascii (std::istream& stream, signed char& c) const
  {
    short s;
    stream >> s;
    c = static_cast<signed char>(s);
  }
  void read_ascii (std::istream& stream, unsigned char& c) const
  {
    unsigned short s;
    stream >> s;
    c = static_cast<unsigned char>(s);
  }

  void read_ascii (std::istream& stream, float& t) const
  {
    stream >> iformat(t);
  }

  void read_ascii (std::istream& stream, double& t) const
  {
    stream >> iformat(t);
  }

  // Default template when Type is not a char type
  template <typename Type>
  void read_ascii (std::istream& stream, Type& t) const
  {
    stream >> t;
  }


  template <typename Type>
  Type read (std::istream& stream) const
  {
    if (m_format == 0) // Ascii
    {
      Type t;
      read_ascii (stream, t);
      return t;
    }
    else // Binary (2 = little endian)
    {
      union
      {
        char uChar[sizeof (Type)];
        Type type;
      } buffer;

      std::size_t size = sizeof (Type);

      stream.read(buffer.uChar, size);

      if (m_format == 2) // Big endian
      {
        for (std::size_t i = 0; i < size / 2; ++ i)
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
  PLY_read_typed_number (std::string name, std::size_t format)
    : PLY_read_number (name, format)
  {
  }
  void get (std::istream& stream) const
  {
    m_buffer = (this->read<Type> (stream));
  }
  const Type& buffer() const
  {
    return m_buffer;
  }
};

template <typename Type>
class PLY_read_typed_list : public PLY_read_number
{
protected:
  mutable std::vector<Type> m_buffer;
public:
  PLY_read_typed_list (std::string name, std::size_t format)
    : PLY_read_number (name, format)
  {
  }
  virtual void get (std::istream& stream) const = 0;

  const std::vector<Type>& buffer() const
  {
    return m_buffer;
  }
};

template <typename SizeType, typename IndexType>
class PLY_read_typed_list_with_typed_size
  : public PLY_read_typed_list<IndexType>
{

public:
  PLY_read_typed_list_with_typed_size (std::string name, std::size_t format)
    : PLY_read_typed_list<IndexType> (name, format)
  {
  }
  void get (std::istream& stream) const
  {
    std::size_t size = static_cast<std::size_t>(this->template read<SizeType>(stream));
    this->m_buffer.resize (size);
    for (std::size_t i = 0; i < size; ++ i)
      this->m_buffer[i] = this->template read<IndexType> (stream);
  }
};

class PLY_element
{
  std::string m_name;
  std::size_t m_number;

  std::vector<PLY_read_number*> m_properties;
public:

  PLY_element (const std::string& name, std::size_t number)
    : m_name (name), m_number (number)
  { }

  PLY_element (const PLY_element& other)
    : m_name (other.m_name), m_number (other.m_number), m_properties (other.m_properties)
  {
    const_cast<PLY_element&>(other).m_properties.clear();
  }

  PLY_element& operator= (const PLY_element& other)
  {
    m_name = other.m_name;
    m_number = other.m_number;
    m_properties = other.m_properties;
    const_cast<PLY_element&>(other).m_properties.clear();
    return *this;
  }

  ~PLY_element()
  {
    for (std::size_t i = 0; i < m_properties.size(); ++ i)
      delete m_properties[i];
  }

  const std::string& name() const { return m_name; }
  std::size_t number_of_items() const { return m_number; }
  std::size_t number_of_properties() const { return m_properties.size(); }

  PLY_read_number* property (std::size_t idx) { return m_properties[idx]; }

  void add_property (PLY_read_number* read_number)
  {
    m_properties.push_back (read_number);
  }

  template <typename Type>
  bool has_property (const char* tag)
  {
    return has_property (tag, Type());
  }
  template <typename Type>
  bool has_property (const char* tag, const std::vector<Type>&)
  {
    for (std::size_t i = 0; i < number_of_properties(); ++ i)
      if (m_properties[i]->name () == tag)
        return (dynamic_cast<PLY_read_typed_list<Type>*>(m_properties[i]) != nullptr);
    return false;
  }

  template <typename Type>
  bool has_property (const char* tag, Type)
  {
    for (std::size_t i = 0; i < number_of_properties(); ++ i)
      if (m_properties[i]->name () == tag)
        return (dynamic_cast<PLY_read_typed_number<Type>*>(m_properties[i]) != nullptr);
    return false;
  }
  bool has_property (const char* tag, double)
  {
    for (std::size_t i = 0; i < number_of_properties(); ++ i)
      if (m_properties[i]->name () == tag)
        return (dynamic_cast<PLY_read_typed_number<double>*>(m_properties[i]) != nullptr
                || dynamic_cast<PLY_read_typed_number<float>*>(m_properties[i]) != nullptr);

    return false;
  }

  template <typename Type>
  void assign (Type& t, const char* tag)
  {
    for (std::size_t i = 0; i < number_of_properties (); ++ i)
      if (m_properties[i]->name () == tag)
      {
        PLY_read_typed_number<Type>*
          property = dynamic_cast<PLY_read_typed_number<Type>*>(m_properties[i]);
        CGAL_assertion (property != nullptr);
        t = property->buffer();
        return;
      }
    t = {};
  }

  template <typename Type>
  void assign (std::vector<Type>& t, const char* tag)
  {
    for (std::size_t i = 0; i < number_of_properties (); ++ i)
      if (m_properties[i]->name () == tag)
      {
        PLY_read_typed_list<Type>*
          property = dynamic_cast<PLY_read_typed_list<Type>*>(m_properties[i]);
        CGAL_assertion (property != nullptr);
        t = property->buffer();
        return;
      }
    t = {};
  }

  void assign (double& t, const char* tag)
  {
    for (std::size_t i = 0; i < number_of_properties (); ++ i)
      if (m_properties[i]->name () == tag)
      {
        PLY_read_typed_number<double>*
          property_double = dynamic_cast<PLY_read_typed_number<double>*>(m_properties[i]);
        if (property_double == nullptr)
        {
          PLY_read_typed_number<float>*
            property_float = dynamic_cast<PLY_read_typed_number<float>*>(m_properties[i]);
          CGAL_assertion (property_float != nullptr);
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

public:
  PLY_reader () { }

  std::size_t number_of_elements() const { return m_elements.size(); }
  PLY_element& element (std::size_t idx)
  {
    return m_elements[idx];
  }

  const std::string& comments() const { return m_comments; }

  template <typename Stream>
  bool init (Stream& stream)
  {
    std::size_t lineNumber = 0; // current line number
    enum Format { ASCII = 0, BINARY_LITTLE_ENDIAN = 1, BINARY_BIG_ENDIAN = 2};
    Format format = ASCII;

    std::string line;
    std::istringstream iss;

    while (getline (stream,line))
    {
      iss.clear();
      iss.str (line);
      ++ lineNumber;

      // Reads file signature on first line
      if (lineNumber == 1)
      {
        std::string signature;
        if (!(iss >> signature) || (signature != "ply"))
        {
          // if wrong file format
          std::cerr << "Error: incorrect file format line " << lineNumber << " of file" << std::endl;
          return false;
        }
      }

      // Reads format on 2nd line
      else if (lineNumber == 2)
      {
        std::string tag, format_string, version;
        if ( !(iss >> tag >> format_string >> version) )
        {
          std::cerr << "Error line " << lineNumber << " of file" << std::endl;
          return false;
        }
        if (format_string == "ascii") format = ASCII;
        else if (format_string == "binary_little_endian") format = BINARY_LITTLE_ENDIAN;
        else if (format_string == "binary_big_endian") format = BINARY_BIG_ENDIAN;
        else
        {
          std::cerr << "Error: unknown file format \"" << format_string << "\" line " << lineNumber << std::endl;
          return false;
        }
      }

      // Comments and vertex properties
      else
      {
        std::string keyword;
        if (!(iss >> keyword))
        {
          std::cerr << "Error line " << lineNumber << " of file" << std::endl;
          return false;
        }

        if (keyword == "property")
        {
          std::string type, name;
          if (!(iss >> type >> name))
          {
            std::cerr << "Error line " << lineNumber << " of file" << std::endl;
            return false;
          }


          if (type == "list") // Special case
          {
            std::string size_type = name;
            std::string index_type;
            name.clear();
            if (!(iss >> index_type >> name))
            {
              std::cerr << "Error line " << lineNumber << " of file" << std::endl;
              return false;
            }

            TRY_TO_GENERATE_LIST_PROPERTY ("char", "int8", boost::int8_t);
            else TRY_TO_GENERATE_LIST_PROPERTY ("uchar", "uint8", boost::uint8_t);
            else TRY_TO_GENERATE_LIST_PROPERTY ("short", "int16", boost::int16_t);
            else TRY_TO_GENERATE_LIST_PROPERTY ("ushort", "uint16", boost::uint16_t);
            else TRY_TO_GENERATE_LIST_PROPERTY ("int", "int32", boost::int32_t);
            else TRY_TO_GENERATE_LIST_PROPERTY ("uint", "uint32", boost::uint32_t);
            else TRY_TO_GENERATE_LIST_PROPERTY ("float", "float32", float);
            else TRY_TO_GENERATE_LIST_PROPERTY ("double", "float64", double);
          }
          else
          {
            TRY_TO_GENERATE_PROPERTY ("char", "int8", boost::int8_t);
            else TRY_TO_GENERATE_PROPERTY ("uchar", "uint8", boost::uint8_t);
            else TRY_TO_GENERATE_PROPERTY ("short", "int16", boost::int16_t);
            else TRY_TO_GENERATE_PROPERTY ("ushort", "uint16", boost::uint16_t);
            else TRY_TO_GENERATE_PROPERTY ("int", "int32", boost::int32_t);
            else TRY_TO_GENERATE_PROPERTY ("uint", "uint32", boost::uint32_t);
            else TRY_TO_GENERATE_PROPERTY ("float", "float32", float);
            else TRY_TO_GENERATE_PROPERTY ("double", "float64", double);
          }

          continue;
        }
        else if (keyword == "comment")
        {
          std::string str = iss.str();
          if (str.size() > 8)
          {
            std::copy (str.begin() + 8, str.end(), std::back_inserter (m_comments));
            m_comments += "\n";
          }
        }
        else if (keyword == "element")
        {
          std::string type;
          std::size_t number;
          if (!(iss >> type >> number))
          {
            std::cerr << "Error line " << lineNumber << " of file" << std::endl;
            return false;
          }

          m_elements.push_back (PLY_element(type, number));
        }
        // When end_header is reached, stop loop and begin reading points
        else if (keyword == "end_header")
          break;
      }
    }
    return true;
  }

  ~PLY_reader ()
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
void process_properties (PLY_element& element, OutputValueType& new_element,
                         std::tuple<PropertyMap, Constructor, PLY_property<T>...>&& current)
{
  typedef typename PropertyMap::value_type PmapValueType;
  std::tuple<T...> values;
  Filler<sizeof...(T)-1>::fill(element, values, current);
  PmapValueType new_value = call_functor<PmapValueType>(std::get<1>(current), values);
  put (std::get<0>(current), new_element, new_value);
}

template <typename OutputValueType,
          typename PropertyMap,
          typename Constructor,
          typename ... T,
          typename NextPropertyBinder,
          typename ... PropertyMapBinders>
void process_properties (PLY_element& element, OutputValueType& new_element,
                         std::tuple<PropertyMap, Constructor, PLY_property<T>...>&& current,
                         NextPropertyBinder&& next,
                         PropertyMapBinders&& ... properties)
{
  typedef typename PropertyMap::value_type PmapValueType;
  std::tuple<T...> values;
  Filler<sizeof...(T)-1>::fill(element, values, current);
  PmapValueType new_value = call_functor<PmapValueType>(std::get<1>(current), values);
  put (std::get<0>(current), new_element, new_value);

  process_properties (element, new_element, std::forward<NextPropertyBinder>(next),
                      std::forward<PropertyMapBinders>(properties)...);
}


template <typename OutputValueType, typename PropertyMap, typename T>
void process_properties (PLY_element& element, OutputValueType& new_element,
                         std::pair<PropertyMap, PLY_property<T> >&& current)
{
  T new_value = T();
  element.assign (new_value, current.second.name);
  put (current.first, new_element, new_value);
}

template <typename OutputValueType, typename PropertyMap, typename T,
          typename NextPropertyBinder, typename ... PropertyMapBinders>
void process_properties (PLY_element& element, OutputValueType& new_element,
                         std::pair<PropertyMap, PLY_property<T> >&& current,
                         NextPropertyBinder&& next,
                         PropertyMapBinders&& ... properties)
{
  T new_value = T();
  element.assign (new_value, current.second.name);
  put (current.first, new_element, new_value);
  process_properties (element, new_element, std::forward<NextPropertyBinder>(next),
                      std::forward<PropertyMapBinders>(properties)...);
}

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
void output_property_header (std::ostream& stream,
                             std::tuple<PropertyMap, PLY_property<T>... >&& current)
{
  Properties_header<sizeof...(T)-1>::write(stream, current);
}


template <typename PropertyMap,
          typename T>
void output_property_header (std::ostream& stream,
                             std::pair<PropertyMap, PLY_property<T> >&& current)
{
  property_header (stream, current.second);
}

template <typename PropertyMap,
          typename T,
          typename NextPropertyHandler,
          typename ... PropertyHandler>
void output_property_header (std::ostream& stream,
                             std::pair<PropertyMap, PLY_property<T> >&& current,
                             NextPropertyHandler&& next,
                             PropertyHandler&& ... properties)
{
  property_header (stream, current.second);
  output_property_header (stream, std::forward<NextPropertyHandler>(next),
                          std::forward<PropertyHandler>(properties)...);
}
template <typename PropertyMap,
          typename ... T,
          typename NextPropertyHandler,
          typename ... PropertyHandler>
void output_property_header (std::ostream& stream,
                             std::tuple<PropertyMap, PLY_property<T>... >&& current,
                             NextPropertyHandler&& next,
                             PropertyHandler&& ... properties)
{
  Properties_header<sizeof...(T)-1>::write(stream, current);
  output_property_header (stream, std::forward<NextPropertyHandler>(next),
                          std::forward<PropertyHandler>(properties)...);
}


template <typename ForwardIterator,
          typename PropertyMap>
void property_write (std::ostream& stream, ForwardIterator it, PropertyMap map)
{
  stream << CGAL::oformat(get (map, *it));
}

template <typename T>
inline T no_char_character (const T& t) { return t; }
inline int no_char_character (const char& t) { return int(t); }
inline int no_char_character (const signed char& t) { return int(t); }
inline int no_char_character (const unsigned char& t) { return int(t); }

template <typename ForwardIterator,
          typename PropertyMap,
          typename T>
void simple_property_write (std::ostream& stream, ForwardIterator it,
                            std::pair<PropertyMap, PLY_property<T> > map)
{
  if (CGAL::get_mode(stream) == IO::ASCII)
    stream << no_char_character(get (map.first, *it));
  else
  {
    typename PropertyMap::value_type value = get(map.first, *it);
    stream.write (reinterpret_cast<char*>(&value), sizeof(value));
  }
}

template <typename ForwardIterator,
          typename PropertyMap,
          typename T>
void simple_property_write (std::ostream& stream, ForwardIterator it,
                            std::pair<PropertyMap, PLY_property<std::vector<T> > > map)
{
  const typename PropertyMap::reference value = get(map.first, *it);

  if (CGAL::get_mode(stream) == IO::ASCII)
  {
    stream << value.size();
    for (std::size_t i = 0; i < value.size(); ++ i)
      stream << " " << no_char_character(value[i]);
  }
  else
  {
    unsigned char size = static_cast<unsigned char>(value.size());
    stream.write (reinterpret_cast<char*>(&size), sizeof(size));
    for (std::size_t i = 0; i < value.size(); ++ i)
    {
      T t = T(value[i]);
      stream.write (reinterpret_cast<char*>(&t), sizeof(t));
    }
  }
}


template <typename ForwardIterator,
          typename PropertyMap,
          typename ... T>
void output_properties (std::ostream& stream,
                        ForwardIterator it,
                        std::tuple<PropertyMap, PLY_property<T>... >&& current)
{
  property_write (stream, it, std::get<0>(current));
  if (get_mode(stream) == IO::ASCII)
    stream << std::endl;
}


template <typename ForwardIterator,
          typename PropertyMap,
          typename T>
void output_properties (std::ostream& stream,
                        ForwardIterator it,
                        std::pair<PropertyMap, PLY_property<T> >&& current)
{
  simple_property_write (stream, it, std::forward<std::pair<PropertyMap, PLY_property<T> > >(current));
  if (get_mode(stream) == IO::ASCII)
    stream << std::endl;
}

template <typename ForwardIterator,
          typename PropertyMap,
          typename T,
          typename NextPropertyHandler,
          typename ... PropertyHandler>
void output_properties (std::ostream& stream,
                        ForwardIterator it,
                        std::pair<PropertyMap, PLY_property<T> >&& current,
                        NextPropertyHandler&& next,
                        PropertyHandler&& ... properties)
{
  simple_property_write (stream, it, current);
  if (get_mode(stream) == IO::ASCII)
    stream << " ";
  output_properties (stream, it, std::forward<NextPropertyHandler>(next),
                     std::forward<PropertyHandler>(properties)...);
}

template <typename ForwardIterator,
          typename PropertyMap,
          typename ... T,
          typename NextPropertyHandler,
          typename ... PropertyHandler>
void output_properties (std::ostream& stream,
                        ForwardIterator it,
                        std::tuple<PropertyMap, PLY_property<T>... >&& current,
                        NextPropertyHandler&& next,
                        PropertyHandler&& ... properties)
{
  property_write (stream, it, std::get<0>(current));
  if (get_mode(stream) == IO::ASCII)
    stream << " ";
  output_properties (stream, it, std::forward<NextPropertyHandler>(next),
                     std::forward<PropertyHandler>(properties)...);
}


// Printer classes used by Point_set_3 and Surface_mesh (translate a
// property map to a PLY property)

template <typename Index>
class Abstract_property_printer
{
public:
  virtual ~Abstract_property_printer() { }
  virtual void print (std::ostream& stream, const Index& index) = 0;
};

template <typename Index, typename PropertyMap>
class Property_printer : public Abstract_property_printer<Index>
{
  PropertyMap m_pmap;
public:
  Property_printer (const PropertyMap& pmap) : m_pmap (pmap)
  {

  }

  virtual void print(std::ostream& stream, const Index& index)
  {
    stream << get(m_pmap, index);
  }
};

template <typename Index, typename PropertyMap,
          typename Type = typename PropertyMap::value_type>
class Simple_property_printer : public Abstract_property_printer<Index>
{
  PropertyMap m_pmap;
public:
  Simple_property_printer (const PropertyMap& pmap) : m_pmap (pmap)
  {

  }

  virtual void print(std::ostream& stream, const Index& index)
  {
    if (get_mode(stream) == IO::ASCII)
      stream << get(m_pmap, index);
    else
    {
      Type t = Type(get (m_pmap, index));
      stream.write (reinterpret_cast<char*>(&t), sizeof(t));
    }
  }
};

template <typename Index, typename PropertyMap>
class Char_property_printer : public Abstract_property_printer<Index>
{
  typedef typename PropertyMap::value_type Type;
  PropertyMap m_pmap;
public:
  Char_property_printer (const PropertyMap& pmap) : m_pmap (pmap)
  {

  }

  virtual void print(std::ostream& stream, const Index& index)
  {
    if (get_mode(stream) == IO::ASCII)
      stream << int(get(m_pmap, index));
    else
    {
      Type t = get (m_pmap, index);
      stream.write (reinterpret_cast<char*>(&t), sizeof(t));
    }
  }
};

} // namespace PLY

} // namespace internal

} // namespace CGAL


#endif // CGAL_IO_PLY_H
