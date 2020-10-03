// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
// Copyright (c) 2011  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_IO_SIGNATURE_H
#define CGAL_MESH_3_IO_SIGNATURE_H

#include <CGAL/license/Triangulation_3.h>

#define CGAL_MESH_3_IO_H // the old include macro, tested by other files

#include <CGAL/Point_3.h>
#include <CGAL/Weighted_point_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_with_weighted_circumcenter_3.h>

#ifdef CGAL_PERIODIC_3_MESH_3_CONFIG_H
#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/Periodic_3_regular_triangulation_3.h>
#endif

#include <boost/variant.hpp>
#include <boost/tuple/tuple.hpp>
#include <utility>

namespace CGAL {

// SFINAE test
template <typename T, typename U>
class has_io_signature
{
private:
  template <U> struct helper;
  template <typename V> static char check(helper<&V::io_signature> *);
  template <typename V> static char (&check(...))[2];

public:
  enum { value = (sizeof(check<T>(0)) == sizeof(char)) };
};

template <class T, bool has_io_signature>
struct Get_io_signature_aux
{
  std::string operator() () const
  {
    return T::io_signature();
  }
}; // end struct template Get_io_signature_aux

template <class T>
struct Get_io_signature_aux<T, false>
{
  std::string operator()() const
  {
    std::cerr << "Type without signature: " << typeid(T).name() << std::endl;
    return std::string();
  }
}; // end template partial specialization Get_io_signature_aux<T, false>


template <class T>
struct Get_io_signature
  : public Get_io_signature_aux<
  T,
  (has_io_signature<T, std::string (T::*)() >::value ||
   has_io_signature<T, std::string (T::*)() const >::value ||
   has_io_signature<T, std::string (*)() >::value )  // signature for
                                                     // static mem func
  >
{
};

template <>
struct Get_io_signature<int>
{
  std::string operator()() {
    return "i";
  }
};

template <>
struct Get_io_signature<unsigned int>
{
  std::string operator()() {
    return "ui";
  }
};

template <>
struct Get_io_signature<char>
{
  std::string operator()() {
    return "c";
  }
};

template <>
struct Get_io_signature<unsigned char>
{
  std::string operator()() {
    return "uc";
  }
};

template <>
struct Get_io_signature<signed char>
{
  std::string operator()() {
    return "sc";
  }
};

template <>
struct Get_io_signature<short int>
{
  std::string operator()() {
    return "s";
  }
};

template <>
struct Get_io_signature<unsigned short int>
{
  std::string operator()() {
    return "us";
  }
};

template <>
struct Get_io_signature<double>
{
  std::string operator()() {
    return "d";
  }
};

template <typename T, typename U>
struct Get_io_signature<boost::variant<T,U> >
{
    std::string operator()() {
      return std::string("boost::variant<") +
        Get_io_signature<T>()() + "," +
        Get_io_signature<U>()() + ">";
  }
};

template <typename T, typename U>
struct Get_io_signature<std::pair<T,U> >
{
    std::string operator()() {
      return std::string("std::pair<") +
        Get_io_signature<T>()() + "," +
        Get_io_signature<U>()() + ">";
  }
};

template <typename T, typename U>
struct Get_io_signature<boost::tuple<T,U> >
{
    std::string operator()() {
      return std::string("std::pair<") +
        Get_io_signature<T>()() + "," +
        Get_io_signature<U>()() + ">";
  }
};

template <typename T, typename U, typename V>
struct Get_io_signature<boost::variant<T,U,V> >
{
    std::string operator()() {
      return std::string("boost::variant<") +
        Get_io_signature<T>()() + "," +
        Get_io_signature<U>()() + "," +
        Get_io_signature<V>()() + ">";
  }
};

template <typename T, typename U,
          typename V, typename W>
struct Get_io_signature<boost::variant<T,U,V,W> >
{
    std::string operator()() {
      return std::string("boost::variant<") +
        Get_io_signature<T>()() + "," +
        Get_io_signature<U>()() + "," +
        Get_io_signature<V>()() + "," +
        Get_io_signature<W>()() + ">";
  }
};

template <class Kernel>
struct Get_io_signature<Point_3<Kernel> >
{
  std::string operator()() {
    return "Point_3";
  }
};

template <class K>
struct Get_io_signature<Weighted_point_3<K> >
{
  std::string operator()() {
    return std::string("Weighted_point<") + Get_io_signature<Point_3<K> >()() + ">";
  }
};

#ifdef CGAL_TRIANGULATION_3_H
template <class Gt, class Tds>
struct
Get_io_signature<Triangulation_3<Gt, Tds > >
{
  std::string operator()() {
    return std::string("Triangulation_3(") +
      Get_io_signature<typename Tds::Vertex::Point>()() +
      ",Vb(" + Get_io_signature<typename Tds::Vertex>()() +
      "),Cb(" + Get_io_signature<typename Tds::Cell>()() +
      "))";
  }
};
#endif

#ifdef CGAL_DELAUNAY_TRIANGULATION_3_H
template <class Gt, class Tds>
struct
Get_io_signature<Delaunay_triangulation_3<Gt, Tds> >
{
  std::string operator()() {
    return Get_io_signature<Triangulation_3<Gt, Tds> >()();
  }
};
#endif

#ifdef CGAL_REGULAR_TRIANGULATION_3_H
template <class Gt, class Tds>
struct
Get_io_signature<Regular_triangulation_3<Gt, Tds> >
{
  std::string operator()() {
    return Get_io_signature<Triangulation_3<Gt, Tds> >()();
  }
};
#endif

#ifdef CGAL_PERIODIC_3_TRIANGULATION_3_H
template <class Gt, class Tds>
struct
Get_io_signature<Periodic_3_triangulation_3<Gt, Tds> >
{
  std::string operator()() {
    return std::string("Periodic_3_triangulation_3(") +
      Get_io_signature<typename Tds::Vertex::Point>()() +
      ",Vb(" + Get_io_signature<typename Tds::Vertex>()() +
      "),Cb(" + Get_io_signature<typename Tds::Cell>()() +
      "))";
  }
};
#endif

#ifdef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_3_H
template <class Gt, class Tds>
struct
Get_io_signature<Periodic_3_regular_triangulation_3<Gt, Tds> >
{
  std::string operator()() {
    return Get_io_signature<Periodic_3_triangulation_3<Gt, Tds> >()();
  }
};
#endif

#ifdef CGAL_TRIANGULATION_VERTEX_BASE_3_H
template <class Gt, class Vb>
struct Get_io_signature<Triangulation_vertex_base_3<Gt, Vb> >
{
  std::string operator()() {
    return "Tvb_3";
  }
};
#endif

#ifdef CGAL_REGULAR_TRIANGULATION_VERTEX_BASE_3_H
template <class Gt, class Vb>
struct Get_io_signature<Regular_triangulation_vertex_base_3<Gt, Vb> >
{
  // identical to Triangulation_vertex_base_3
  std::string operator()() {
    return "Tvb_3";
  }
};
#endif

#ifdef CGAL_TRIANGULATION_VERTEX_BASE_WITH_INFO_3_H
template <class Info, class Gt, class Vb>
struct
Get_io_signature<Triangulation_vertex_base_with_info_3<Info, Gt, Vb> >
{
  std::string operator()() {
    return Get_io_signature<Vb>()();
  }
};
#endif

#ifdef CGAL_TRIANGULATION_CELL_BASE_3_H
template <class Gt, class Cb>
struct
Get_io_signature<Triangulation_cell_base_3<Gt, Cb> >
{
  std::string operator()() {
    return "Tcb_3";
  }
};
#endif

#ifdef CGAL_TRIANGULATION_CELL_BASE_WITH_INFO_3_H
template <class Info, class Gt, class Cb>
struct
Get_io_signature<Triangulation_cell_base_with_info_3<Info, Gt, Cb> >
{
  std::string operator()() {
    return Get_io_signature<Cb>()();
  }
};
#endif

#ifdef CGAL_REGULAR_TRIANGULATION_CELL_BASE_3_H
template <class Gt, class Cb, class Container>
struct
Get_io_signature<Regular_triangulation_cell_base_3<Gt, Cb, Container> >
{
  std::string operator()() {
    return "RTcb_3";
  }
};
#endif

#ifdef CGAL_REGULAR_TRIANGULATION_CELL_BASE_WITH_CIRCUMCENTER_3_H
template <class Gt, class Cb>
struct
Get_io_signature<Regular_triangulation_cell_base_with_weighted_circumcenter_3<Gt, Cb> >
{
  std::string operator()() {
    return "RTWWCcb_3";
  }
};
#endif

} // end namespace CGAL


#endif // CGAL_MESH_3_IO_SIGNATURE_H
