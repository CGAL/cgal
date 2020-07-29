// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_GENERIC_WRITER_H
#define CGAL_IO_GENERIC_WRITER_H 1

#include <CGAL/basic.h>
#include <iterator>

namespace CGAL {

template <class Writer>
class I_Generic_writer_vertex_proxy {
    Writer&  m_writer;
public:
    typedef typename Writer::Point Point;
    I_Generic_writer_vertex_proxy( Writer& w) : m_writer(w) {}
    void operator= ( const Point& p) { m_writer.write_vertex(p); }
};

template <class Writer>
class I_Generic_writer_vertex_iterator
{
    Writer&  m_writer;
public:
    typedef std::output_iterator_tag  iterator_category;
    typedef typename Writer::Point    value_type;
    typedef std::ptrdiff_t            difference_type;
    typedef value_type*               pointer;
    typedef value_type&               reference;

    typedef I_Generic_writer_vertex_proxy< Writer>    Proxy;
    typedef I_Generic_writer_vertex_iterator< Writer> Self;

    I_Generic_writer_vertex_iterator( Writer& w) : m_writer(w) {}
    Self& operator++()      { return *this; }
    Self& operator++(int)   { return *this; }
    Proxy operator*() const { return Proxy( m_writer); }
};

template <class Writer>
class I_Generic_writer_facet_proxy {
    Writer&  m_writer;
public:
    I_Generic_writer_facet_proxy( Writer& w) : m_writer(w) {}
    void operator= ( std::size_t i) { m_writer.write_facet_index(i); }
};

template <class Writer>
class I_Generic_writer_facet_iterator
{
    Writer& m_writer;
public:
    typedef std::output_iterator_tag  iterator_category;
    typedef std::size_t               value_type;
    typedef std::ptrdiff_t            difference_type;
    typedef value_type*               pointer;
    typedef value_type&               reference;

    typedef  I_Generic_writer_facet_proxy<Writer>    Proxy;
    typedef  I_Generic_writer_facet_iterator<Writer> Self;

    I_Generic_writer_facet_iterator( Writer& w) : m_writer(w) {}
    Self& operator++()      { return *this; }
    Self& operator++(int)   { return *this; }
    Proxy operator*() const { return Proxy( m_writer); }
};

// The Generic_writer class contains also the state necessary for all
// its iterators. Since these iterators are of the category
// output_iterator, they could not be used more than once and their
// positional state can be kept in Generic_writer, which simplifies
// the situation where the iterators are copied by value.

template < class Writer, class Pt >
class Generic_writer {
    Writer       m_writer;
    std::size_t  m_vertices;
    std::size_t  m_halfedges;
    std::size_t  m_facets;

    std::size_t  m_vcnt;
    std::size_t  m_fcnt;
    std::size_t  m_icnt;
public:
    typedef Pt                                      Point;
    typedef Generic_writer< Writer, Pt>             Self;
    typedef I_Generic_writer_vertex_iterator<Self>  Vertex_iterator;
    typedef I_Generic_writer_facet_iterator<Self>   Facet_iterator;

    Generic_writer( const Writer& writer, std::ostream& out,
                    std::size_t v, std::size_t h, std::size_t f)
        : m_writer( writer), m_vertices(v), m_halfedges(h), m_facets(f),
          m_vcnt(0), m_fcnt(0), m_icnt(0)
    {
        m_writer.write_header( out, v, h, f);
    }
    const Writer& writer()            const  { return m_writer;    }
    std::size_t   size_of_vertices()  const  { return m_vertices;  }
    std::size_t   size_of_halfedges() const  { return m_halfedges; }
    std::size_t   size_of_facets()    const  { return m_facets;    }

    Vertex_iterator  vertices_begin() { return Vertex_iterator( *this); }
    Facet_iterator   facets_begin() {
        if ( m_vcnt != m_vertices) {
            std::cerr << "error: Generic_writer: wrong number of "
                         "vertices written, " << m_vcnt << " instead of "
                      << m_vertices << "." << std::endl;
            m_writer.out().clear( std::ios::badbit);
        }
        m_writer.write_facet_header();
        if ( m_facets == 0)
            m_writer.write_footer();
        return Facet_iterator( *this);
    }

    // Interface used by the iterators and their proxies.
    void write_vertex( const Point& p) {
        ++m_vcnt;
        m_writer.write_vertex( ::CGAL::to_double( p.x()),
                               ::CGAL::to_double( p.y()),
                               ::CGAL::to_double( p.z()));
    }
    void write_facet_index( std::size_t i) {
        if ( m_fcnt > m_facets) {
            std::cerr << "error: Generic_writer: too many facets written."
                      << std::endl;
            m_writer.out().clear( std::ios::badbit);
        }
        if ( m_icnt == 0) {
            m_writer.write_facet_begin( i);
            m_icnt = i;
        } else {
            m_writer.write_facet_vertex_index( i);
            m_icnt --;
            if ( m_icnt == 0) {
                m_writer.write_facet_end();
                m_fcnt ++;
                if (m_fcnt == m_facets)
                    m_writer.write_footer();
            }
        }
    }
};

} //namespace CGAL
#endif // CGAL_IO_GENERIC_WRITER_H //
// EOF //
