// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : Generic_writer.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// Generic STL compliant interface to write boundary rep file formats.
// ============================================================================

#ifndef CGAL_IO_GENERIC_WRITER_H
#define CGAL_IO_GENERIC_WRITER_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_PROTECT_ITERATOR_H
#include <iterator.h>
#define CGAL_PROTECT_ITERATOR_H
#endif // CGAL_PROTECT_ITERATOR_H

template < class Writer, class Pt >
class CGAL__Vertex_output_proxy {
    Writer&  _writer;
public:
    CGAL__Vertex_output_proxy( Writer& writer) : _writer(writer) {}
    void
    operator= ( const Pt& p) {
#ifdef CGAL_IO_GENERIC_WRITER
    _writer.write_vertex( CGAL_to_double(p.x()),
                          CGAL_to_double(p.y()),
                          CGAL_to_double(p.z()));
#else
        _writer.write_vertex( double(p.x()), double(p.y()), double(p.z()));
#endif
    }
};

template < class Writer, class Pt >
class CGAL__Vertex_output_iterator : public output_iterator {
    Writer&  _writer;
    size_t   _cnt;
    size_t   _max;
public:
    typedef Pt   Point;
    typedef CGAL__Vertex_output_proxy< Writer, Pt>  Proxy;

    CGAL__Vertex_output_iterator( Writer& writer, size_t v)
        : _writer( writer), _cnt(0), _max(v) {}
    CGAL__Vertex_output_iterator< Writer, Pt>&
    operator++() {
        CGAL_assertion( _cnt < _max);
        _cnt++;
        return *this;
    }
    CGAL__Vertex_output_iterator< Writer, Pt>&
    operator++(int) {
        CGAL_assertion( _cnt < _max);
        ++(*this);
        return *this;
    }
    Proxy
    operator*() const {
        CGAL_assertion( _cnt <= _max);
        return Proxy( _writer);
    }
};

template < class Writer >  class CGAL__Facet_output_iterator;

template < class Writer >
class CGAL__Facet_output_proxy {
    const CGAL__Facet_output_iterator< Writer> &  _foi;
public:
    CGAL__Facet_output_proxy( const CGAL__Facet_output_iterator<Writer>&
                             foi)
        : _foi(foi) {}
    inline
    void  operator= ( size_t i);
};


template < class Writer >
class CGAL__Facet_output_iterator : public output_iterator {
    friend class CGAL__Facet_output_proxy<Writer>;

    Writer& _writer;
    size_t  _fcnt;
    size_t  _icnt;

public:
    typedef  CGAL__Facet_output_proxy<Writer>  Proxy;
    CGAL__Facet_output_iterator( Writer& writer, size_t f)
    : _writer( writer), _fcnt(f), _icnt(0)
    {
        if ( f == 0)
            _writer.footer();
    }
    CGAL__Facet_output_iterator< Writer >&
    operator++()    { return *this; }
    CGAL__Facet_output_iterator< Writer >&
    operator++(int) { return *this; }
    Proxy
    operator*() const {
        CGAL_assertion( _fcnt > 0);
        return Proxy( *this);
    }
};

template < class Writer >
inline
void  CGAL__Facet_output_proxy< class Writer >::
operator= ( size_t i) {
    if (_foi._icnt == 0) {
        _foi._writer.write_facet_begin( i);
        ((CGAL__Facet_output_iterator< Writer >&)(_foi))._icnt = i;
    } else {
        _foi._writer.write_facet_vertex_index( i);
        ((CGAL__Facet_output_iterator< Writer >&)(_foi))._icnt --;
        if (_foi._icnt == 0) {
            _foi._writer.write_facet_end();
            ((CGAL__Facet_output_iterator< Writer >&)(_foi))._fcnt --;
            if (_foi._fcnt == 0)
                _foi._writer.footer();
        }
    }
}


template < class Writer, class Pt >
class CGAL_Generic_writer {
    Writer  _writer;
    size_t  _vertices;
    size_t  _halfedges;
    size_t  _facets;

public:
    typedef Pt                                      Point;
    typedef CGAL__Vertex_output_iterator<Writer,Pt>  Vertex_iterator;
    typedef CGAL__Facet_output_iterator<Writer>      Facet_iterator;

    CGAL_Generic_writer( const Writer& writer, ostream& out,
                        size_t v, size_t h, size_t f)
    : _writer( writer), _vertices(v), _halfedges(h), _facets(f) {
        _writer.header( out, v, h, f);
    }
    const Writer& writer()     const  { return _writer;    }
    size_t size_of_vertices()  const  { return _vertices;  }
    size_t size_of_halfedges() const  { return _halfedges; }
    size_t size_of_facets()    const  { return _facets;    }

    Vertex_iterator  vertices_begin() {
                         return Vertex_iterator( _writer, _vertices);
    }
    Facet_iterator   facets_begin() {
                         _writer.write_facet_header();
                         return Facet_iterator( _writer, _facets);
    }
};
#endif // CGAL_IO_GENERIC_WRITER_H //
// EOF //
