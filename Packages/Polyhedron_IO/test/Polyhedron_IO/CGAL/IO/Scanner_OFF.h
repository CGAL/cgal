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
// file          : CGAL_Scanner_OFF.h
// package       : $CGAL_Package: $
// source        : polyhedron_io.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : Herve Bronnimann  <Herve.Bronnimann@sophia.inria.fr>
//
// STL compliant file scanner for an object in object file format (OFF).
// ============================================================================

#ifndef CGAL_IO_SCANNER_OFF_H
#define CGAL_IO_SCANNER_OFF_H 1
#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif
#ifndef CGAL_PROTECT_ITERATOR_H
#include <iterator.h>
#define CGAL_PROTECT_ITERATOR_H
#endif // CGAL_PROTECT_ITERATOR_H
#ifndef CGAL_PROTECT_PAIR_H
#include <pair.h>
#define CGAL_PROTECT_PAIR_H
#endif // CGAL_PROTECT_PAIR_H
#ifndef CGAL_IO_FILE_SCANNER_OFF_H
#include <CGAL/IO/File_scanner_OFF.h>
#endif // CGAL_IO_FILE_SCANNER_OFF_H

template < class Pt >
class CGAL__Vertex_iterator_OFF : public input_iterator< Pt, int> {
    CGAL_File_scanner_OFF& _scanner;
    int                   _cnt;
    Pt                    _point;
public:
    typedef Pt                            Point;

    CGAL__Vertex_iterator_OFF( CGAL_File_scanner_OFF& scanner, int i)
        : _scanner( scanner), _cnt(i) {
        if ( _cnt < _scanner.size_of_vertices()) {
            CGAL_file_scan_vertex( _scanner, _point);
            _scanner.skip_to_next_vertex( _cnt);
        }
    }
    int
    count() const { return _cnt; }
    bool
    operator==( const CGAL__Vertex_iterator_OFF<Pt>& i) const {
        return _cnt == i._cnt;
    }
    bool
    operator!=( const CGAL__Vertex_iterator_OFF<Pt>& i) const {
        return _cnt != i._cnt;
    }
    CGAL__Vertex_iterator_OFF<Pt>&
    operator++() {
        CGAL_assertion( _cnt >= 0 && _cnt < _scanner.size_of_vertices());
        _cnt++;
        if ( _cnt < _scanner.size_of_vertices()) {
            CGAL_file_scan_vertex( _scanner, _point);
            _scanner.skip_to_next_vertex( _cnt);
        }
        return *this;
    }
    CGAL__Vertex_iterator_OFF<Pt>
    operator++(int) {
        CGAL_assertion( _cnt >= 0 && _cnt < _scanner.size_of_vertices());
        CGAL__Vertex_iterator_OFF<Pt> tmp = *this;
        ++(*this);
        return tmp;
    }
    const Pt&
    operator*() const {
        CGAL_assertion( _cnt >= 0 && _cnt < _scanner.size_of_vertices());
        return _point;
    }
// #ifdef  CGAL_CFG_ARROW_OPERATOR
    const Pt*
    operator->() const { return & operator*(); }
// #endif
};

class CGAL__Facet_iterator_OFF;

class CGAL__Facet_index_iterator_OFF : public input_iterator< int, int> {
    friend class CGAL__Facet_iterator_OFF;

    CGAL_File_scanner_OFF& _scanner;
    int                   _cnt;
    int                   _max;
    int                   _facet;
    int                   _index;

    void scan_next() {
        if ( _cnt < _max)
            _scanner.scan_facet_vertex_index( _index, _facet);
    }

public:
    CGAL__Facet_index_iterator_OFF( CGAL_File_scanner_OFF& scanner,
                                   int i, int max, int f)
        : _scanner( scanner), _cnt(i), _max(max), _facet(f) {}
    int
    count() const { return _cnt; }
    bool
    operator==( const CGAL__Facet_index_iterator_OFF& i) const {
        return _cnt == i._cnt;
    }
    bool
    operator!=( const CGAL__Facet_index_iterator_OFF& i) const {
        return _cnt != i._cnt;
    }
    CGAL__Facet_index_iterator_OFF&
    operator++() {
        CGAL_assertion( _cnt >= 0 && _cnt < _max);
        _cnt++;
        scan_next();
        return *this;
    }
    CGAL__Facet_index_iterator_OFF
    operator++(int) {
        CGAL_assertion( _cnt >= 0 && _cnt < _max);
        CGAL__Facet_index_iterator_OFF tmp = *this;
        ++(*this);
        return tmp;
    }
    const int&
    operator*() const {
        CGAL_assertion( _cnt >= 0 && _cnt < _max);
        return _index;
    }
};

class CGAL__Facet_iterator_OFF : public
        input_iterator< pair< CGAL__Facet_index_iterator_OFF,
                              CGAL__Facet_index_iterator_OFF>, int> {
public:
    typedef pair< CGAL__Facet_index_iterator_OFF,
                  CGAL__Facet_index_iterator_OFF>  Pair;
private:
    CGAL_File_scanner_OFF& _scanner;
    int                   _cnt;
    Pair                  _pair;

public:
    CGAL__Facet_iterator_OFF( CGAL_File_scanner_OFF& scanner, int i)
    : _scanner( scanner), _cnt(i), _pair(
        CGAL__Facet_index_iterator_OFF( _scanner, 0, 0, _cnt),
        CGAL__Facet_index_iterator_OFF( _scanner, 0, 0, _cnt)
    ){
        if ( _cnt < _scanner.size_of_facets()) {
            int _no;
            _scanner.scan_facet( _no, _cnt);
            _pair.first._max  = _no;
            _pair.second._cnt = _no;
            _pair.second._max = _no;
            _pair.first.scan_next();
        }
    }
    int
    count() const { return _cnt; }
    bool
    operator==( const CGAL__Facet_iterator_OFF& i) const {
        return _cnt == i._cnt;
    }
    bool
    operator!=( const CGAL__Facet_iterator_OFF& i) const {
        return _cnt != i._cnt;
    }
    CGAL__Facet_iterator_OFF&
    operator++() {
        CGAL_assertion( _cnt >= 0 && _cnt < _scanner.size_of_facets());
        _cnt++;
        if ( _cnt < _scanner.size_of_facets()) {
            int _no;
            _scanner.skip_to_next_facet( _cnt);
            _scanner.scan_facet( _no, _cnt);
            _pair.first._max  = _no;
            _pair.second._cnt = _no;
            _pair.second._max = _no;
            _pair.first.scan_next();
        }
        return *this;
    }
    CGAL__Facet_iterator_OFF
    operator++(int) {
        CGAL_assertion( _cnt >= 0 && _cnt < _scanner.size_of_facets());
        CGAL__Facet_iterator_OFF tmp = *this;
        ++(*this);
        return tmp;
    }
    const Pair&
    operator*() const {
        CGAL_assertion( _cnt >= 0 && _cnt < _scanner.size_of_facets());
        return _pair;
    }
// #ifdef  CGAL_CFG_ARROW_OPERATOR
    const Pair*
    operator->() const { return & operator*(); }
// #endif
};


// The distance function is implemented to work in constant time
// for all these three iterators.

#ifdef CGAL_PARTIAL_SPECIALISATION
template < class Pt >
void distance( const CGAL__Vertex_iterator_OFF<Pt>& first,
               const CGAL__Vertex_iterator_OFF<Pt>& last, int& n) {
    n = last.count() - first.count();
}
#endif
void distance( const CGAL__Facet_iterator_OFF& first,
               const CGAL__Facet_iterator_OFF& last, int& n) {
    n = last.count() - first.count();
}
void distance( const CGAL__Facet_index_iterator_OFF& first,
               const CGAL__Facet_index_iterator_OFF& last, int& n) {
    n = last.count() - first.count();
}

void distance( const CGAL__Facet_iterator_OFF& first,
               const CGAL__Facet_iterator_OFF& last, size_t& n) {
    n = size_t(last.count() - first.count());
}
void distance( const CGAL__Facet_index_iterator_OFF& first,
               const CGAL__Facet_index_iterator_OFF& last, size_t& n) {
    n = size_t(last.count() - first.count());
}


template < class Pt >
class CGAL_Scanner_OFF {
    CGAL_File_scanner_OFF  _scanner;
public:
    typedef Pt                            Point;
    typedef CGAL__Vertex_iterator_OFF<Pt>  Vertex_iterator;
    typedef CGAL__Facet_iterator_OFF       Facet_iterator;
    typedef CGAL__Facet_index_iterator_OFF Facet_index_iterator;

    CGAL_Scanner_OFF( istream& in, bool verbose = false)
        : _scanner( in, verbose) {}
    CGAL_Scanner_OFF( istream& in, const CGAL_File_header_OFF& header)
        : _scanner( in, header) {}
    bool verbose() const            { return _scanner.verbose();    }
    bool is_SKEL() const            { return _scanner.is_SKEL();    }
    bool is_OFF()  const            { return _scanner.is_OFF();     }
    int  size_of_vertices()  const  { return _scanner.size_of_vertices(); }
    int  size_of_halfedges() const  { return _scanner.size_of_halfedges();}
    int  size_of_facets()    const  { return _scanner.size_of_facets();   }
    bool has_colors() const         { return _scanner.has_colors(); }
    bool has_normals() const        { return _scanner.has_normals();}
    bool is_binary() const          { return _scanner.is_binary();  }
    CGAL_File_info&
         file_info()                { return _scanner.file_info();  }
    const CGAL_File_info&
         file_info() const          { return _scanner.file_info();  }

    Vertex_iterator  vertices_begin() {
                         return Vertex_iterator( _scanner, 0);
    }
    Vertex_iterator  vertices_end() {
                         return Vertex_iterator( _scanner,
                                                 size_of_vertices());
    }
    Facet_iterator   facets_begin() {
                         return Facet_iterator( _scanner, 0);
    }
    Facet_iterator   facets_end() {
                         return Facet_iterator( _scanner,
                                                size_of_facets());
    }
};
#endif // CGAL_IO_SCANNER_OFF_H //
// EOF //
