// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_HALFEDGEDS_CONST_DECORATOR_H
#define CGAL_HALFEDGEDS_CONST_DECORATOR_H 1

#include <CGAL/HalfedgeDS_items_decorator.h>
#include <vector>
#include <CGAL/IO/Verbose_ostream.h>

namespace CGAL {

template < class p_HDS >
class HalfedgeDS_const_decorator
    : public HalfedgeDS_items_decorator<p_HDS> {

// TYPES (inherited from Items_decorator, but have to be repeated)
// ---------------------------------------------------------------
public:
    typedef p_HDS                                 HDS;
    typedef p_HDS                                 HalfedgeDS;
    typedef typename HDS::Traits                  Traits;
    typedef typename HDS::Vertex                  Vertex;
    typedef typename HDS::Halfedge                Halfedge;
    typedef typename HDS::Face                    Face;

    typedef typename HDS::Vertex_handle           Vertex_handle;
    typedef typename HDS::Vertex_const_handle     Vertex_const_handle;
    typedef typename HDS::Vertex_iterator         Vertex_iterator;
    typedef typename HDS::Vertex_const_iterator   Vertex_const_iterator;

    typedef typename HDS::Halfedge_handle         Halfedge_handle;
    typedef typename HDS::Halfedge_const_handle   Halfedge_const_handle;
    typedef typename HDS::Halfedge_iterator       Halfedge_iterator;
    typedef typename HDS::Halfedge_const_iterator Halfedge_const_iterator;

    typedef typename HDS::Face_handle             Face_handle;
    typedef typename HDS::Face_const_handle       Face_const_handle;
    typedef typename HDS::Face_iterator           Face_iterator;
    typedef typename HDS::Face_const_iterator     Face_const_iterator;

    typedef typename HDS::size_type               size_type;
    typedef typename HDS::difference_type         difference_type;
    typedef typename HDS::iterator_category       iterator_category;

// The following types are equal to either `Tag_true' or `Tag_false',
// dependent whether the named feature is supported or not.

    typedef typename HDS::Supports_vertex_halfedge
                                                  Supports_vertex_halfedge;
    typedef typename HDS::Supports_halfedge_prev  Supports_halfedge_prev;
    typedef typename HDS::Supports_halfedge_vertex
                                                  Supports_halfedge_vertex;
    typedef typename HDS::Supports_halfedge_face  Supports_halfedge_face;
    typedef typename HDS::Supports_face_halfedge  Supports_face_halfedge;

    typedef typename HDS::Supports_removal        Supports_removal;


    using HalfedgeDS_items_decorator<p_HDS>::insert_tip;
    using HalfedgeDS_items_decorator<p_HDS>::get_prev;
    using HalfedgeDS_items_decorator<p_HDS>::set_prev;
    using HalfedgeDS_items_decorator<p_HDS>::get_face;
    using HalfedgeDS_items_decorator<p_HDS>::set_face;
    using HalfedgeDS_items_decorator<p_HDS>::get_vertex;
    using HalfedgeDS_items_decorator<p_HDS>::set_vertex;
    using HalfedgeDS_items_decorator<p_HDS>::get_vertex_halfedge;
    using HalfedgeDS_items_decorator<p_HDS>::set_vertex_halfedge;
    using HalfedgeDS_items_decorator<p_HDS>::get_face_halfedge;
    using HalfedgeDS_items_decorator<p_HDS>::set_face_halfedge;
    using HalfedgeDS_items_decorator<p_HDS>::set_vertex_in_vertex_loop;
    using HalfedgeDS_items_decorator<p_HDS>::set_face_in_face_loop;
    using HalfedgeDS_items_decorator<p_HDS>::insert_halfedge;

protected:
    typedef typename Vertex::Base                 VBase;
    typedef typename Halfedge::Base               HBase;
    typedef typename Halfedge::Base_base          HBase_base;
    typedef typename Face::Base                   FBase;


// PRIVATE MEMBER VARIABLES
// ----------------------------------
private:
    const p_HDS*  hds;

// CREATION
    // ----------------------------------
public:
    // No default constructor, keeps always a reference to a HDS!

    HalfedgeDS_const_decorator( const p_HDS& h) : hds(&h) {}
        // keeps internally a const reference to `hds'.

    bool is_valid( bool verb = false, int level = 0) const;
        // returns `true' if the halfedge data structure `hds' is valid
        // with respect to the `level' value as defined in the reference
        // manual. If `verbose' is `true', statistics are written to
        // `cerr'.

    bool normalized_border_is_valid( bool verb = false) const;
        // returns `true' if the border halfedges are in normalized
        // representation, which is when enumerating all halfedges with
        // the halfedge iterator the following holds: The non-border edges
        // precede the border edges. For border edges, the second halfedge
        // is a border halfedge. (The first halfedge may or may not be a
        // border halfedge.) The halfedge iterator `border_halfedges_begin
        // ()' denotes the first border edge. If `verbose' is `true',
        // statistics are written to `cerr'.
};

template < class p_HDS >
bool
HalfedgeDS_const_decorator<p_HDS>::
normalized_border_is_valid( bool verb) const {
    bool valid = true;
    Verbose_ostream verr(verb);
    verr << "begin CGAL::HalfedgeDS_const_decorator<HDS>::"
            "normalized_border_is_valid( verb=true):" << std::endl;

    Halfedge_const_iterator e = hds->halfedges_begin();
    size_type count = 0;
    while( e != hds->halfedges_end() && ! e->is_border() && !
           e->opposite()->is_border()) {
        ++e;
        ++e;
        ++count;
    }
    verr << "    non-border edges: " << count << std::endl;
    if ( e != hds->border_halfedges_begin()) {
        verr << "    first border edge does not start at "
                "border_halfedges_begin()" << std::endl;
        valid = false;
    } else {
        count = 0;
        while( valid && e != hds->halfedges_end() &&
               e->opposite()->is_border()) {
            ++e;
            ++e;
            ++count;
        }
        verr << "    border     edges: " << count << std::endl;
        verr << "    total      edges: " << hds->size_of_halfedges()/2
             << std::endl;
        if ( e != hds->halfedges_end()) {
            if ( e->is_border()) {
                verr << "    border edge " << count
                     << ": wrong orientation." << std::endl;
            }
            verr << "    the sum of full + border equals not total edges."
                 << std::endl;
            valid = false;
        }
    }
    verr << "end of CGAL::HalfedgeDS_const_decorator<HDS>::normalized_"
            "border_is_valid(): structure is "
         << ( valid ? "valid." : "NOT VALID.") << std::endl;
    return valid;
}

template < class p_HDS >
bool
HalfedgeDS_const_decorator<p_HDS>::
is_valid( bool verb, int level) const {
    Verbose_ostream verr(verb);
    verr << "begin CGAL::HalfedgeDS_const_decorator<HDS>::is_valid("
            " verb=true, level = " << level << "):" << std::endl;

    bool valid = ( 1 != (hds->size_of_halfedges() & 1));
    if ( ! valid)
        verr << "number of halfedges is odd." << std::endl;

    // All halfedges.
    Halfedge_const_iterator begin = hds->halfedges_begin();
    Halfedge_const_iterator end   = hds->halfedges_end();
    size_type  n = 0;
    size_type nb = 0;
    for( ; valid && (begin != end); begin++) {
        verr << "halfedge " << n << std::endl;
        if ( begin->is_border())
            verr << "    is border halfedge" << std::endl;
        // Pointer integrity.
        valid = valid && ( begin->next() != Halfedge_const_handle());
        valid = valid && ( begin->opposite() != Halfedge_const_handle());
        if ( ! valid) {
            verr << "    pointer integrity corrupted (ptr==0)."
                 << std::endl;
            break;
        }
        // opposite integrity.
        valid = valid && ( begin->opposite() != begin);
        valid = valid && ( begin->opposite()->opposite() == begin);
        if ( ! valid) {
            verr << "    opposite pointer integrity corrupted."
                 << std::endl;
            break;
        }
        // previous integrity.
        valid = valid && ( ! check_tag( Supports_halfedge_prev()) ||
                           get_prev(begin->next()) == begin);
        if ( ! valid) {
            verr << "    previous pointer integrity corrupted."
                 << std::endl;
            break;
        }
        if ( level > 0) {
            // vertex integrity.
            valid = valid && ( ! check_tag( Supports_halfedge_vertex())
                          || get_vertex( begin) != Vertex_const_handle());
            if ( ! valid) {
                verr << "    vertex pointer integrity corrupted."
                     << std::endl;
                break;
            }
            valid = valid && ( get_vertex( begin) ==
                               get_vertex( begin->next()->opposite()));
            if ( ! valid) {
                verr << "    vertex pointer integrity2 corrupted."
                     << std::endl;
                break;
            }
            // face integrity.
            valid = valid && ( ! check_tag( Supports_halfedge_face())
                               || begin->is_border()
                               || get_face(begin) != Face_const_handle());
            if ( ! valid) {
                verr << "    face pointer integrity corrupted."
                     << std::endl;
                break;
            }
            valid = valid && ( get_face(begin) == get_face(begin->next()));
            if ( ! valid) {
                verr << "    face pointer integrity2 corrupted."
                     << std::endl;
                break;
            }
        }
        ++n;
        if ( begin->is_border())
            ++nb;
    }
    verr << "summe border halfedges (2*nb) = " << 2 * nb << std::endl;
    if ( valid && n != hds->size_of_halfedges())
        verr << "counting halfedges failed." << std::endl;
    if ( valid && level >= 4 && (nb != hds->size_of_border_halfedges()))
        verr << "counting border halfedges failed." << std::endl;
    valid = valid && ( n  == hds->size_of_halfedges());
    valid = valid && ( level < 4 ||
                       (nb == hds->size_of_border_halfedges()));

    // All vertices.
    Vertex_const_iterator vbegin = hds->vertices_begin();
    Vertex_const_iterator vend   = hds->vertices_end();
    size_type v = 0;
    n = 0;
    for( ; valid && (vbegin != vend); ++vbegin) {
        verr << "vertex " << v << std::endl;
        // Pointer integrity.
        if ( get_vertex_halfedge( vbegin) != Halfedge_const_handle())
            valid = valid && ( ! check_tag(
                Supports_halfedge_vertex()) ||
                get_vertex( get_vertex_halfedge( vbegin)) == vbegin);
        else
            valid = valid && (! check_tag(
                                Supports_vertex_halfedge()));
        if ( ! valid) {
            verr << "    halfedge pointer in vertex corrupted."
                 << std::endl;
            break;
        }
        // cycle-around-vertex test.
        Halfedge_const_handle h = get_vertex_halfedge( vbegin);
        if ( level >= 2 && h != Halfedge_const_handle()) {
            Halfedge_const_handle g = h;
            do {
                verr << "    halfedge " << n << std::endl;
                ++n;
                h = h->next()->opposite();
                valid = valid && ( n <= hds->size_of_halfedges() && n!=0);
                if ( ! valid)
                    verr << "    too many halfedges around vertices."
                         << std::endl;
            } while ( valid && (h != g));
        }
        ++v;
    }
    if ( valid && v != hds->size_of_vertices())
        verr << "counting vertices failed." << std::endl;
    if ( valid && level >= 2 && ( check_tag( Supports_vertex_halfedge())
         && n  != hds->size_of_halfedges()))
        verr << "counting halfedges via vertices failed." << std::endl;
    valid = valid && ( v == hds->size_of_vertices());
    valid = valid && ( level < 2 ||
                       ! check_tag( Supports_vertex_halfedge()) ||
                       n  == hds->size_of_halfedges());

    // All faces.
    Face_const_iterator fbegin = hds->faces_begin();
    Face_const_iterator fend   = hds->faces_end();
    size_type f = 0;
    n = 0;
    for( ; valid && (fbegin != fend); ++fbegin) {
        verr << "face " << f << std::endl;
        // Pointer integrity.
        if ( get_face_halfedge( fbegin) != Halfedge_const_handle())
            valid = valid && ( ! check_tag(
                Supports_halfedge_face()) ||
                get_face( get_face_halfedge( fbegin)) == fbegin);
        else
            valid = valid && (! check_tag(
                Supports_face_halfedge()) || begin->is_border());
        if ( ! valid) {
            verr << "    halfedge pointer in face corrupted." << std::endl;
            break;
        }
        // cycle-around-face test.
        Halfedge_const_handle h = get_face_halfedge( fbegin);
        if ( level >= 3 && h != Halfedge_const_handle()) {
            Halfedge_const_handle g = h;
            do {
                verr << "    halfedge " << n << std::endl;
                ++n;
                h = h->next();
                valid = valid && ( n <= hds->size_of_halfedges() && n!=0);
                if ( ! valid)
                    verr << "    too many halfedges around faces."
                         << std::endl;
            } while ( valid && (h != g));
        }
        ++f;
    }
    if ( valid && f != hds->size_of_faces())
        verr << "counting faces failed." << std::endl;
    if ( valid && level >= 3 && check_tag( Supports_face_halfedge()) &&
         n + nb  != hds->size_of_halfedges())
        verr << "counting halfedges via faces failed." << std::endl;
    valid = valid && ( f == hds->size_of_faces());
    valid = valid && ( level < 3 ||
                       ! check_tag( Supports_face_halfedge()) ||
                       n + nb  == hds->size_of_halfedges());

    if ( level >= 4) {
        verr << "level 4: normalized_border_is_valid( verbose = true)"
             << std::endl;
        valid = valid && ( normalized_border_is_valid( verb));
    }
    verr << "end of CGAL::HalfedgeDS_const_decorator<HDS>::"
            "is_valid(): structure is " << ( valid ? "valid." :
            "NOT VALID.") << std::endl;
    return valid;
}

} //namespace CGAL

#endif // CGAL_HALFEDGEDS_CONST_DECORATOR_H //
// EOF //
