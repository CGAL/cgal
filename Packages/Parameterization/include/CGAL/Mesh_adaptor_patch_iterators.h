// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_MESH_ADAPTOR_PATCH_ITERATORS_H
#define CGAL_MESH_ADAPTOR_PATCH_ITERATORS_H

#include <CGAL/Mesh_adaptor_patch_vertex.h>
#include <CGAL/parameterization_assertions.h>

#include <list>

CGAL_BEGIN_NAMESPACE


/// Class Mesh_adaptor_patch_vertex_list
/// Type of the list of all vertices of a Mesh_adaptor_patch_3<PatchableMeshAdaptor_3> mesh
template<class PatchableMeshAdaptor_3>
class Mesh_adaptor_patch_vertex_list 
    : public std::list< Mesh_adaptor_patch_vertex<PatchableMeshAdaptor_3> >
{
// Public types
public:

    /// Export template parameter type
    typedef PatchableMeshAdaptor_3           Adaptor;
};


/// Class Mesh_adaptor_patch_vertex_list_iterator
/// Same behavior as Mesh_adaptor_patch_vertex_list<PatchableMeshAdaptor_3>::iterator 
/// + conversion to Mesh_adaptor_patch_vertex_handle
template<class PatchableMeshAdaptor_3>
class Mesh_adaptor_patch_vertex_list_iterator 
    : public Mesh_adaptor_patch_vertex_list<PatchableMeshAdaptor_3>::iterator
{
// Private types
private:

    typedef typename Mesh_adaptor_patch_vertex_list<PatchableMeshAdaptor_3>::iterator    
                                            Base;
    typedef Mesh_adaptor_patch_vertex_list_iterator 
                                            Self;
    typedef Mesh_adaptor_patch_vertex<PatchableMeshAdaptor_3> 
                                            Vertex;

// Public types
public:

    /// Export template parameter type
    typedef PatchableMeshAdaptor_3           Adaptor;

// Public operations
public:

    Mesh_adaptor_patch_vertex_list_iterator() {}
    Mesh_adaptor_patch_vertex_list_iterator(const Base& toCopy) : Base(toCopy) {}

    Mesh_adaptor_patch_vertex_list_iterator(const Self& toCopy) : Base(toCopy) {}
    Self& operator=(const Self& toCopy) { Base::operator=(toCopy); return *this; }

    Self & operator++()     { Base::operator++(); return *this; }
    Self & operator--()     { Base::operator--(); return *this; }
    Self operator++(int)    { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int)    { Self tmp(*this); --(*this); return tmp; }

    bool operator==(const Self& it) const { return Base::operator==(it); }
    bool operator!=(const Self& it) const { return ! (*this == it); }

    /// Comparison to NULL pointer
    bool operator==(CGAL_NULL_TYPE ptr) const { 
        CGAL_parameterization_assertion(ptr == NULL);
        return Base::operator==( Base() ); 
    }
    bool operator!=(CGAL_NULL_TYPE ptr) const { return ! (*this == ptr); }

    /// Conversion to handle
    operator Mesh_adaptor_patch_vertex_handle<Adaptor>() const {
        return &*(*this);
    }
    operator Mesh_adaptor_patch_vertex_const_handle<Adaptor>() const {
        return &*(*this);
    }
};


/// Class Mesh_adaptor_patch_vertex_list_const_iterator
/// Same behavior as Mesh_adaptor_patch_vertex_list<PatchableMeshAdaptor_3>::const_iterator 
/// + conversion to Mesh_adaptor_patch_vertex_const_handle
template<class PatchableMeshAdaptor_3>
class Mesh_adaptor_patch_vertex_list_const_iterator 
    : public Mesh_adaptor_patch_vertex_list<PatchableMeshAdaptor_3>::const_iterator
{
// Private types
private:

    typedef typename Mesh_adaptor_patch_vertex_list<PatchableMeshAdaptor_3>::const_iterator    
                                            Base;
    typedef Mesh_adaptor_patch_vertex_list_const_iterator 
                                            Self;
    typedef Mesh_adaptor_patch_vertex<PatchableMeshAdaptor_3> 
                                            Vertex;

// Public types
public:

    /// Export template parameter type
    typedef PatchableMeshAdaptor_3           Adaptor;

// Public operations
public:

    Mesh_adaptor_patch_vertex_list_const_iterator() {}
    Mesh_adaptor_patch_vertex_list_const_iterator(const Base& toCopy) : Base(toCopy) {}

    Mesh_adaptor_patch_vertex_list_const_iterator(const Self& toCopy) : Base(toCopy) {}
    Self& operator=(const Self& toCopy) { Base::operator=(toCopy); return *this; }

    Self & operator++()     { Base::operator++(); return *this; }
    Self & operator--()     { Base::operator--(); return *this; }
    Self operator++(int)    { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int)    { Self tmp(*this); --(*this); return tmp; }

    bool operator==(const Self& it) const { return Base::operator==(it); }
    bool operator!=(const Self& it) const { return ! (*this == it); }

    /// Comparison to NULL pointer
    bool operator==(CGAL_NULL_TYPE ptr) const { 
        CGAL_parameterization_assertion(ptr == NULL);
        return Base::operator==( Base() ); 
    }
    bool operator!=(CGAL_NULL_TYPE ptr) const { return ! (*this == ptr); }

    /// Conversion to handle
    operator Mesh_adaptor_patch_vertex_const_handle<Adaptor>() const {
        return &*(*this);
    }
};


CGAL_END_NAMESPACE

#endif //CGAL_MESH_ADAPTOR_PATCH_ITERATORS_H

