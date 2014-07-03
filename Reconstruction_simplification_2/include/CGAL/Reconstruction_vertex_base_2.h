// Copyright (c) 2014  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Fernando de Goes, Pierre Alliez

#ifndef RECONSTRUCTION_VERTEX_BASE_2_H_
#define RECONSTRUCTION_VERTEX_BASE_2_H_

#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/sample.h>

/// The Reconstruction_vertex_base_2 class (corresponding to
/// Reconstruction_vertex_base_2 in prototype) is the default
/// vertex class of the Reconstruction_triangulation_2 class.
///
/// - Each vertex stores a CSample as well as the corresponding relocated point
///
namespace CGAL {
/// @param Kernel   Geometric traits class
/// @param Vb   Vertex base class, model of TriangulationVertexBase_2.
template < class Kernel, class Vb = Triangulation_vertex_base_2<Kernel> >
class Reconstruction_vertex_base_2 : public Vb
{
public:
    typedef Vb Base;
    typedef Sample<Kernel> Sample;
    typedef typename Kernel::Point_2 Point;
    typedef typename Base::Face_handle Face_handle;

    template < typename TDS2 >
    struct Rebind_TDS {
        typedef typename Base::template Rebind_TDS<TDS2>::Other Vb2;
        typedef Reconstruction_vertex_base_2<Kernel,Vb2> Other;
    };

private:
    int   m_id;
    bool  m_pinned;
    Sample* m_sample;
    Point m_relocated;

public:
    Reconstruction_vertex_base_2()
    : Base()
    {
        m_id = -1;
        m_pinned = false;
        m_sample = NULL;
    }

    Reconstruction_vertex_base_2(const Point & p)
    : Base(p)
    {
        m_id = -1;
        m_pinned = false;
        m_sample = NULL;
    }

    Reconstruction_vertex_base_2(Face_handle f)
    : Base(f)
    {
        m_id = -1;
        m_pinned = false;
        m_sample = NULL;
    }

    Reconstruction_vertex_base_2(const Point & p, Face_handle f)
    : Base(p, f)
    {
        m_id = -1;
        m_pinned = false;
        m_sample = NULL;
    }

    virtual ~Reconstruction_vertex_base_2() { }

    int  id() const { return m_id; }
    int& id() { return m_id; }

    bool  pinned() const { return m_pinned; }
    bool& pinned() { return m_pinned; }

    Sample* get_sample() const { return m_sample; }
    void set_sample(Sample* sample) { m_sample = sample; }

    const Point& relocated() const { return m_relocated; }
    Point& relocated() { return m_relocated; }
};
//---------------STRUCT LESS VERTEX_HANDLE---------------------
template <class T>
struct less_Vertex_handle
{
    bool operator() (const T& a, const T& b) const
    {
        return (a->id() < b->id());
    }
};

} //end namespace

#endif /* RECONSTRUCTION_VERTEX_BASE_2_H_ */
