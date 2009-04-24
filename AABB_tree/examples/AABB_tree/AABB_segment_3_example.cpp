// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#include <iostream>
#include <list>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Plane_3 Plane;
typedef K::Segment_3 Segment;

template <class Iterator>
class Segment_primitive
{
        // type
public:
        typedef K::Segment_3 Object; // object type
        typedef Iterator Id; // Id type

        // member data
private:
        Id m_it; // iterator
        Object m_object; // 3D segment

public:
        Segment_primitive(Id it)
                : m_it(it)
        {
                m_object = *it; // copy segmetn
        }
public:
        const Object& object() const { return m_object; }
        Object& object() { return m_object; }
        Id id() { return m_it; }
};

typedef std::list<typename Segment>::iterator Iterator;
typedef Segment_primitive<Iterator> Primitive;
typedef CGAL::AABB_traits<typename K, typename Primitive> Traits;
typedef CGAL::AABB_tree<typename Traits> Polyhedron_tree;

int main(void)
{
        Point a(1.0, 0.0, 0.0);
        Point b(0.0, 1.0, 0.0);
        Point c(0.0, 0.0, 1.0);
        Point d(0.0, 0.0, 0.0);

        std::list<Segment> segments;
        segments.push_back(Segment(a,b));
        segments.push_back(Segment(a,c));
        segments.push_back(Segment(c,d));

        Polyhedron_tree tree(segments.begin(),segments.end());
        Plane plane(a,b,d);
        std::cout << tree.number_of_intersections(plane) 
                << " intersections(s) with plane" << std::endl;
        return 0;
}