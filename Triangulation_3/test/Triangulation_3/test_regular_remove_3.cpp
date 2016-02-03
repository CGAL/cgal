// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Christophe Delage (Christophe.Delage@sophia.inria.fr)

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>

#include <iostream>
#include <fstream>
#include <cassert>
//#include <cstdlib>
#include <set>

#include <boost/random/linear_congruential.hpp>
#include <boost/cstdint.hpp>

#include <CGAL/_test_types.h>
//#include <CGAL/_test_cls_regular_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel FK;
typedef CGAL::Regular_triangulation_euclidean_traits_3<FK> bare_traits;

int degeneracy_counter = 0;

// This traits class counts the number of times a power_test has returned 0.
// This gives a rough idea of how degenerate a data set is.
struct traits : public bare_traits
{
    struct Power_test_3 : public bare_traits::Power_test_3
    {
        typedef bare_traits::Power_test_3 P3;
        Oriented_side operator() (const Weighted_point &p0,
                                  const Weighted_point &p) const
        {
            Oriented_side result = P3::operator()(p0, p);
            if (result == 0) ++degeneracy_counter;
            return result;
        }
        Oriented_side operator() (const Weighted_point &p0,
                                  const Weighted_point &p1,
                                  const Weighted_point &p) const
        {
            Oriented_side result = P3::operator()(p0, p1, p);
            if (result == 0) ++degeneracy_counter;
            return result;
        }
        Oriented_side operator() (const Weighted_point &p0,
                                  const Weighted_point &p1,
                                  const Weighted_point &p2,
                                  const Weighted_point &p) const
        {
            Oriented_side result = P3::operator()(p0, p1, p2, p);
            if (result == 0) ++degeneracy_counter;
            return result;
        }
        Oriented_side operator() (const Weighted_point &p0,
                                  const Weighted_point &p1,
                                  const Weighted_point &p2,
                                  const Weighted_point &p3,
                                  const Weighted_point &p) const
        {
            Oriented_side result = P3::operator()(p0, p1, p2, p3, p);
            if (result == 0) ++degeneracy_counter;
            return result;
        }
    };

    Power_test_3 power_test_3_object() const
    { return Power_test_3(); }
};


// Explicit instantiation of the whole class :
template class CGAL::Regular_triangulation_3<traits>;


typedef CGAL::Regular_triangulation_3<traits>                 Cls;

typedef traits::Bare_point Point;
typedef traits::Weighted_point Weighted_point;

typedef Cls::Vertex_handle                         Vertex_handle;

// We don't want to use compare_xyz because it thinks two weighted points
// located at the same place with different weights are identical.
struct less_xyzw 
{
    bool operator() (const Weighted_point &p,
                     const Weighted_point &q) const
    {
        if (p.x() < q.x()) return true;
        if (p.x() > q.x()) return false;
        if (p.y() < q.y()) return true;
        if (p.y() > q.y()) return false;
        if (p.z() < q.z()) return true;
        if (p.z() > q.z()) return false;
        if (p.weight() < q.weight()) return true;
        return false;
    }
};

typedef std::set<Weighted_point,less_xyzw> point_set;
typedef point_set::iterator set_iterator;

// Base class for a weighted point generator.
class point_iterator
{
    Weighted_point _wp;
    int _i;
public:
    point_iterator (int i = 0) : _i (i) {}

    void operator++ () { --_i; }
    const Weighted_point &operator* () const { return _wp; }
    bool operator== (const point_iterator &i) const { return _i == i._i; }
    bool operator!= (const point_iterator &i) const { return ! (*this == i); }
    
protected:
    void set_point (int x, int y, int z, int w)
    {
        _wp = Weighted_point (Point (x, y, z), w);
    }
};

static boost::minstd_rand randgen;

// point_iterator_x generates points randomly on a grid (thus making lots of
// degenerate cases), staying in dimension x.
struct point_iterator_0 : public point_iterator
{
    point_iterator_0 () {}
    point_iterator_0 (int i) : point_iterator(i + 1) { ++*this; }

    void operator++ ()
    {
        int w = randgen() % 10;
        set_point (0, 0, 0, w);
        point_iterator::operator++();
    }
};

struct point_iterator_1 : public point_iterator
{
    point_iterator_1 () {}
    point_iterator_1 (int i) : point_iterator(i + 1) { ++*this; }

    void operator++ ()
    {
        int x = randgen() % 10;
        int w = randgen() % 10;
        set_point (x, 0, 0, w);
        point_iterator::operator++();
    }
};

struct point_iterator_2 : public point_iterator
{
    point_iterator_2 () {}
    point_iterator_2 (int i) : point_iterator(i + 1) { ++*this; }

    void operator++ ()
    {
        int x = randgen() % 10;
        int y = randgen() % 10;
        int w = randgen() % 10;
        set_point (x, y, 0, w);
        point_iterator::operator++();
    }
};

struct point_iterator_3 : public point_iterator
{
    point_iterator_3 () {}
    point_iterator_3 (int i) : point_iterator(i + 1) { ++*this; }

    void operator++ ()
    {
        int x = randgen() % 10;
        int y = randgen() % 10;
        int z = randgen() % 10;
        int w = randgen() % 10;
        set_point (x, y, z, w);
        point_iterator::operator++();
    }
};

class point_reader 
{
    std::istream *in;
    Weighted_point wp;
    int nb;
public:
    point_reader () : in (0), wp(), nb(0) {}
    point_reader (std::istream &is) : in(&is)
    {
        if (*in >> nb) {
            ++nb;
            ++*this;
        } else
            nb = 0;
    }
    point_reader &operator++ ()
    {
        if (nb > 0) {
            --nb;
            if (nb > 0) *in >> wp;
        }
        return *this;
    }
    bool operator== (const point_reader& p) const { return nb == p.nb; }
    bool operator!= (const point_reader& p) const { return nb != p.nb; }
    const Weighted_point &operator*() const { return wp; }
};

// Inserts number points in the triangulation and the multiset, using PI as the
// point generator.
template < class PI >
void insert (Cls &T, point_set &points, int number)
{
    int i = 1;
    for (PI pi (number), pend; pi != pend; ++pi, ++i) {
        points.insert (*pi);
        T.insert (*pi);
        std::cout << "\r number of inserted points: " << i << std::flush;
    }
    std::cout << std::endl;
    assert(T.is_valid());
    std::cout << " number of vertices: " << T.number_of_vertices()
        << std::endl;
    
    std::cout << " number of degeneracies: " << degeneracy_counter << std::endl;
    degeneracy_counter = 0;
}

// Removes number points from T, and removes each one from points also.
// Checks that each point of T that is removed exists in points.
void remove (Cls &T, point_set &points, int number)
{
    for (int i = 1; !points.empty(); ++i)
    {
        number--;
        assert(T.number_of_vertices() != 0);
        Vertex_handle v = T.finite_vertices_begin();
        set_iterator pos = points.find (v->point());
        assert(pos != points.end());
        T.remove (v);
        points.erase (pos);
        std::cout << "\r number of removed points:  " << i << std::flush;
    }
    std::cout << std::endl;
    
    assert(number >= 0);
    assert(T.number_of_vertices() == 0);
    assert(T.is_valid());
    assert(points.empty());
    std::cout << " number of degeneracies: " << degeneracy_counter << std::endl;
    degeneracy_counter = 0;
}

// Adds p which is supposed to increase the dimension of T from dim to dim+1,
// then removes it.
void dim_jump (Cls &T, const Point &p, int dim)
{
    assert(T.dimension() == dim);

    Vertex_handle v = T.insert (Weighted_point (p, 0));
    assert(T.is_valid());
    assert(v != Vertex_handle());
    assert(T.dimension() == dim + 1);
    
    T.remove (v);
    assert(T.is_valid());
    assert(T.dimension() == dim);
}


bool test_case (std::istream &is)
{
    point_reader pi (is), pend;
    if (pi == pend) return false;

    point_set points;
    Cls T;
    int number = 0;
    
    do {
        ++number;
        points.insert (*pi);
        T.insert (*pi);
    } while (++pi != pend);
    assert(T.is_valid());
    
    for (int i = 0; !points.empty(); ++i) {
        assert(T.number_of_vertices() != 0);
        Vertex_handle v = T.finite_vertices_begin();
        set_iterator pos = points.find (v->point());
        assert(pos != points.end());
        T.remove (v);
        points.erase(pos);
    }
    assert(T.number_of_vertices() == 0);
    assert(points.empty());
    
    return true;
}

int main(int argc, char **argv)
{
    std::cout << " with CGAL::Regular_triangulation_euclidean_traits_3: "
            << std::endl;

    // Test various data sets that crashed the code at some point in the past.
    // File format is:
    // number of points of first data set
    // point
    // ...
    // number of points of second data set
    // point
    // ...
    {
        std::ifstream fin ("data/regular_remove_3");
        assert(fin);
        std:: cout << " test `data/regular_remove_3'" << std::endl;
        while (test_case (fin)) 
            // semicolon
            ;
    }
    
    // Hardcoded seeds so that the test-suite is deterministic.
    boost::int32_t seed0 = 42, seed1 = 43, seed2 = 42, seed3 = 42;

    // You can also pass seeds on the command line.
    if (argc > 1) std::sscanf (argv[1], "%d", &seed0);
    if (argc > 2) std::sscanf (argv[2], "%d", &seed1);
    if (argc > 3) std::sscanf (argv[3], "%d", &seed2);
    if (argc > 4) std::sscanf (argv[4], "%d", &seed3);
    
    Cls T;
    point_set points;

    degeneracy_counter = 0;

    std::cout << " test dimension 0" << std::endl;
    randgen.seed(seed0);

    insert<point_iterator_0> (T, points, 10);
    dim_jump (T, Point (0, 0, 1), 0);
    remove (T, points, 10);

    std::cout << " test dimension 1" << std::endl;
    randgen.seed(seed1);

    insert<point_iterator_1> (T, points, 20);
    dim_jump (T, Point (0, 0, 1), 1);
    remove (T, points, 20);

    std::cout << " test dimension 2" << std::endl;
    randgen.seed(seed2);

    insert<point_iterator_2> (T, points, 100);
    dim_jump (T, Point (0, 0, 1), 2);
    remove (T, points, 100);

    std::cout << " test dimension 3" << std::endl;
    randgen.seed(seed3);

    insert<point_iterator_3> (T, points, 500);
    assert(T.dimension() == 3);
    remove (T, points, 500);

    std::cout << " quit" << std::endl;

    return 0;
}

