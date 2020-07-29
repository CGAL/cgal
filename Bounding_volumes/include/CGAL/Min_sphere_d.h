// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner

#ifndef CGAL_MIN_SPHERE_D_H
#define CGAL_MIN_SPHERE_D_H

#include <CGAL/license/Bounding_volumes.h>


// Class declarations
// ==================

// Class interface and implementation
// ==================================
// includes

#  include <CGAL/basic.h>

#  include <CGAL/Optimisation/assertions.h>

#  include <CGAL/Optimisation/basic.h>

#  include <CGAL/Min_sphere_d/Optimisation_sphere_d.h>


#include <list>

#include <iostream>

namespace CGAL {


template <class Traits>
class Min_sphere_d
{


    public:
        typedef typename Traits::Rep_tag        Rep_tag;
        typedef typename Traits::RT             RT;
        typedef typename Traits::FT             FT;
        typedef typename Traits::Point_d        Point; // Point type

        typedef  typename Traits::Access_dimension_d
          Access_dimension_d;
        typedef  typename Traits::Access_coordinates_begin_d
          Access_coordinates_begin_d;
        typedef  typename Traits::Construct_point_d
          Construct_point_d;

        typedef typename std::list<Point>::const_iterator
                Point_iterator;
        typedef typename std::list<Point>::const_iterator
                Support_point_iterator;

    private:
        typedef typename std::list<Point>::iterator             It;


    private:
        int                                     d;            // ambient dim
        std::list<Point>                        points;       // keeps P = Q_n
        Traits                                  tco;          // traits object
        Optimisation_sphere_d<Rep_tag, FT, RT, Point,Traits>
                                                ms_basis; // keeps  miniball
        It                                      support_end;  // delimites S

public:
    Min_sphere_d ()
    : d(-1), tco( Traits()), ms_basis (tco),
      support_end(points.begin())
        {}


    Min_sphere_d (const Traits& traits)
    : d(-1), tco( traits), ms_basis (tco),
      support_end(points.begin())
        {}




    // STL-like constructor (member template)
    template <class InputIterator>
    Min_sphere_d( InputIterator first,
                       InputIterator last)
      : d(-1),
#if ( _MSC_VER != 1300)
      points( first, last),
#endif
      tco( Traits()),
      ms_basis (tco)
#if ( _MSC_VER != 1300)
      ,support_end(points.begin())
#endif
    {
#if ( _MSC_VER == 1300)
      std::copy(first,last,std::back_inserter(points));
      support_end = points.begin();
#endif
        if (points.size()>0) {
            d = tco.access_dimension_d_object() (points.front());
            CGAL_optimisation_precondition ((d>=0) && all_points_have_dim(d));
            ms_basis.get_sphere(Rep_tag()).set_size (d);
            pivot_ms();
        }
    }

    template <class InputIterator>
    Min_sphere_d( InputIterator first,
                       InputIterator last,
                       const Traits& traits)
      : d(-1),
#if ( _MSC_VER != 1300)
      points( first, last),
#endif
      tco( traits),
      ms_basis (tco)
#if ( _MSC_VER != 1300)
      ,support_end(points.begin())
#endif

    {
#if ( _MSC_VER == 1300)
      std::copy(first,last,std::back_inserter(points));
      support_end = points.begin();
#endif
        if (points.size()>0) {
            d = tco.access_dimension_d_object() (points.front());
            CGAL_optimisation_precondition ((d>=0) && all_points_have_dim(d));
            ms_basis.get_sphere(Rep_tag()).set_size (d);
            pivot_ms();
        }
    }

    Min_sphere_d (const Min_sphere_d& msph)
    : d(msph.ambient_dimension()),
      points (msph.points_begin(), msph.points_end()), tco (msph.tco),
          ms_basis (tco), support_end (points.begin())

    {
        if (d != -1) {
            ms_basis.get_sphere(Rep_tag()).set_size (d);
            pivot_ms();
        }
    }

    Min_sphere_d& operator=(const Min_sphere_d& msph)
    {
        if (this != &msph) {
            points.erase (points.begin(), points.end());
            d = msph.ambient_dimension();
            points.insert
              (points.begin(), msph.points_begin(), msph.points_end());
              ms_basis.get_sphere(Rep_tag()).set_tco(msph.tco);
            support_end = points.begin();
            tco = msph.tco;
            if (d != -1) {
                ms_basis.get_sphere(Rep_tag()).set_size (d);
                pivot_ms();
            }
        }
        return *this;
    }


    int number_of_points() const
    {
        return static_cast<int>(points.size());
    }

    int number_of_support_points() const
    {
        return ms_basis.get_sphere(Rep_tag()).number_of_support_points();
    }

    Point_iterator points_begin () const
    {
        return points.begin();
    }

    Point_iterator points_end () const
    {
        return points.end();
    }

    Support_point_iterator support_points_begin () const
    {
        return points.begin();
    }

    Support_point_iterator support_points_end () const
    {
        return support_end;
    }

    int ambient_dimension () const
    {
        return d;
    }

    Point center () const
    {
        CGAL_optimisation_precondition (!is_empty());
        return ms_basis.get_sphere(Rep_tag()).center();
    }

    FT squared_radius () const
    {
        CGAL_optimisation_precondition (!is_empty());
        return ms_basis.get_sphere(Rep_tag()).squared_radius();
    }


    Bounded_side bounded_side (const Point& p) const
    {
        if (d == -1)
           return ON_UNBOUNDED_SIDE;
        else {
          CGAL_optimisation_precondition
           (d == tco.access_dimension_d_object()(p));
           return (Bounded_side
               (-CGAL::sign (ms_basis.get_sphere(Rep_tag()).excess (p))));
        }
    }

    bool has_on_bounded_side (const Point& p) const
    {
        if (d == -1)
           return false;
        else {
          CGAL_optimisation_precondition
           (d == tco.access_dimension_d_object()(p));
           return (CGAL_NTS is_negative (ms_basis.get_sphere(Rep_tag()).excess (p)));
        }
    }

    bool has_on_unbounded_side (const Point& p) const
    {
        if (d == -1)
           return true;
        else {
          CGAL_optimisation_precondition
          (d == tco.access_dimension_d_object()(p));
           return (CGAL_NTS is_positive (ms_basis.get_sphere(Rep_tag()).excess (p)));
        }
    }

    bool has_on_boundary (const Point& p) const
    {
        if (d == -1)
           return false;
        else {
          CGAL_optimisation_precondition
          (d == tco.access_dimension_d_object()(p));
           return (CGAL_NTS is_zero (ms_basis.get_sphere(Rep_tag()).excess (p)));
        }
    }

    bool is_empty () const
    {
        return (d == -1);
    }

    bool is_degenerate () const
    {
        return (ms_basis.get_sphere(Rep_tag()).number_of_support_points() < 2);
    }


    void clear ()
    {
         d = -1;
         points.erase (points.begin(), points.end());
         ms_basis.get_sphere(Rep_tag()).set_size (-1);
         support_end = points.begin();
    }

    // STL-like set(member template)
    template <class InputIterator>
    void set ( InputIterator first,
               InputIterator last)
    {
        points.erase (points.begin(), points.end());
        points.insert (points.begin(), first, last);
        support_end = points.begin();
        if (points.size()>0) {
            d = tco.access_dimension_d_object() (points.front());
            CGAL_optimisation_precondition ((d>=0) && all_points_have_dim (d));
            ms_basis.get_sphere(Rep_tag()).set_size (d);
            pivot_ms();
        } else {
            d = -1;
            ms_basis.get_sphere(Rep_tag()).set_size (-1);
        }
    }

    void insert (const Point& p)
    {
        if (has_on_unbounded_side (p)) {
            if (is_empty()) {
                d = tco.access_dimension_d_object() (p);
                CGAL_optimisation_precondition (d>=0);
                ms_basis.get_sphere(Rep_tag()).set_size (d);
            }
            // ensure precondition of pivot_ms
            ms_basis.get_sphere(Rep_tag()).push (p);
            pivot_ms ();
            ms_basis.get_sphere(Rep_tag()).pop ();
            points.push_front (p);  // ensure postcondition of insert
        } else
            points.push_back (p);   // just append p
        if (support_end == points.end()) --support_end;
    }

    template <class InputIterator>
    void insert (InputIterator first, InputIterator last)
    {
        for (InputIterator i=first; i!=last; ++i)
            insert (*i);
    }


  bool is_valid (bool verbose = false, int /* level */ = 0) const
    {
        Verbose_ostream verr (verbose);

        // sphere verification
        verr << "  (a) sphere verification..." << std::flush;
        if (ms_basis.get_sphere(Rep_tag()).is_valid (verbose))
            verr << "passed." << std::endl;
        else
            return false;

        // containment check
        verr << "  (b) containment check..." << std::flush;

        // non-support-points
        Point_iterator i;
        for (i=support_end; i!=points.end(); ++i)
            if (has_on_unbounded_side (*i))
                return (_optimisation_is_valid_fail (verr,
                  "sphere does not contain all points"));

        // support points
        for (i=points.begin(); i!=support_end; ++i)
            if (!has_on_boundary (*i))
                return (_optimisation_is_valid_fail (verr,
                  "sphere does not have all support points on boundary"));

        verr << "passed." << std::endl;
        verr << "object is valid!" << std::endl;
        return true;
    }


    const Traits& traits() const
    {
         return tco;
    }


private:
    void mtf_ms (It k)
    {
        support_end = points.begin();
        if (ms_basis.get_sphere(Rep_tag()).size_of_basis()==d+1) return;
        for (It i = points.begin(); i!=k;) {
            It j = i++;
            if (CGAL_NTS is_positive (ms_basis.get_sphere(Rep_tag()).excess(*j))) {
                ms_basis.get_sphere(Rep_tag()).push (*j);
                mtf_ms (j);
                ms_basis.get_sphere(Rep_tag()).pop();
                move_to_front (j);
            }
        }
    }


    void pivot_ms ()
    {
        It t = points.begin();
        std::advance (t, (std::min) (d+1, (int)points.size()));
        mtf_ms (t);

        RT excess, e;
        do {
            excess = RT(0);
            It pivot;
            for (It i=t; i!=points.end(); ++i) {
                e = ms_basis.get_sphere(Rep_tag()).excess(*i);
                if (e > excess) {
                   excess = e;
                   pivot = i;
                }
            }
            if (CGAL_NTS is_positive (excess)) {
                t = support_end;
                if (t==pivot) ++t; //  inserted from the esa code
                ms_basis.get_sphere(Rep_tag()).push (*pivot);
                mtf_ms (support_end);
                ms_basis.get_sphere(Rep_tag()).pop();
                move_to_front (pivot);
            }
        } while (CGAL_NTS is_positive (excess));
    }


    void move_to_front (It j)
    {
        if (support_end == j)
           ++support_end;
        points.splice (points.begin(), points, j);
    }


    bool all_points_have_dim (int dim) const
    {
        for (Point_iterator i=points.begin(); i!=points.end(); ++i)
            if (tco.access_dimension_d_object()(*i) != dim)
                return false;
        return true;
    }


};

// Function declarations
// =====================
// I/O
// ---

template < class Traits >
std::ostream&
operator << ( std::ostream& os, const Min_sphere_d<Traits>& ms);

template < class Traits >
std::istream&
operator >> ( std::istream& is, Min_sphere_d<Traits> & ms);

} //namespace CGAL

#include <CGAL/Min_sphere_d/Min_sphere_d_impl.h>

#endif // CGAL_MIN_SPHERE_D_H

// ===== EOF ==================================================================
