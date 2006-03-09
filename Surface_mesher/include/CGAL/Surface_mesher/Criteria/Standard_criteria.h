// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Steve OUDOT


#ifndef CGAL__STANDARD_CRITERIA_H
#define CGAL__STANDARD_CRITERIA_H

#include <cmath>
#include <vector>


namespace CGAL {

  namespace Surface_mesher {


  template <class Criterion>
  class Standard_criteria {

  protected:
    typedef typename Criterion::Quality FT;
    typedef std::vector<Criterion*> Criteria;
    Criteria criteria;

  public:
    typedef typename Criterion::Facet Facet;
    typedef std::vector<FT> Quality;

    Standard_criteria() {}

    Standard_criteria (const Criteria& c) : criteria (c) {}

    void set_criteria(const Criteria& c)
    {
      criteria = c;
    }

    bool is_bad (const Facet& f, Quality& q ) {
      bool bad = false;
      int i = 0;
      q.resize(criteria.size());
      for (typename Criteria::iterator cit = criteria.begin(); cit !=
	     criteria.end(); ++cit)
	if ((*cit)->is_bad (f, q[i++]))
          bad = true;
      return bad;
    }
  };

  // abstract basic Criterion class
  template <class Tr>
  class Refine_criterion {
  public:
    typedef typename Tr::Facet Facet;
    typedef typename Tr::Geom_traits::FT Quality;
    virtual bool is_bad (const Facet&, Quality& ) const = 0;
    virtual ~Refine_criterion() {}
  };

  // Aspect_ratio Criterion class
  template <class Tr>
  class Aspect_ratio_criterion : public Refine_criterion <Tr> {
  public:
    typedef typename Refine_criterion <Tr>::Quality Quality;
    typedef typename Tr::Geom_traits::FT FT;

  private:
    typedef typename Refine_criterion <Tr>::Facet Facet;
    typedef typename Tr::Point Point;

    Quality B;

  public:
    // Nb: the default bound of the criterion is such that the criterion
    // is always fulfilled
    Aspect_ratio_criterion(const FT angle_min = 0.)
    { // TODO: document that FT must constructible from a double!
      set_angle_min(angle_min);
    }

    inline
    Quality bound() const { return B; }

    inline
    Quality angle_min() const { return std::asin (std::sqrt(B)); }

    inline
    void set_bound(const Quality b) { B = b; };

    inline
    void set_angle_min(const FT angle_min) {
      B = std::sin (CGAL_PI * CGAL::to_double(angle_min) / 180);
      B = B * B;
    };

    bool is_bad (const Facet& fh, Quality& q) const {
      CGAL_assertion (fh.first->is_facet_on_surface (fh.second));

      typedef typename Tr::Geom_traits Geom_traits;
      Geom_traits gt;
      typename Geom_traits::Triangle_3 t;

      Point p1 = fh.first->vertex ((fh.second+1)&3)->point();
      Point p2 = fh.first->vertex ((fh.second+2)&3)->point();
      Point p3 = fh.first->vertex ((fh.second+3)&3)->point();
      t = gt.construct_triangle_3_object()(p1,p2,p3);

      Quality d12,d13,d23;
      d12=gt.compute_squared_distance_3_object()(p1,p2);
      d13=gt.compute_squared_distance_3_object()(p1,p3);
      d23=gt.compute_squared_distance_3_object()(p2,p3);

      Quality aspect_ratio = 4 * gt.compute_squared_area_3_object()(t)
	* min_3(d12,d13,d23) / (d12*d13*d23);

      CGAL_assertion (aspect_ratio >= 0 && aspect_ratio <= 1);
      q = aspect_ratio;
      return aspect_ratio < B;
    }


  private:
    static 
    Quality min_3 (const Quality a, const Quality b, const Quality c) {
      if (a<=b && a<=c)
	return(a);

      else if (b<=c)
	return(b);

      else
	return(c);
    }
  };  // end Aspect_ratio_criterion

  // Curvature_adapted size Criterion class
  template <class Tr>
  class Curvature_size_criterion : public Refine_criterion <Tr> {
  public:
    typedef typename Refine_criterion <Tr>::Quality Quality;

  private:
    typedef typename Refine_criterion <Tr>::Facet Facet;
    typedef typename Tr::Point Point;

    Quality B;

  public:
    // Nb: the default bound of the criterion is such that the criterion
    // is always fulfilled
    Curvature_size_criterion(const Quality b = 1000) : B(b * b) {}

    inline
    Quality bound() const { return std::sqrt (B); }

    inline
    void set_bound(const Quality b) { B = b * b; };


    bool is_bad (const Facet& fh, Quality& q) const {
      CGAL_assertion (fh.first->is_facet_on_surface (fh.second));

      typedef typename Tr::Geom_traits Geom_traits;
      Geom_traits gt;

      Point p1 = fh.first->vertex ((fh.second+1)&3)->point();
      Point p2 = fh.first->vertex ((fh.second+2)&3)->point();
      Point p3 = fh.first->vertex ((fh.second+3)&3)->point();

      q = gt.compute_squared_distance_3_object()
	(gt.construct_circumcenter_3_object()(p1,p2,p3),
	 fh.first->get_facet_surface_center(fh.second));
      return q > B;
    }

  };  // end Curvature_size_criterion

  // Uniform size Criterion class
  template <class Tr>
  class Uniform_size_criterion : public Refine_criterion <Tr> {
  public:
    typedef typename Refine_criterion <Tr>::Quality Quality;

  private:
    typedef typename Refine_criterion <Tr>::Facet Facet;
    typedef typename Tr::Point Point;

    Quality B;

  public:
    // Nb: the default bound of the criterion is such that the criterion
    // is always fulfilled
    Uniform_size_criterion(const Quality b = 1000) : B(b * b) {}

    inline
    Quality bound() const { return std::sqrt (B); }

    inline
    void set_bound(const Quality b) { B = b * b; };


    bool is_bad (const Facet& fh, Quality& q) const {
      CGAL_assertion (fh.first->is_facet_on_surface (fh.second));

      typedef typename Tr::Geom_traits Geom_traits;
      Geom_traits gt;

      Point p1 = fh.first->vertex ((fh.second+1)&3)->point();

      q =  gt.compute_squared_distance_3_object()
	(p1, fh.first->get_facet_surface_center (fh.second));
      return q > B;
    }
  };  // end Uniform_size_criterion



  // Edge size Criterion class
  template <class Tr>
  class Edge_size_criterion : public Refine_criterion <Tr> {
  public:
    typedef typename Refine_criterion <Tr>::Quality Quality;

  private:
    typedef typename Refine_criterion <Tr>::Facet Facet;
    typedef typename Tr::Point Point;

    Quality B;

  public:
    // Nb: the default bound of the criterion is such that the criterion
    // is always fulfilled
    Edge_size_criterion(const Quality b = 1000) : B(b * b) {}

    inline
    Quality bound() const { return std::sqrt (B); }

    inline
    void set_bound(const Quality b) { B = b * b; };


    bool is_bad (const Facet& fh, Quality& q) const {
      typedef typename Tr::Geom_traits Geom_traits;
      typedef typename Geom_traits::FT FT;
      Geom_traits gt;

      Point p1 = fh.first->vertex ((fh.second+1)&3)->point();
      Point p2 = fh.first->vertex ((fh.second+2)&3)->point();
      Point p3 = fh.first->vertex ((fh.second+3)&3)->point();

      FT d12 = gt.compute_squared_distance_3_object() (p1, p2);
      FT d13 = gt.compute_squared_distance_3_object() (p1, p3);
      FT d23 = gt.compute_squared_distance_3_object() (p2, p3);

      if (d12 > d13) {
	if (d12 > d23)
	  q =  d12;
	else
	  q = d23;
      }
      else
	q = d13;
      return q > B;
    }

  };  // end Edge_size_criterion


  }  // namespace Surface_mesher

}  // namespace CGAL


#endif  // end CGAL__STANDARD_CRITERIA_H

