// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
//
// Every use of CGAL requires a license. 
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation. 
//
// Commercial licenses
// - Please check the CGAL web site http://www.cgal.org/index2.html for 
//   availability.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// release       : CGAL-2.4
// release_date  : 2002, May 16
//
// chapter       : $CGAL_Chapter: Optimisation $
// file          : include/CGAL/Min_sphere_of_spheres_pair.h
// package       : Min_sphere_of_spheres_d (1.00)
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Kaspar Fischer
//
// coordinator   : ETH Zurich (Kaspar Fischer)
//
// implementation: dD Smallest Enclosing Sphere of Spheres
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

#ifndef CGAL_MIN_SPHERE_OF_SPHERES_PAIR_H
#define CGAL_MIN_SPHERE_OF_SPHERES_PAIR_H

#include <utility>

namespace CGAL {

  // another namespace to "hide" implementation details:
  namespace Min_sphere_of_spheres_impl {

    template<typename FT>
    class Pair : public std::pair<FT,FT> {
    private:
      typedef std::pair<FT,FT> Base;
      
    public: // construction:
      inline Pair() : Base() {}
    
      inline Pair(const FT& a,const FT& b) : Base(a,b) {}
      
      inline Pair(int i) : Base(i,0) {}
      
      inline Pair& operator=(const FT& x) {
	first  = x;
	second = 0;
	return *this;
      }
      
    public:  // arithmetic and comparision:
      inline Pair operator+(const Pair& a) const {
	return Pair(first+a.first,second+a.second);
      }
      
      inline Pair operator-(const Pair& a) const {
	return Pair(first-a.first,second-a.second);
      }
      
      inline Pair operator-(const FT& a) const {
	return Pair(first-a,second);
      }
      
      inline Pair operator*(const FT& a) const {
	return Pair(first*a,second*a);
      }
      
      inline Pair operator/(const FT& a) const {
	CGAL_assertion(a != 0);
	return Pair(first/a,second/a);
      }
      
      inline Pair& operator+=(const Pair& p) {
	first  += p.first;
	second += p.second;
	return *this;
      }
      
      inline Pair& operator-=(const Pair& p) {
	first  -= p.first;
	second -= p.second;
	return *this;
      }
      
      inline bool operator!=(const Pair& p) const {
	return first!=p.first || second!=p.second;
      }
    };
    
    template<typename FT>
    inline Pair<FT> operator+(const FT& a,const Pair<FT>& p) {
      return Pair<FT>(a+p.first,p.second);
    }

    template<typename FT>
    inline Pair<FT> operator-(const FT& a,const Pair<FT>& p) {
      return Pair<FT>(a-p.first,-p.second);
    }
      
    template<typename FT>
    inline bool isNeg(const FT& p,const FT&) {
      return p < 0;
    }
    
    template<typename FT>
    inline bool isNeg(const Pair<FT> p,const FT& d) {
      const bool aneg = p.first<0, bneg = p.second<0;
      
      if (aneg && bneg)
	return true;
      if (!aneg && !bneg)
	return false;
      
      // So what remains are the cases (i) a<0,b>=0 and (ii) a>=0,b<0:
      //   (i)  We need to test b*sqrt(d)<-a with b,-a>=0.
      //   (ii) We need to test a<(-b)*sqrt(d) with a,-b>=0.
      // Hence:
      const FT x = CGAL_NTS square(p.second)*d, y = CGAL_NTS square(p.first);
      return aneg? x<y : y<x;
    }
    
    template<typename FT>
    inline bool isZero(const Pair<FT> p,const FT& d) {
      if (d != 0)
        // check whether the sides of a=-b*sqrt(d) (*)
        // have different signs:
        if ((p.first>0) ^ (p.second<0))
          return false;

      // Here we have either:
      //   (i)   d=0, or
      //   (ii)  a>0,b<0,d!=0, or
      //   (iii) a<=0,b>=0,d!=0
      // Hence both sides of (*) are either positive or negative.
      return CGAL_NTS square(p.first) == CGAL_NTS square(p.second)*d;
    }
    
    template<typename FT>
    inline bool isNegOrZero(const FT& p,const FT& d) {
      return p <= 0;
    }
    
    template<typename FT>
    inline bool isNegOrZero(const Pair<FT> p,const FT& d) {
      return isNeg(p,d) || isZero(p,d);
    }

  } // namespace Min_sphere_of_spheres_impl
    
} // namespace CGAL

#endif // MIN_SPHERE_OF_SPHERES_PAIR_H
