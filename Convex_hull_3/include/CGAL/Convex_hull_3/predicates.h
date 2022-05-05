// Copyright (c) 2022 INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Samuel Hornus and Sébastien Loriot
//

#ifndef CGAL_CONVEX_HULL_3_PREDICATES_H
#define CGAL_CONVEX_HULL_3_PREDICATES_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <vector>

namespace CGAL {

namespace Convex_hull_3 {

namespace predicates_impl
{
// vec.h
//----------------------------

/*
 *
 *  Vec3
 *
 */

template< typename NT >
struct Vec3
{
    typedef Vec3<NT>    Self;
    typedef NT          NumberType;

#define vecset(x,y,z) data_[0]=(x); data_[1]=(y); data_[2]=(z);
    Vec3() { vecset(NT(0), NT(0), NT(0)) }
    Vec3(const Self & s) = default;
    Vec3(const NT & x, const NT & y, const NT & z) { vecset(x,y,z) }
    inline void set(const NT & x, const NT & y, const NT & z) {
      vecset(x,y,z)
    }
#undef vecset

    inline NT operator()(const int i) const { assert(3 > i && 0 <= i); return data_[i]; }
    inline NT operator[](const int i) const { return data_[i]; }
    inline NT & operator[](const int i) { return data_[i]; }
    inline NT x() const { return data_[0]; }
    inline NT y() const { return data_[1]; }
    inline NT z() const { return data_[2]; }
    inline NT & operator()(const int i) { assert(3 > i && 0 <= i); return data_[i]; }
    inline NT & x() { return data_[0]; }
    inline NT & y() { return data_[1]; }
    inline NT & z() { return data_[2]; }
    const NT * data() const { return data_; }

    inline NT squaredLength() const { return x()*x()+y()*y()+z()*z(); }
    inline NT l1norm() const { return std::max(std::abs(x()), std::max(std::abs(y()), std::abs(z()))); }
    inline NT norm2() const { return x()*x()+y()*y()+z()*z(); }
    inline NT length() const { return ::sqrt(squaredLength()); }
    inline void normalize() { const NT l = length(); x() /= l; y() /= l; z() /= l; }
    inline Self normalized() const { const NT l = length(); return Self(x()/l, y()/l, z()/l); }
    inline void LInfNormalize() { const NT l = l1norm(); x() /= l; y() /= l; z() /= l; }
    inline Self LInfNormalized() const { const NT l = l1norm(); return Self(x()/l, y()/l, z()/l); }
    inline static Self cross(const Self & lhs, const Self & rhs) {
      return Self(lhs.y() * rhs.z() - lhs.z() * rhs.y(),
          lhs.z() * rhs.x() - lhs.x() * rhs.z(),
          lhs.x() * rhs.y() - lhs.y() * rhs.x());
    }

    inline Self operator-(const Self & rhs) const { return Self(x()-rhs.x(), y()-rhs.y(), z()-rhs.z()); }
    inline Self operator-() const { return Self(-data_[0], data_[1], -data_[2]); }
    inline Self operator+(const Self & rhs) const { return Self(x()+rhs.x(), y()+rhs.y(), z()+rhs.z()); }
    inline NT operator|(const Self & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z(); }
    inline void negate() { data_[0] = - data_[0]; data_[1] = - data_[1]; data_[2] = - data_[2]; }
    inline Self negated() const { return Self(-x(), -y(), -z()); }

    NT data_[3];
};

template< typename NT >
Vec3<NT> operator/(const Vec3<NT> & lhs, const NT & rhs)
{
    return Vec3<NT>(lhs.x()/rhs, lhs.y()/rhs, lhs.z()/rhs);
}

template< typename NT >
Vec3<NT> operator*(const NT & lhs, const Vec3<NT> & rhs)
{
    return Vec3<NT>(lhs*rhs.x(), lhs*rhs.y(), lhs*rhs.z());
}

template< typename NT >
bool operator==(const Vec3<NT> & lhs, const Vec3<NT> & rhs)
{
    return ( lhs.x() == rhs.x() ) && ( lhs.y() == rhs.y() ) && ( lhs.z() == rhs.z() );
}

/*
 *
 *  Vec4
 *
 */

template< typename NT >
struct Vec4
{
    typedef Vec4<NT>    Self;
    typedef Vec3<NT>    V3;
    typedef NT          NumberType;

#define vecset(x,y,z,w) data_[0]=(x); data_[1]=(y); data_[2]=(z); data_[3] = w;
    Vec4() { vecset(NT(0), NT(0), NT(0), NT(0)) }
    Vec4(const Self & s) { vecset(s.x(), s.y(), s.z(), s.w()) }
    Vec4(const V3 & s) { vecset(s.x(), s.y(), s.z(), NT(1)) }
    Vec4(const NT & x, const NT & y, const NT & z, const NT & w) { vecset(x,y,z,w) }
    inline void set(const NT & x, const NT & y, const NT & z, const NT & w) {
      vecset(x,y,z,w)
    }
#undef vecset

    inline NT operator()(const int i) const { assert(4 > i && 0 <= i); return data_[i]; }
    inline NT operator[](const int i) const { return data_[i]; }
    inline NT & operator[](const int i) { return data_[i]; }
    inline NT x() const { return data_[0]; }
    inline NT y() const { return data_[1]; }
    inline NT z() const { return data_[2]; }
    inline NT w() const { return data_[3]; }
    inline NT & operator()(const int i) { assert(4 > i && 0 <= i); return data_[i]; }
    inline NT & x() { return data_[0]; }
    inline NT & y() { return data_[1]; }
    inline NT & z() { return data_[2]; }
    inline NT & w() { return data_[3]; }
    const NT * data() const { return data_; }

    V3 toVec3() const { return V3(x(), y(), z()); }

    inline NT squaredLength() const { return x()*x()+y()*y()+z()*z()+w()*w(); }
    inline NT length() const { return ::sqrt(squaredLength()); }
    inline void normalize() { NT l = length(); x() /= l; y() /= l; z() /= l; w() /= l; }

    inline Self operator-(const Self & rhs) const { return Self(x()-rhs.x(), y()-rhs.y(), z()-rhs.z(), w()-rhs.w()); }
    inline Self operator+(const Self & rhs) const { return Self(x()+rhs.x(), y()+rhs.y(), z()+rhs.z(), w()+rhs.w()); }
    inline NT operator|(const Self & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z() + w()*rhs.w(); }
    inline NT operator|(const V3 & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z(); }
    inline NT dotAs3(const Self & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z(); }
    //inline NT operator|(const V3 & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z() + w(); }
    inline NT operator()(const V3 & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z() + w(); }
    inline void negate() { data_[0] = - data_[0]; data_[1] = - data_[1]; data_[2] = - data_[2]; data_[3] = - data_[3]; }

    inline static V3 crossAs3(const Self & lhs, const Self & rhs) {
          return V3(lhs.y() * rhs.z() - lhs.z() * rhs.y(),
              lhs.z() * rhs.x() - lhs.x() * rhs.z(),
              lhs.x() * rhs.y() - lhs.y() * rhs.x());
        }

    NT data_[4];
};

template< typename NT >
Vec4<NT> operator*(const NT & lhs, const Vec4<NT> & rhs)
{
    return Vec4<NT>(lhs*rhs.x(), lhs*rhs.y(), lhs*rhs.z(), lhs*rhs.w());
}
//----------------------------


// spherical.h
//----------------------------

template <class NT>
struct SphericalPolygonElement {
  Vec3<NT> vertex_; // A vertex of the spherical polygon
  Vec3<NT> north_; // The north pole of the equatorial arc/edge leading OUT OF that vertex_ (arcs are oriented west-to-east, or CCW in a right-handed frame.
  // In the spherical polygon (v0, n0), (v1, n1), (v2, n2), ... we have
  // v1 = cross(n0, n1),  more generally: v_{i+1} = cross(n_i, n_{i+1})  and
  // n1 = cross(v1, v2),  more generally: n_i     = cross(v_i, v_{i+1}).
  SphericalPolygonElement(){}
  SphericalPolygonElement(const Vec3<NT> & n) : north_(n.normalized()) {}
  SphericalPolygonElement(const Vec3<NT> & v, const Vec3<NT> & n) : vertex_(v), north_(n) {}
};

template <class NT>
struct SphericalPolygon : public std::vector<SphericalPolygonElement<NT>> {

  typedef std::vector<SphericalPolygonElement<NT>> Base;
  typedef typename Base::iterator iterator;
  typedef typename Base::const_iterator const_iterator;
  SphericalPolygon() {
    this->reserve(16);
  }

  Vec3<NT> averageDirection() const {
    // PRECONDITION : all northes are normalized.
    switch( this->size() ) {
      case 0 : return Vec3<NT>(0.0f, 0.0f, 0.0f); break;
      case 1 : return this->begin()->north_; break;
      case 2 : return (*this)[0].north_ + (*this)[1].north_; break;
      default : {
                  Vec3<NT> avg;
                  for( const SphericalPolygonElement<NT> & v : *this )
                    avg = avg + v.vertex_;
                  return avg;
                } break;
    }
  }

  void clip(const Vec3<NT> & OrigVertex, SphericalPolygon<NT> & result, bool doClean=true) const {
    // PRECONDITION : clipNorth, and all northes are normalized.
#define _ray_spherical_eps 1e-6f
    const int n = this->size();
    result.clear();
    switch( n ) {
      case 0 : break; // 0 means empty, so nothing to do
      case 1 : {
                 result = (*this);
                 Vec3<NT> clipNorth = OrigVertex.normalized();
                 NT dot = this->begin()->north_ | clipNorth;
                 if( dot < -0.99984769515 ) { // about one degree
                   // intersection of two almost opposite hemispheres ==> empty
                   result.clear();
                   break;
                 } else if( dot > 0.99984769515 ) {
                   break;
                 }
                 Vec3<NT> v(Vec3<NT>::cross(clipNorth, this->begin()->north_).LInfNormalized());
                 result.begin()->vertex_ = v;
                 result.emplace_back(v.negated(), clipNorth);
                 break;
               }
      case 2 : {
                 result = (*this);
                 Vec3<NT> clipNorth = OrigVertex.LInfNormalized();
                 iterator next = result.begin();
                 iterator cur = next++;
                 NT vDot = this->begin()->vertex_ | clipNorth;
                 if( vDot >= _ray_spherical_eps ) {
                   // we'll get a triangle
                   next->vertex_ = Vec3<NT>::cross(clipNorth, next->north_).LInfNormalized();
                   Vec3<NT> v(Vec3<NT>::cross(cur->north_, clipNorth).LInfNormalized());
                   result.emplace(next, v, clipNorth);
                 } else if( vDot <= - _ray_spherical_eps ) {
                   // we'll get a triangle
                   cur->vertex_ = Vec3<NT>::cross(clipNorth, cur->north_).LInfNormalized();
                   Vec3<NT> v(Vec3<NT>::cross(next->north_, clipNorth).LInfNormalized());
                   result.emplace_back(v, clipNorth);
                 } else {
                   // we keep a moon crescent
                   NT curTest(clipNorth | Vec3<NT>::cross(cur->north_, cur->vertex_));
                   Vec3<NT> nextTest(Vec3<NT>::cross(next->north_, next->vertex_));
                   if( curTest > 0.0f ) {
                     if( (clipNorth | nextTest) <= 0.0f ) {
                       next->north_ = clipNorth;
                       cur->vertex_ = Vec3<NT>::cross(next->north_, cur->north_);
                       cur->vertex_.LInfNormalize();
                       next->vertex_ = cur->vertex_;
                       next->vertex_.negate();
                     } else {
                       // the crescent is unchanged
                       //std::cerr << "kept a crescent\n";
                     }
                   } else {
                     if( (clipNorth | nextTest) > 0.0f ) {
                       cur->north_ = clipNorth;
                       next->vertex_ = Vec3<NT>::cross(cur->north_, next->north_);
                       next->vertex_.LInfNormalize();
                       cur->vertex_ = next->vertex_;
                       cur->vertex_.negate();
                     } else {
                       //std::cerr << "killed a crescent\n";
                       result.clear();
                     }
                   }
                 }
                 break;
               }
      default : { // n >= 3
                  int nbKept(0);
                  const_iterator cur = this->begin();
                  Vec3<NT> clipNorth = OrigVertex.LInfNormalized();
                  NT nextDot, curDot = clipNorth | cur->vertex_;
                  while( cur != this->end() ) {
                    if( cur+1 == this->end() )
                      nextDot = clipNorth | this->begin()->vertex_;
                    else
                      nextDot = clipNorth | (cur+1)->vertex_;
                    if( curDot >= _ray_spherical_eps ) { // cur is "IN"
                      ++nbKept;
                      result.push_back(*cur);
                      if( nextDot <= -_ray_spherical_eps ) { // next is "OUT"
                        result.emplace_back(Vec3<NT>::cross(cur->north_, clipNorth).LInfNormalized(), clipNorth);
                      }
                    } else if( curDot > -_ray_spherical_eps ) { // cur is "ON" the clipping plane
                      ++nbKept;
                      if ( nextDot <= -_ray_spherical_eps ) // next is "OUT"
                        result.emplace_back(cur->vertex_, clipNorth);
                      else
                        result.push_back(*cur);
                    } else { // cur is "OUT"
                      if ( nextDot >= _ray_spherical_eps ) { // next is "IN"
                        result.emplace_back(Vec3<NT>::cross(clipNorth, cur->north_).LInfNormalized(), cur->north_);
                      }
                    }
                    curDot = nextDot;
                    ++cur;
                  }
                  if( (result.size() < 3/*too small*/) || ((nbKept == n)/*no change*/ && doClean) ) {
                    result.clear();
                  }
                  //if( nbKept == n ) {
                    //std::cerr << "**";
                  //}
                  break;
                }
    }
  }
};

//----------------------------

template <class NT, class PointRange>
Vec3<NT> vertex(const PointRange& p, int i)
{
  return Vec3<NT>(p[i].x(), p[i].y(), p[i].z());
}

template <class NT, class Convex>
bool differenceCoversZeroInDir(const Convex& A, const Convex& B, int & vA, int & vB, const Vec3<NT> & dir) {
  // difference above is: A - B
  NT maxOverA = vertex<NT>(A, 0) | dir;
  NT minOverB = vertex<NT>(B, 0) | dir;
  vA = vB = 0;
  const int na = A.size();
  for( int i = 1; i < na; ++i ) {
    NT tempA = vertex<NT>(A, i) | dir;
    if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
  }
  const int nb = B.size();
  for( int i = 1; i < nb; ++i ) {
    NT tempB = vertex<NT>(B, i) | dir;
    if( tempB < minOverB ) { minOverB = tempB; vB = i; }
  }
  return maxOverA >= minOverB;
}

template<typename NT, typename Convex >
bool sphericalDisjoint(const Convex & a, const Convex & b, unsigned long INTER_MAX_ITER) {
  SphericalPolygon<NT> positiveBound, tempPoly;
  int vA, vB;
  Vec3<NT> dir(vertex<NT>(b, 0) - vertex<NT>(a, 0));
  if( ! differenceCoversZeroInDir(a, b, vA, vB, dir) ) return true;
  positiveBound.clear();
  positiveBound.emplace_back(dir);
  positiveBound.clip(vertex<NT>(b, vB) - vertex<NT>(a, vA), tempPoly); positiveBound.swap(tempPoly);
  if( positiveBound.empty() ) return false;
  unsigned long planeStatPerPair = 0;
  do {
    if( ! differenceCoversZeroInDir(a, b, vA, vB, positiveBound.averageDirection()) ) return true;
    if(INTER_MAX_ITER!=0 && (++planeStatPerPair >= INTER_MAX_ITER))
    {
      return false;
    }
    positiveBound.clip(vertex<NT>(b, vB) - vertex<NT>(a, vA), tempPoly); positiveBound.swap(tempPoly);
    if( positiveBound.empty() ) return false;
  } while( true );
}

} // end of predicates_impl namespace

/**
* \ingroup PkgConvexHull3Predicates
*
* indicates if the convex hulls of two point sets intersect or not.
*
* @tparam PointRange1 a model of the concept `RandomAccessContainer`
* whose value type is a point type from a CGAL %Kernel
* @tparam PointRange2 a model of the concept `RandomAccessContainer`
* whose value type is the point type from a CGAL %Kernel
* @tparam NamedParameters_1 a sequence of \ref bgl_namedparameters "Named Parameters"
* @tparam NamedParameters_2 a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param r1 first range of points whose convex hull is considered in the do-intersect test
* @param r2 first range of points whose convex hull is considered in the do-intersect test
* @param np_1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
* @param np_2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of the range `r1` (`r2)}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `np1` and `np2`}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{`np1` only}
*   \cgalParamNEnd
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{if not `0` (no limit), indicates the maximum number of iterations that the algorithm is allowed to do.
*                           If this value is not `0`, then an intersection might be reported even if the convex hulls does not intersect.
                            However, if the convex hulls are reported not to intersect, this is guaranteed.}
*     \cgalParamType{an positive integer convertible to `std::size_t`}
*     \cgalParamDefault{`0`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
*/
template <class PointRange1, class PointRange2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
bool do_intersect(const PointRange1& r1, const PointRange2& r2,
                  const NamedParameters1 np1 = parameters::default_values(),
                  const NamedParameters2 np2 = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef Point_set_processing_3_np_helper<PointRange1, NamedParameters1> NP_helper1;
  typedef typename NP_helper1::Const_point_map PointMap1;
  typedef typename NP_helper1::Geom_traits Geom_traits;

  typedef Point_set_processing_3_np_helper<PointRange2, NamedParameters2> NP_helper2;
  typedef typename NP_helper2::Const_point_map PointMap2;

  typedef typename boost::property_traits<PointMap1>::value_type Point_3;
  CGAL_static_assertion((std::is_same<Point_3, typename boost::property_traits<PointMap2>::value_type>::value));

  PointMap1 point_map1 = NP_helper1::get_const_point_map(r1, np1);
  PointMap2 point_map2 = NP_helper2::get_const_point_map(r2, np2);

  // TODO: avoid doing a copy
  std::vector<typename Geom_traits::Point_3> a, b;
  a.reserve(r1.size());
  b.reserve(r2.size());
  for (const auto& p : r1)
    a.push_back(get(point_map1, p));
  for (const auto& p : r2)
    b.push_back(get(point_map2, p));

  unsigned int max_nb_iterations = choose_parameter(get_parameter(np1, internal_np::number_of_iterations), 0);

  return !predicates_impl::sphericalDisjoint<typename Geom_traits::FT>(a, b, max_nb_iterations);
}

template <class PointRange1, class PointRange2,
          class NamedParameters1 = parameters::Default_named_parameters,
          class NamedParameters2 = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
FT
#else
typename Point_set_processing_3_np_helper<PointRange1, NamedParameters1>::Geom_traits::FT
#endif
separation_distance(const PointRange1& r1, const PointRange2& r2,
                    const NamedParameters1 np1 = parameters::default_values(),
                    const NamedParameters2 np2 = parameters::default_values());



// TODO: add OBB dedicated code in OBB package (do-intersect for sure and check for separation_distance).

}} // CGAL::Convex_hull_3 namespace

#endif // CGAL_CONVEX_HULL_3_PREDICATES_H
