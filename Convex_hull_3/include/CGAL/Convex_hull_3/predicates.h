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

namespace CGAL {

namespace Convex_hull_3 {

#include <CGAL/Named_function_parameters.h>

namespace predicates_impl
{
// vec.h
//----------------------------
template< typename NT >
struct Vec2
{
    typedef Vec2<NT>    Self;
    typedef NT          NumberType;

#define vecset(x,y) data_[0]=(x); data_[1]=(y);
    Vec2() { vecset(NT(0), NT(0)) }
    Vec2(const Self & s) { vecset(s.x(), s.y()) }
    explicit Vec2(const NT & x) { vecset(x,x) }
    Vec2(const NT & x, const NT & y) { vecset(x,y) }
    inline void set(const NT & x, const NT & y) {
      vecset(x,y)
    }
#undef vecset

    inline NT operator()(const int i) const { assert(2 > i && 0 <= i); return data_[i]; }
    inline NT operator[](const int i) const { return data_[i]; }
    inline NT x() const { return data_[0]; }
    inline NT y() const { return data_[1]; }
    inline NT & operator()(const int i) { assert(2 > i && 0 <= i); return data_[i]; }
    inline NT & x() { return data_[0]; }
    inline NT & y() { return data_[1]; }
    const NT * data() const { return data_; }

    inline NT squaredLength() const { return x()*x()+y()*y(); }

    inline Self operator-(const Self & rhs) const { return Self(x()-rhs.x(), y()-rhs.y()); }
    inline Self operator+(const Self & rhs) const { return Self(x()+rhs.x(), y()+rhs.y()); }
    inline NT operator|(const Self & rhs) const { return x()*rhs.x() + y()*rhs.y(); }
    inline void negate() { data_[0] = - data_[0]; data_[1] = - data_[1]; }
    inline void rotate() { NT temp(data_[0]); data_[0] = - data_[1]; data_[1] = temp; }

    NT data_[2];
};

template< typename NT >
Vec2<NT> operator*(const NT & lhs, const Vec2<NT> & rhs)
{
    return Vec2<NT>(lhs*rhs.x(), lhs*rhs.y());
}

template< typename NT >
Vec2<NT> vecmin(const Vec2<NT> & lhs, const Vec2<NT> & rhs)
{
    return Vec2<NT>(std::min(lhs.x(), rhs.x()), std::min(lhs.y(), rhs.y()));
}

template< typename NT >
Vec2<NT> vecmax(const Vec2<NT> & lhs, const Vec2<NT> & rhs)
{
    return Vec2<NT>(std::max(lhs.x(), rhs.x()), std::max(lhs.y(), rhs.y()));
}

template< typename NT >
bool operator==(const Vec2<NT> & lhs, const Vec2<NT> & rhs)
{
    return ( lhs.x() == rhs.x() ) && ( lhs.y() == rhs.y() );
}

template< typename NT >
bool operator!=(const Vec2<NT> & lhs, const Vec2<NT> & rhs)
{
    return ! ( lhs == rhs );
}

template< typename Out, typename NT >
Out & operator<<(Out & out, const Vec2<NT> & v)
{
    out << v.x() << ' ' << v.y();
    return out;
}

template< typename NT >
bool leftTurn(const Vec2<NT> & a, const Vec2<NT> & b, const Vec2<NT> & c)
{
    Vec2<NT> p(a.y()-b.y(), b.x()-a.x());
    return (p | (c-a)) >= NT(0);
}

extern const Vec2<float> vec2_zero;

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
    Vec3(const Self & s) { vecset(s.x(), s.y(), s.z()) }
    Vec3(const Vec2<NT> & v2, const NT & z) { vecset(v2.x(),v2.y(),z) }
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
    inline Self operator+(const Self & rhs) const { return Self(x()+rhs.x(), y()+rhs.y(), z()+rhs.z()); }
    inline NT operator|(const Self & rhs) const { return x()*rhs.x() + y()*rhs.y() + z()*rhs.z(); }
    inline void negate() { data_[0] = - data_[0]; data_[1] = - data_[1]; data_[2] = - data_[2]; }
    inline Self negated() const { return Self(-x(), -y(), -z()); }

    NT data_[3];
};

extern const Vec3<float> vec3_zero, unit_x, unit_y, unit_z;

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
Vec3<NT> vecmin(const Vec3<NT> & lhs, const Vec3<NT> & rhs)
{
    return Vec3<NT>(std::min(lhs.x(), rhs.x()), std::min(lhs.y(), rhs.y()), std::min(lhs.z(), rhs.z()));
}

template< typename NT >
Vec3<NT> vecmax(const Vec3<NT> & lhs, const Vec3<NT> & rhs)
{
    return Vec3<NT>(std::max(lhs.x(), rhs.x()), std::max(lhs.y(), rhs.y()), std::max(lhs.z(), rhs.z()));
}

template< typename Out, typename NT >
Out & operator<<(Out & out, const Vec3<NT> & v)
{
    out << v.x() << ' ' << v.y() << ' ' << v.z();
    return out;
}

template< typename NT >
bool operator==(const Vec3<NT> & lhs, const Vec3<NT> & rhs)
{
    return ( lhs.x() == rhs.x() ) && ( lhs.y() == rhs.y() ) && ( lhs.z() == rhs.z() );
}

template< typename NT >
Vec3<NT> orthonormalVector(const Vec3<NT> & v) {
  Vec3<NT> r = ( v.x() * v.x() + v.y() * v.y() > 1e-4 ) ?
    Vec3<NT>(v.y(), -v.x(), NT(0))
    :
    Vec3<NT>(NT(0), v.z(), -v.y());
  r.normalize();
  return r;
}

/*
 *
 *  Vec4
 *
 */

template< typename NT >
struct alignas(16) Vec4
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

    NT data_[4];
};

template< typename NT >
Vec4<NT> operator*(const NT & lhs, const Vec4<NT> & rhs)
{
    return Vec4<NT>(lhs*rhs.x(), lhs*rhs.y(), lhs*rhs.z(), lhs*rhs.w());
}

template< typename NT >
Vec4<NT> vecmin(const Vec4<NT> & lhs, const Vec4<NT> & rhs)
{
    return Vec4<NT>(std::min(lhs.x(), rhs.x()),
                    std::min(lhs.y(), rhs.y()),
                    std::min(lhs.z(), rhs.z()),
                    std::min(lhs.w(), rhs.w()));
}

template< typename NT >
Vec4<NT> vecmax(const Vec4<NT> & lhs, const Vec4<NT> & rhs)
{
    return Vec4<NT>(std::max(lhs.x(), rhs.x()),
                    std::max(lhs.y(), rhs.y()),
                    std::max(lhs.z(), rhs.z()),
                    std::max(lhs.w(), rhs.w()));
}

template< typename Out, typename NT >
Out & operator<<(Out & out, const Vec4<NT> & v)
{
    out << v.x() << ' ' << v.y() << ' ' << v.z() << ' ' << v.w();
    return out;
}

/*
 *
 *  Common vector types
 *
 */

typedef Vec2<double>         Vec2d;
typedef Vec2<float>          Vec2f;
typedef Vec2<int>            Vec2i;
typedef Vec2<unsigned int>   Vec2ui;
typedef Vec3<double>         Vec3d;
typedef Vec3<float>          Vec3f;
typedef Vec3<int>            Vec3i;
typedef Vec3<unsigned int>   Vec3ui;
typedef Vec4<double>         Vec4d;
typedef Vec4<float>          Vec4f;
typedef Vec4<int>            Vec4i;
typedef Vec4<unsigned int>   Vec4ui;
typedef Vec4<unsigned short> Vec4us;

//----------------------------


// spherical.h
//----------------------------
struct SphericalPolygonElement {
  Vec3f vertex_;
  Vec3f north_;
  Vec3f silVertex_;
  SphericalPolygonElement(){}
  SphericalPolygonElement(const Vec3f & sil)
    : north_(sil.normalized()), silVertex_(sil) {}
  //SphericalPolygonElement(const Vec3f & v, const Vec3f & n)
  //  : vertex_(v), north_(n), silVertex_(n)  {}
  SphericalPolygonElement(const Vec3f & v, const Vec3f & n, const Vec3f & sil)
    : vertex_(v), north_(n), silVertex_(sil) {}
};

struct SphericalPolygon : public std::vector<SphericalPolygonElement> {

  typedef std::vector<SphericalPolygonElement> Base;
  typedef Base::iterator iterator;
  typedef Base::const_iterator const_iterator;
  SphericalPolygon() {
    reserve(16);
  }

  Vec3f averageDirection() const {
    // PRECONDITION : all northes are normalized.
    switch( size() ) {
      case 0 : return Vec3f(0.0f, 0.0f, 0.0f); break;
      case 1 : return begin()->north_; break;
      case 2 : return (*this)[0].north_ + (*this)[1].north_; break;
      default : {
                  Vec3f avg;
                  for( const SphericalPolygonElement & v : *this )
                    avg = avg + v.vertex_;
                  return avg;
                } break;
    }
  }

  void set_to_triangle(const Vec3f pts[3]) {
    clear();
    emplace_back(Vec3f::cross(pts[0], pts[1]).LInfNormalized(), pts[1].LInfNormalized(), pts[1]);
    emplace_back(Vec3f::cross(pts[1], pts[2]).LInfNormalized(), pts[2].LInfNormalized(), pts[2]);
    emplace_back(Vec3f::cross(pts[2], pts[0]).LInfNormalized(), pts[0].LInfNormalized(), pts[0]);
  }

  void clip(const Vec3f & OrigVertex, const Vec3f & silVertex, SphericalPolygon & result, bool doClean=true) const {
    // PRECONDITION : clipNorth, and all northes are normalized.
#define _ray_spherical_eps 1e-6f
    const int n = size();
    result.clear();
    switch( n ) {
      case 0 : break;
      case 1 : {
                 result = (*this);
                 Vec3f clipNorth = silVertex.normalized();
                 float dot = begin()->north_ | clipNorth;
                 if( dot < -0.99984769515 ) { // about one degree
                   // intersection of two almost opposite hemispheres ==> empty
                   result.clear();
                   break;
                 } else if( dot > 0.99984769515 ) {
                   break;
                 }
                 Vec3f v(Vec3f::cross(clipNorth, begin()->north_).LInfNormalized());
                 result.begin()->vertex_ = v;
                 result.emplace_back(v.negated(), clipNorth, OrigVertex);
                 break;
               }
      case 2 : {
                 result = (*this);
                 Vec3f clipNorth = silVertex.LInfNormalized();
                 iterator next = result.begin();
                 iterator cur = next++;
                 float vDot = begin()->vertex_ | clipNorth;
                 if( vDot >= _ray_spherical_eps ) {
                   // we'll get a triangle
                   next->vertex_ = Vec3f::cross(clipNorth, next->north_).LInfNormalized();
                   Vec3f v(Vec3f::cross(cur->north_, clipNorth).LInfNormalized());
                   result.emplace(next, v, clipNorth, OrigVertex);
                 } else if( vDot <= - _ray_spherical_eps ) {
                   // we'll get a triangle
                   cur->vertex_ = Vec3f::cross(clipNorth, cur->north_).LInfNormalized();
                   Vec3f v(Vec3f::cross(next->north_, clipNorth).LInfNormalized());
                   result.emplace_back(v, clipNorth, OrigVertex);
                 } else {
                   // we keep a moon crescent
                   float curTest(clipNorth | Vec3f::cross(cur->north_, cur->vertex_));
                   Vec3f nextTest(Vec3f::cross(next->north_, next->vertex_));
                   if( curTest > 0.0f ) {
                     if( (clipNorth | nextTest) <= 0.0f ) {
                       next->north_ = clipNorth;
                       next->silVertex_ = OrigVertex;
                       cur->vertex_ = Vec3f::cross(next->north_, cur->north_);
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
                       cur->silVertex_ = OrigVertex;
                       next->vertex_ = Vec3f::cross(cur->north_, next->north_);
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
                  const_iterator cur = begin();
                  Vec3f clipNorth = silVertex.LInfNormalized();
                  float nextDot, curDot = clipNorth | cur->vertex_;
                  while( cur != end() ) {
                    if( cur+1 == end() )
                      nextDot = clipNorth | begin()->vertex_;
                    else
                      nextDot = clipNorth | (cur+1)->vertex_;
                    if( curDot >= _ray_spherical_eps ) { // cur is "IN"
                      ++nbKept;
                      result.push_back(*cur);
                      if( nextDot <= -_ray_spherical_eps ) { // next is "OUT"
                        result.emplace_back(Vec3f::cross(cur->north_, clipNorth).LInfNormalized(), clipNorth, OrigVertex);
                      }
                    } else if( curDot > -_ray_spherical_eps ) { // cur is "ON" the clipping plane
                      ++nbKept;
                      if ( nextDot <= -_ray_spherical_eps ) // next is "OUT"
                        result.emplace_back(cur->vertex_, clipNorth, OrigVertex);
                      else
                        result.push_back(*cur);
                    } else { // cur is "OUT"
                      if ( nextDot >= _ray_spherical_eps ) { // next is "IN"
                        result.emplace_back(Vec3f::cross(clipNorth, cur->north_).LInfNormalized(), cur->north_, cur->silVertex_);
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
  void clip(const Vec3f & silVertex, SphericalPolygon & result, bool doClean=true) const {
    clip(silVertex, silVertex, result, doClean);
  }
};

//----------------------------

template <class PointRange>
Vec3f vertex(const PointRange& p, int i)
{
  return Vec3f(p[i].x(), p[i].y(), p[i].z());
}

template <class Convex>
bool differenceCoversZeroInDir(const Convex& A, const Convex& B, int & vA, int & vB, const Vec3f & dir) {
  // difference above is: A - B
  float maxOverA = vertex(A, 0) | dir;
  float minOverB = vertex(B, 0) | dir;
  vA = vB = 0;
  const int na = A.size();
  for( int i = 1; i < na; ++i ) {
    float tempA = vertex(A, i) | dir;
    if( tempA > maxOverA ) { maxOverA = tempA; vA = i; }
  }
  const int nb = B.size();
  for( int i = 1; i < nb; ++i ) {
    float tempB = vertex(B, i) | dir;
    if( tempB < minOverB ) { minOverB = tempB; vB = i; }
  }
  return maxOverA >= minOverB;
}

template< typename Convex >
bool sphericalDisjoint(const Convex & a, const Convex & b, unsigned long INTER_MAX_ITER) {
  static SphericalPolygon positiveBound, tempPoly;
  int vA, vB;
  Vec3f dir(vertex(b, 0) - vertex(a, 0));
  if( ! differenceCoversZeroInDir(a, b, vA, vB, dir) ) return true;
  positiveBound.clear();
  positiveBound.emplace_back(dir);
  positiveBound.clip(vertex(b, vB) - vertex(a, vA), tempPoly); positiveBound.swap(tempPoly);
  if( positiveBound.empty() ) return false;
  unsigned long planeStatPerPair = 0;
  do {
    if( ! differenceCoversZeroInDir(a, b, vA, vB, positiveBound.averageDirection()) ) return true;
    if(INTER_MAX_ITER!=0 && (++planeStatPerPair >= INTER_MAX_ITER))
    {
      return false;
    }
    positiveBound.clip(vertex(b, vB) - vertex(a, vA), tempPoly); positiveBound.swap(tempPoly);
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

  return !predicates_impl::sphericalDisjoint(a, b, max_nb_iterations);
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
