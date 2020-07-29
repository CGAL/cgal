// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/LindstromTurk_params.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

#include <CGAL/Cartesian/MatrixC33.h>

#include <limits>
#include <vector>

// This should be in
//
// Implementation of the collapsing cost and placement strategy from:
//
//  "Fast and Memory Efficient Polygonal Symplification"
//  Peter Lindstrom, Greg Turk

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

template<class TM_, class Profile_>
class LindstromTurkCore
{
public:
  typedef TM_                                                            TM;
  typedef Profile_                                                       Profile;

  typedef boost::graph_traits<TM>                                        GraphTraits;
  typedef typename GraphTraits::vertex_descriptor                        vertex_descriptor;
  typedef typename GraphTraits::halfedge_descriptor                      halfedge_descriptor;

  typedef internal::LindstromTurk_params                                 LT_params;

  typedef typename Profile::VertexPointMap                               Vertex_point_pmap;
  typedef typename Profile::VertexPointMap                               Vertex_point_map;
  typedef typename boost::property_traits<Vertex_point_pmap>::value_type Point;
  typedef typename boost::property_traits<Vertex_point_pmap>::reference  Point_reference;

  typedef typename Profile::Geom_traits                                  Geom_traits;
  typedef typename Geom_traits::FT                                       FT;
  typedef typename Geom_traits::Vector_3                                 Vector;

  typedef boost::optional<FT>                                            Optional_FT;
  typedef boost::optional<Point>                                         Optional_point;
  typedef boost::optional<Vector>                                        Optional_vector;

  typedef MatrixC33<Geom_traits>                                         Matrix;

  typedef typename Profile::Triangle                                     Triangle;
  typedef std::vector<vertex_descriptor>                                 vertex_descriptor_vector;
  typedef std::vector<halfedge_descriptor>                               halfedge_descriptor_vector;

  typedef typename Profile::Triangle_vector::const_iterator              const_triangle_iterator;
  typedef typename halfedge_descriptor_vector::const_iterator            const_border_edge_iterator;

public:
  LindstromTurkCore(const LT_params& aParams, const Profile& aProfile);

  Optional_point compute_placement();
  Optional_FT compute_cost(const Optional_point& p);

private :
  struct Triangle_data
  {
    Triangle_data(const Vector& aNormalV, const FT& aNormalL)
      : NormalV(aNormalV), NormalL(aNormalL)
    {}

    Vector NormalV;
    FT NormalL;
  };

  struct Boundary_data
  {
    Boundary_data(const Point& s_, const Point& t_, const Vector& v_, const Vector& n_)
      : s(s_), t(t_), v(v_), n(n_) {}

    const Point s, t;
    const Vector v, n;
  };

  typedef std::vector<Triangle_data>                                       Triangle_data_vector;
  typedef std::vector<Boundary_data>                                       Boundary_data_vector;

private :
  void extract_triangle_data();
  void extract_boundary_data();

  void add_boundary_preservation_constraints(const Boundary_data_vector& aBdry);
  void add_volume_preservation_constraints(const Triangle_data_vector& triangles);
  void add_boundary_and_volume_optimization_constraints(const Boundary_data_vector& aBdry,
                                                        const Triangle_data_vector& triangles);
  void add_shape_optimization_constraints(const vertex_descriptor_vector& link);

  FT compute_boundary_cost(const Vector& v, const Boundary_data_vector& aBdry);
  FT compute_volume_cost(const Vector& v, const Triangle_data_vector& triangles);
  FT compute_shape_cost(const Point& p, const vertex_descriptor_vector& link);

  Point_reference get_point(const vertex_descriptor v) const
  {
    return get(mProfile.vertex_point_map(), v);
  }

  const Geom_traits& geom_traits() const { return mProfile.geom_traits(); }
  const TM& surface() const { return mProfile.surface(); }

  static Vector point_cross_product(const Point& a, const Point& b)
  {
    return cross_product(a-ORIGIN, b-ORIGIN);
  }

  // This is the (uX)(Xu) product described in the Lindstrom-Turk paper
  static Matrix LT_product(const Vector& u)
  {
    FT a00 = (u.y()*u.y()) + (u.z()*u.z());
    FT a01 = -u.x()*u.y();
    FT a02 = -u.x()*u.z();

    FT a10 = a01;
    FT a11 = (u.x()*u.x()) + (u.z()*u.z());
    FT a12 = - u.y() * u.z();

    FT a20 = a02;
    FT a21 = a12;
    FT a22 = (u.x()*u.x()) + (u.y()*u.y());

    return Matrix(a00, a01, a02,
                  a10, a11, a12,
                  a20, a21, a22);
  }

  static FT big_value() { return static_cast<FT>((std::numeric_limits<double>::max)()); }

  static bool is_finite(const FT& n) { return CGAL_NTS is_finite(n); }
  static bool is_finite(const Point& p) { return is_finite(p.x())  && is_finite(p.y())  && is_finite(p.z()); }
  static bool is_finite(const Vector& v) { return is_finite(v.x())  && is_finite(v.y())  && is_finite(v.z()); }
  static bool is_finite(const Matrix& m) { return is_finite(m.r0()) && is_finite(m.r1()) && is_finite(m.r2()); }

  template<class T>
  static boost::optional<T> filter_infinity(const T& n) { return is_finite(n) ? boost::optional<T>(n) : boost::optional<T>(); }


private:
  const LT_params& mParams;
  const Profile& mProfile;

  void add_constraint_if_alpha_compatible(const Vector& Ai, const FT& bi);
  void add_constraint_from_gradient(const Matrix& H, const Vector& c);

private:
  Triangle_data_vector mTriangle_data;
  Boundary_data_vector mBdry_data;

  int mConstraints_n;
  Matrix mConstraints_A;
  Vector mConstraints_b;

  FT mSquared_cos_alpha;
  FT mSquared_sin_alpha;
};

template<class TM, class K>
LindstromTurkCore<TM,K>::
LindstromTurkCore(const LT_params& aParams, const Profile& aProfile)
  :
    mParams(aParams),
    mProfile(aProfile),
    mConstraints_n(0),
    mConstraints_A(NULL_MATRIX),
    mConstraints_b(NULL_VECTOR)
{
  double alpha = 1.0 * CGAL_PI / 180.0;

  FT cos_alpha = std::cos(alpha);
  FT sin_alpha = std::sin(alpha);

  mSquared_cos_alpha = cos_alpha * cos_alpha;
  mSquared_sin_alpha = sin_alpha * sin_alpha;

  extract_triangle_data();
  extract_boundary_data();
}

template<class TM, class K>
void
LindstromTurkCore<TM,K>::
extract_boundary_data()
{
  mBdry_data.reserve(mProfile.border_edges().size());
  for(halfedge_descriptor border_edge : mProfile.border_edges())
  {
    halfedge_descriptor face_edge = opposite(border_edge, surface());

    vertex_descriptor sv = source(face_edge, surface());
    vertex_descriptor tv = target(face_edge, surface());

    const Point_reference sp = get_point(sv);
    const Point_reference tp = get_point(tv);

    Vector v = tp - sp;
    Vector n = point_cross_product(tp,sp);

    CGAL_SMS_LT_TRACE(1, "  Boundary edge. S:" << xyz_to_string(sp) << " T:" << xyz_to_string(tp)
                          << " V:" << xyz_to_string(v) << " N:" << xyz_to_string(n));

    mBdry_data.push_back(Boundary_data(sp,tp,v,n));
  }
}

template<class TM, class K>
void
LindstromTurkCore<TM,K>::
extract_triangle_data()
{
  mTriangle_data.reserve(mProfile.triangles().size());

  for(const Triangle& tri : mProfile.triangles())
  {
    const Point_reference p0 = get_point(tri.v0);
    const Point_reference p1 = get_point(tri.v1);
    const Point_reference p2 = get_point(tri.v2);

    Vector v01 = p1 - p0;
    Vector v02 = p2 - p0;

    Vector lNormalV = cross_product(v01,v02);
    FT lNormalL = point_cross_product(p0,p1) * (p2 - ORIGIN);

    CGAL_SMS_LT_TRACE(1, "  Extracting triangle v" << tri.v0 << "->v" << tri.v1 << "->v" << tri.v2
                           << " N:" << xyz_to_string(lNormalV) << " L:" << n_to_string(lNormalL));

    mTriangle_data.push_back(Triangle_data(lNormalV,lNormalL));
  }
}

template<class TM, class K>
typename LindstromTurkCore<TM,K>::Optional_point
LindstromTurkCore<TM,K>::
compute_placement()
{
  Optional_point rPlacementP;
  Optional_vector lPlacementV;

  CGAL_SMS_LT_TRACE(0, "Computing LT placement for E" << mProfile.v0_v1()
                          << " (V" << mProfile.v0() << "->V" << mProfile.v1() << ")");

  // Each vertex constraint is an equation of the form: Ai * v = bi
  // Where 'v' is a vector representing the vertex,
  // 'Ai' is a (row) vector and 'bi' a scalar.
  //
  // The vertex is completely determined with 3 such constraints,
  // so is the solution to the folloing system:
  //
  //  A.r0(). * v = b0
  //  A1 * v = b1
  //  A2 * v = b2
  //
  // Which in matrix form is :  A * v = b
  //
  // (with 'A' a 3x3 matrix and 'b' a vector)
  //
  // The member variable mConstrinas contains A and b. Indidivual constraints (Ai,bi) can be added to it.
  // Once 3 such constraints have been added 'v' is directly solved a:
  //
  //  v = b*inverse(A)
  //
  // A constraint (Ai,bi) must be alpha-compatible with the previously added constraints (see Paper); if it's not, is discarded.
  //
  if(mBdry_data.size() > 0)
    add_boundary_preservation_constraints(mBdry_data);

  if(mConstraints_n < 3)
    add_volume_preservation_constraints(mTriangle_data);

  if(mConstraints_n < 3)
    add_boundary_and_volume_optimization_constraints(mBdry_data,mTriangle_data);

  if(mConstraints_n < 3)
    add_shape_optimization_constraints(mProfile.link());

  // It might happen that there were not enough alpha-compatible constraints.
  // In that case there is simply no good vertex placement
  if(mConstraints_n == 3)
  {
    // If the matrix is singular it's inverse cannot be computed so an 'absent' value is returned.
    boost::optional<Matrix> lOptional_Ai = inverse_matrix(mConstraints_A);
    if(lOptional_Ai)
    {
      const Matrix& lAi = *lOptional_Ai;

      CGAL_SMS_LT_TRACE(2,"       b: " << xyz_to_string(mConstraints_b));
      CGAL_SMS_LT_TRACE(2,"  inv(A): " << matrix_to_string(lAi));

      lPlacementV = filter_infinity(mConstraints_b * lAi);

      CGAL_SMS_LT_TRACE(0,"  New vertex point: " << xyz_to_string(*lPlacementV));
    }
    else
    {
      CGAL_SMS_LT_TRACE(1,"  Can't solve optimization, singular system.");
    }
  }
  else
  {
    CGAL_SMS_LT_TRACE(1,"  Can't solve optimization, not enough alpha-compatible constraints.");
  }

  if(lPlacementV)
    rPlacementP = ORIGIN + (*lPlacementV);

  return rPlacementP;
}

template<class TM, class K>
typename LindstromTurkCore<TM,K>::Optional_FT
LindstromTurkCore<TM,K>::
compute_cost(const Optional_point& aP)
{
  Optional_FT rCost;

  if(aP)
  {
    CGAL_SMS_LT_TRACE(0,"Computing LT cost for E" << mProfile.v0_v1());
    Vector lV = (*aP) - ORIGIN;

    FT lSquaredLength = squared_distance(mProfile.p0(), mProfile.p1());

    FT lBdryCost = compute_boundary_cost(lV, mBdry_data);
    FT lVolumeCost = compute_volume_cost(lV, mTriangle_data);
    FT lShapeCost  = compute_shape_cost(*aP, mProfile.link());

    FT lTotalCost = FT(mParams.m_volume_weight) * lVolumeCost
                    + FT(mParams.m_boundary_weight) * lBdryCost * lSquaredLength
                    + FT(mParams.m_shape_weight) * lShapeCost  * lSquaredLength * lSquaredLength;

    rCost = filter_infinity(lTotalCost);

    CGAL_SMS_LT_TRACE(0, "    Squared edge length: " << n_to_string(lSquaredLength)
                           << "\n    Boundary cost: " << n_to_string(lBdryCost) << " weight: " << mParams.m_boundary_weight
                           << "\n    Volume cost: " << n_to_string(lVolumeCost) << " weight: " << mParams.m_volume_weight
                           << "\n    Shape cost: " << n_to_string(lShapeCost)   << " weight: " << mParams.m_shape_weight
                           << "\n  TOTAL COST: " << n_to_string(lTotalCost));

  }

  return rCost;
}

template<class TM, class K>
void
LindstromTurkCore<TM,K>::
add_boundary_preservation_constraints(const Boundary_data_vector& aBdry)
{
  if(aBdry.size() > 0)
  {
    Vector e1 = NULL_VECTOR;
    Vector e2 = NULL_VECTOR;

    FT e1x = FT(0),
       e1y = FT(0),
       e1z = FT(0);

    for(typename Boundary_data_vector::const_iterator it = aBdry.begin(); it != aBdry.end(); ++it)
    {
      e1 = e1 + it->v;
      e2 = e2 + it->n;

      FT vx = it->v.x();
      FT vy = it->v.y();
      FT vz = it->v.z();

      e1x = e1x + vx;
      e1y = e1y + vy;
      e1z = e1z + vz;

      CGAL_SMS_LT_TRACE(1,"    vx:" << n_to_string(vx) << " vy:" << n_to_string(vy) << " vz:" << n_to_string(vz) << " e1x:"
                        << n_to_string(e1x) << " e1y:" << n_to_string(e1y) << " e1z:" << n_to_string(e1z));
    }

    Matrix H = LT_product(e1);
    Vector c = cross_product(e1,e2);

    CGAL_SMS_LT_TRACE(1, "  Adding boundary preservation constraint."
                           << "\n      SumV:" << xyz_to_string(e1)
                           << "\n      SumN:" << xyz_to_string(e2)
                           << "\n      H:" << matrix_to_string(H)
                           << "\n      c:" << xyz_to_string(c));

    add_constraint_from_gradient(H,c);
  }
}

template<class TM, class K>
void
LindstromTurkCore<TM,K>::
add_volume_preservation_constraints(const Triangle_data_vector& triangles)
{
  CGAL_SMS_LT_TRACE(1,"  Adding volume preservation constraint. " << triangles.size() << " triangles.");

  Vector lSumV = NULL_VECTOR;
  FT lSumL(0);

  for(typename Triangle_data_vector::const_iterator it = triangles.begin(), eit = triangles.end(); it != eit; ++it)
  {
    lSumV = lSumV + it->NormalV;
    lSumL = lSumL + it->NormalL;
  }

  CGAL_SMS_LT_TRACE(1, "      SumV:" << xyz_to_string(lSumV) << "\n      SumL:" << n_to_string(lSumL));

  add_constraint_if_alpha_compatible(lSumV,lSumL);
}

template<class TM, class K>
void
LindstromTurkCore<TM,K>::
add_boundary_and_volume_optimization_constraints(const Boundary_data_vector& aBdry,
                                                 const Triangle_data_vector& triangles)
{
  CGAL_SMS_LT_TRACE(1,"  Adding boundary and volume optimization constraints. ");

  Matrix H = NULL_MATRIX;
  Vector c = NULL_VECTOR;

  // Volume optimization
  for(const Triangle_data& lTri : triangles)
  {
    H += direct_product(lTri.NormalV, lTri.NormalV);
    c = c - (lTri.NormalL * lTri.NormalV);
  }

  CGAL_SMS_LT_TRACE(2,"      Hv:" << matrix_to_string(H) << "\n      cv:" << xyz_to_string(c));

  if(aBdry.size() > 0)
  {
    // Boundary optimization
    Matrix Hb = NULL_MATRIX;
    Vector cb = NULL_VECTOR;

    for(typename Boundary_data_vector::const_iterator it = aBdry.begin(); it != aBdry.end(); ++it)
    {
      Matrix H = LT_product(it->v);
      Vector c = cross_product(it->v, it->n);

      Hb += H;
      cb = cb + c;
    }

    CGAL_SMS_LT_TRACE(2,"      Hb:" << matrix_to_string(Hb) << "\n      cb:" << xyz_to_string(cb));

    // Weighted average
    FT scaled_boundary_weight = FT(9) * FT(mParams.m_boundary_weight) * squared_distance(mProfile.p0(),mProfile.p1()) ;

    H *= mParams.m_volume_weight;
    c = c * mParams.m_volume_weight;

    H += scaled_boundary_weight * Hb;
    c = c + (scaled_boundary_weight * cb);

    CGAL_SMS_LT_TRACE(2, "      H:" << matrix_to_string(H) << "\n      c:" << xyz_to_string(c));
    CGAL_SMS_LT_TRACE(2, "      VolW:" << mParams.m_volume_weight << " BdryW:" << mParams.m_boundary_weight << " ScaledBdryW:" << scaled_boundary_weight);
  }

  add_constraint_from_gradient(H,c);
}

template<class TM, class K>
void
LindstromTurkCore<TM,K>::
add_shape_optimization_constraints(const vertex_descriptor_vector& link)
{
  FT s(double(link.size()));

  Matrix H(s, 0.0, 0.0,
           0.0, s, 0.0,
           0.0, 0.0, s);

  Vector c = NULL_VECTOR;
  for(const vertex_descriptor v : link)
    c = c + (ORIGIN - get_point(v));

  CGAL_SMS_LT_TRACE(1,"  Adding shape optimization constraint. Shape vector: " << xyz_to_string(c));

  add_constraint_from_gradient(H,c);
}

template<class TM, class K>
typename LindstromTurkCore<TM,K>::FT
LindstromTurkCore<TM,K>::
compute_boundary_cost(const Vector& v,
                      const Boundary_data_vector& aBdry)
{
  FT rCost(0);
  for(typename Boundary_data_vector::const_iterator it = aBdry.begin(); it != aBdry.end(); ++it)
  {
    Vector u = (it->t - ORIGIN) - v;
    Vector c = cross_product(it->v, u);
    rCost += c*c;
  }

  return rCost / FT(4);
}

template<class TM, class K>
typename LindstromTurkCore<TM,K>::FT
LindstromTurkCore<TM,K>::
compute_volume_cost(const Vector& v,
                    const Triangle_data_vector& triangles)
{
  FT rCost(0);

  for(const Triangle_data& lTri : triangles)
  {
    FT lF = lTri.NormalV * v - lTri.NormalL;
    rCost += (lF * lF);
  }

  return rCost / FT(36);
}

template<class TM, class K>
typename LindstromTurkCore<TM,K>::FT
LindstromTurkCore<TM,K>::
compute_shape_cost(const Point& p,
                   const vertex_descriptor_vector& link)
{
  FT rCost(0);
  for(const vertex_descriptor v : link)
    rCost += squared_distance(p, get_point(v));

  return rCost;
}

template<class TM, class K>
void
LindstromTurkCore<TM,K>::
add_constraint_if_alpha_compatible(const Vector& Ai,
                                   const FT& bi)
{
  CGAL_SMS_LT_TRACE(3, "    Adding new constraints if alpha-compatible.\n      Ai: " << xyz_to_string(Ai) << "\n      bi:" << n_to_string(bi) << ")");

  FT slai = Ai*Ai;
  CGAL_SMS_LT_TRACE(3, "\n      slai: " << n_to_string(slai) << ")");

  if(is_finite(slai))
  {
    FT l = CGAL_NTS sqrt(slai);
    CGAL_SMS_LT_TRACE(3, "      l: " << n_to_string(l));

    if(!CGAL_NTS is_zero(l))
    {
      Vector Ain = Ai / l;
      FT bin = bi / l;
      CGAL_SMS_LT_TRACE(3, "      Ain: " << xyz_to_string(Ain) << " bin:" << n_to_string(bin));

      bool lAddIt = true;

      if(mConstraints_n == 1)
      {
        FT d01 = mConstraints_A.r0() * Ai ;
        FT sla0 = mConstraints_A.r0() * mConstraints_A.r0();
        FT sd01 = d01 * d01;
        FT max = sla0 * slai * mSquared_cos_alpha;

        CGAL_SMS_LT_TRACE(3, "      Second constraint. d01: " << n_to_string(d01) << " sla0:" << n_to_string(sla0)
                                                 << " sd01:" << n_to_string(sd01) << " max:" << n_to_string(max));

        if(sd01 > max)
          lAddIt = false;
      }
      else if(mConstraints_n == 2)
      {
        Vector N = cross_product(mConstraints_A.r0(),mConstraints_A.r1());

        FT dc012 = N * Ai;
        FT slc01 = N * N;
        FT sdc012 = dc012 * dc012;

        FT min = slc01 * slai * mSquared_sin_alpha;

        CGAL_SMS_LT_TRACE(3, "      Third constraint. N: " << xyz_to_string(N)
                               << " dc012:" << n_to_string(dc012) << " slc01:" << n_to_string(slc01)
                               << " sdc012:" << n_to_string(sdc012) << " min:" << n_to_string(min));

        if(sdc012 <= min)
          lAddIt = false;
      }

      if(lAddIt)
      {
        switch(mConstraints_n)
        {
          case 0:
            mConstraints_A.r0() = Ain;
            mConstraints_b = Vector(bin,mConstraints_b.y(),mConstraints_b.z());
            break;
          case 1:
            mConstraints_A.r1() = Ain;
            mConstraints_b = Vector(mConstraints_b.x(),bin,mConstraints_b.z());
            break;
          case 2:
            mConstraints_A.r2() = Ain;
            mConstraints_b = Vector(mConstraints_b.x(),mConstraints_b.y(),bin);
            break;
        }

        CGAL_SMS_LT_TRACE(3, "      Accepting # " << mConstraints_n << " A:" << matrix_to_string(mConstraints_A) << " b:" << xyz_to_string(mConstraints_b));

        ++mConstraints_n;
      }
      else
      {
        CGAL_SMS_LT_TRACE(3, "      INCOMPATIBLE. Discarded");
      }
    }
    else
    {
      CGAL_SMS_LT_TRACE(3, "        l is ZERO. Discarded");
    }
  }
  else
  {
    CGAL_SMS_LT_TRACE(3, "      OVERFLOW. Discarded");
  }
}

template<class K>
int index_of_max_component(const typename K::Vector_3& v)
{
  typedef typename K::FT                     FT;

  int i = 0;
  FT max = v.x();

  if(max < v.y())
  {
    max = v.y();
    i = 1;
  }

  if(max < v.z())
  {
    max = v.z();
    i = 2;
  }

  return i;
}

template<class TM, class K>
void
LindstromTurkCore<TM,K>::
add_constraint_from_gradient(const Matrix& H,
                             const Vector& c)
{
  CGAL_SMS_LT_TRACE(3, "    Adding constraint from gradient. Current n=" << mConstraints_n);

  CGAL_precondition(mConstraints_n >= 0 && mConstraints_n <= 2);

  switch(mConstraints_n)
  {
    case 0:
      add_constraint_if_alpha_compatible(H.r0(), -c.x());
      add_constraint_if_alpha_compatible(H.r1(), -c.y());
      add_constraint_if_alpha_compatible(H.r2(), -c.z());
      break;
    case 1:
    {
      const Vector& A0 = mConstraints_A.r0();

      CGAL_assertion(A0 != NULL_VECTOR);

      Vector AbsA0(CGAL::abs(A0.x()),
                   CGAL::abs(A0.y()),
                   CGAL::abs(A0.z()));

      Vector Q0;
      switch(index_of_max_component<Geom_traits>(AbsA0))
      {
        // Since A0 is guaranteed to be non-zero, the denominators here are known to be non-zero too.

        case 0:
          Q0 = Vector(- A0.z()/A0.x(), 0, 1);
          break;
        case 1:
          Q0 = Vector(0,- A0.z()/A0.y(), 1);
          break;
        case 2:
          Q0 = Vector(1, 0, - A0.x()/A0.z());
          break;
        default:
          Q0 = NULL_VECTOR; // This should never happen!
      }

      CGAL_SMS_LT_TRACE(3, "      Q0:" << xyz_to_string(Q0));

      CGAL_assertion(Q0 != NULL_VECTOR);

      Vector Q1 = cross_product(A0, Q0);
      Vector A1 = H * Q0;
      Vector A2 = H * Q1;

      FT b1 = - (Q0 * c);
      FT b2 = - (Q1 * c);

      CGAL_SMS_LT_TRACE(3, "      Q1:" << xyz_to_string(Q1) << "\n      A1: " << xyz_to_string(A1) << "\n      A2:" << xyz_to_string(A2)
                             << "\n      b1:" << n_to_string(b1) << "\n      b2:" << n_to_string(b2));

      add_constraint_if_alpha_compatible(A1,b1);
      add_constraint_if_alpha_compatible(A2,b2);
    }
      break;

    case 2:
    {
      Vector Q = cross_product(mConstraints_A.r0(),mConstraints_A.r1());
      Vector A2 = H * Q;
      FT b2 = - (Q * c);
      CGAL_SMS_LT_TRACE(3, "      Q:" << xyz_to_string(Q) << "\n      A2: " << xyz_to_string(A2) << "\n      b2:" << n_to_string(b2));

      add_constraint_if_alpha_compatible(A2,b2);
    }
      break;
  }
}

} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_DETAIL_LINDSTROM_TURK_CORE_H
