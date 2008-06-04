// Copyright (c) 2007-2008  ETH Zurich (France).
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
// Author(s)     : Gael Guennebaud, Laurent Saboret


#ifndef CGAL_APSS_IMPLICIT_FUNCTION_H
#define CGAL_APSS_IMPLICIT_FUNCTION_H

#include <vector>
#include <algorithm>

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Surface_mesher/Implicit_surface_oracle_3.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>
#include <CGAL/surface_reconstruction_assertions.h>

CGAL_BEGIN_NAMESPACE


/// APSS_implicit_function computes an implicit function
/// that defines a Point Set Surface (PSS) based on
/// moving least squares (MLS) fitting of algebraic spheres.
/// See "Algebraic Point Set Surfaces" by Guennebaud and Gross (2007).
///
/// @heading Is Model for the Concepts:
/// Model of the Reconstruction_implicit_function concept.
///
/// @heading Design Pattern:
/// A model of ReconstructionImplicitFunction is a
/// Strategy [GHJV95]: it implements a strategy of surface mesh reconstruction.
///
/// @heading Parameters:
/// @param Gt Geometric traits class.
/// @param PointWithNormal_3 Model of PointWithNormal_3 concept.

template <class Gt, class PointWithNormal_3>
class APSS_implicit_function
{
// Public types
public:

  typedef Gt Geom_traits; ///< Kernel's geometric traits

  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;
  typedef typename Geom_traits::Sphere_3 Sphere;

  typedef PointWithNormal_3 Point_with_normal;       ///< Model of PointWithNormal_3
  typedef typename Point_with_normal::Normal Normal; ///< Model of OrientedNormal_3 concept.
  typedef typename Geom_traits::Vector_3 Vector;

// Private types
private:

  // Item in the Kd-tree: position (Point_3) + normal + index
  class KdTreeElement : public Point_with_normal
  {
  public:
    unsigned int index;

    KdTreeElement() {}
    KdTreeElement(const Point_with_normal& pwn, unsigned int id=0)
      : Point_with_normal(pwn), index(id)
    {}
    KdTreeElement(const Point& p, unsigned int id=0)
      : Point_with_normal(p), index(id)
    {}
  };

  // Helper class for the Kd-tree
  class KdTreeGT : public Geom_traits
  {
  public:
    typedef KdTreeElement Point_3;
  };

  typedef Search_traits_3<KdTreeGT> TreeTraits;
  typedef Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

// Public methods
public:

  /// Create an APSS implicit function from a point set.
  ///
  /// Precondition: the value type of InputIterator must be convertible to Point_with_normal.
  ///
  /// @param first First point of point set.
  /// @param beyond Past-the-end point of point set.
  template < class InputIterator >
  APSS_implicit_function(InputIterator first, InputIterator beyond)
  {
    m = new Private;

    unsigned int i=0;
    m->treeElements.reserve(std::distance(first, beyond));
    for (InputIterator it=first ; it != beyond ; ++it,++i)
    {
      m->treeElements.push_back(KdTreeElement(*it,i));
    }
    m->tree = new Tree(m->treeElements.begin(), m->treeElements.end());

    // Compute the radius of each point = (distance max to KNN)/2.
    // The union of these ball defines the surface definition domain.
    i = 0;
    for (InputIterator it=first ; it != beyond ; ++it,++i)
    {
      KdTreeElement query(*it);
      Neighbor_search search(*(m->tree), query, 16);
      FT maxdist2 = 0.;
      for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
      {
        if (maxdist2<it->second)
          maxdist2 = it->second;
      }
      m->radii.push_back(2.*sqrt(maxdist2/16.));
    }
    
    // Compute barycenter, bounding box, bounding sphere and standard deviation.
    update_bounding_box(first, beyond);

    // Find a point inside the surface.
    find_inner_point();

    m->sqError = 0.0000001 * Gt().compute_squared_radius_3_object()(m->bounding_sphere);
  }

  APSS_implicit_function(const APSS_implicit_function& other) {
    m = other.m;
    m->count++;
  }

  APSS_implicit_function& operator = (const APSS_implicit_function& other) {
    m = other.m;
    m->count++;
  }

  ~APSS_implicit_function() {
    if (--(m->count)==0)
      delete m;
  }

  void setNofNeighbors(unsigned int k) {m->nofNeighbors = k;}

  /// Get the bounding box.
  Iso_cuboid bounding_box() const
  {
    return m->bounding_box;
  }

  /// Get the surface's bounding sphere.
  const Sphere& bounding_sphere() const
  {
    return m->bounding_sphere;
  }

  /// Get the region of interest, ignoring the outliers.
  /// This method is used to define the OpenGL arcball sphere.
  Sphere region_of_interest() const
  {
    // A good candidate is a sphere containing the dense region of the point cloud:
    // - center point is barycenter
    // - Radius is 2 * standard deviation
    FT radius = (FT)2 * (FT)m->diameter_standard_deviation;
    return Sphere(m->barycenter, radius*radius);
  }

  /// You should call compute_implicit_function() once when points insertion is over.
  /// Return false on error.
  bool compute_implicit_function()
  {
    // Nothing to do
    return true;
  }

private:

  /** Fit an algebraic sphere on a set of neigbors in a Moving Least Square sense.
  The weight function is scaled such that the weight of the furthest neighbor is 0.
  */
  void fit(Neighbor_search& neighbors) const
  {
    typename Neighbor_search::iterator last = neighbors.end();
    last--;
    FT r2 = last->second;

    Vector sumP(0,0,0);
    Vector sumN(0,0,0);
    FT sumDotPP = 0.;
    FT sumDotPN = 0.;
    FT sumW = 0.;

    r2 *= 1.001;
    FT invr2 = 1./r2;
    for (typename Neighbor_search::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
    {
      //const Point& p = m->points[it->first.index];
      Vector p = it->first - CGAL::ORIGIN;
      const Vector& n = it->first.normal().get_vector();
      FT w = 1. - it->second*invr2;
      w = w*w; w = w*w;

      sumP = add(sumP,mul(w,p));
      sumN = add(sumN,mul(w,n));
      sumDotPP += w * dot(p,p);
      sumDotPN += w * dot(p,n);
      sumW += w;
    }

    FT invSumW = 1./sumW;
    m->as.u4 = 0.5 * (sumDotPN - invSumW*dot(sumP,sumN))/(sumDotPP - invSumW*dot(sumP,sumP));
    m->as.u13 = mul(invSumW,add(sumN,mul(-2.*m->as.u4,sumP)));
    m->as.u0 = -invSumW*(dot(m->as.u13,sumP)+m->as.u4*sumDotPP);
    m->as.finalize();
  }

  /** Check whether the point p is close to the input points or not.
  */
  inline bool isValid(const Neighbor_search& neighbors, const Point& /* p */) const
  {
    for (typename Neighbor_search::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
    {
      FT r = 2. * m->radii[it->first.index];
      if (r*r > it->second)
        return true;
    }
    return false;
  }

public:

  /// [ImplicitFunction interface]
  ///
  /// Evaluate implicit function for any 3D point.
  FT operator()(const Point& p) const
  {
    KdTreeElement query(p);
    Neighbor_search search(*(m->tree), query, m->nofNeighbors);
    //unsigned int closestId = search.begin()->first.index;

    if (!isValid(search,p))
    {
      Vector n;
      Point pp = p;
      project(pp,n,1);
      Vector h = sub(p,pp);
      return length(h) * ( dot(n,h)>0. ? 1. : -1.);
    }

    fit(search);
    return m->as.euclideanDistance(p);
  }

  /// Get point inside the surface.
  Point get_inner_point() const
  {
    return m->inner_point;
  }

// Private methods:
private:

  /** Projects the point p onto the MLS surface.
  */
  void project(Point& p, unsigned int maxNofIterations = 20) const
  {
    Vector n;
    project(p,n,maxNofIterations);
  }

  /** Projects the point p onto the MLS surface, and returns an approximate normal.
  */
  void project(Point& p, Vector& n, unsigned int maxNofIterations = 20) const
  {
    Point source = p;

    FT delta2 = 0.;
    unsigned int countIter = 0;
    do {

      Neighbor_search search(*(m->tree), p, m->nofNeighbors);

      // if p is far away the input point cloud, start with the closest point.
      if (!isValid(search,p))
      {
        p = search.begin()->first;
        n = search.begin()->first.normal().get_vector();
        delta2 = search.begin()->second;
      }
      else
      {
        fit(search);

        Point oldP = p;
        if (m->as.state==AlgebraicSphere::SPHERE)
        {
          // projection onto a sphere.
          Vector dir = normalize(sub(source,m->as.center));
          p = add(m->as.center,mul(m->as.radius,dir));
          FT flipN = m->as.u4<0. ? -1. : 1.;
          if (!isValid(search,p))
          {
            // if the closest intersection is far away the input points,
            // then we take the other one.
            p = add(m->as.center,mul(-m->as.radius,dir));
            flipN = -flipN;
          }
          n = mul(flipN,dir);

          if (!isValid(search,p))
          {
            std::cout << "Invalid projection\n";
          }
        }
        else if (m->as.state==AlgebraicSphere::PLANE)
        {
          // projection onto a plane.
          p = sub(source, mul(dot(m->as.normal,source-CGAL::ORIGIN)+m->as.d,m->as.normal));
        }
        else
        {
          // iterative projection onto the algebraic sphere.
          p = m->as.iProject(source);
        }

        Vector diff = sub(oldP,p);
        delta2 = dot(diff,diff);
      }

    } while ( ((++countIter)<maxNofIterations) && (delta2<m->sqError) );
  }

  inline static FT dot(const Vector& a, const Vector& b) {
    return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
  }
  inline static FT length(const Vector& a) {
    return sqrt(dot(a,a));
  }
  inline static Vector mul(FT s, const Vector& p) {
    return Vector(p.x()*s, p.y()*s, p.z()*s);
  }
  inline static Point add(FT s, const Point& p) {
    return Point(p.x()+s, p.y()+s, p.z()+s);
  }
  inline static Point add(const Point& a, const Vector& b) {
    return Point(a.x()+b.x(), a.y()+b.y(), a.z()+b.z());
  }
  inline static Vector add(const Vector& a, const Vector& b) {
    return Vector(a.x()+b.x(), a.y()+b.y(), a.z()+b.z());
  }
  inline static Point sub(const Point& a, const Vector& b) {
    return Point(a.x()-b.x(), a.y()-b.y(), a.z()-b.z());
  }
  inline static Vector sub(const Point& a, const Point& b) {
    return Vector(a.x()-b.x(), a.y()-b.y(), a.z()-b.z());
  }
  inline static Vector normalize(const Vector& p) {
    FT s = 1. / length(p);
    return mul(s,p);
  }
  inline static Vector cross(const Vector& a, const Vector& b) {
    return Vector(a.y()*b.z() - a.z()*b.y(),
      a.z()*b.x() - a.x()*b.z(),
      a.x()*b.y() - a.y()*b.x());
  }

  /// Compute barycenter, bounding box, bounding sphere and standard deviation.
  ///
  /// @param first First point of point set.
  /// @param beyond Past-the-end point of point set.
  template < class InputIterator >
  void update_bounding_box(InputIterator first, InputIterator beyond)
  {
    if (first == beyond)
      return;

    // Update bounding box and barycenter.
    // TODO: we should use the functions in PCA component instead.
    FT xmin,xmax,ymin,ymax,zmin,zmax;
    xmin = ymin = zmin =  1e38;
    xmax = ymax = zmax = -1e38;
    Vector v = CGAL::NULL_VECTOR;
    FT norm = 0;
    for (InputIterator it = first; it != beyond; it++)
    {
      Point p = *it;

      // update bbox
      xmin = (std::min)(p.x(),xmin);
      ymin = (std::min)(p.y(),ymin);
      zmin = (std::min)(p.z(),zmin);
      xmax = (std::max)(p.x(),xmax);
      ymax = (std::max)(p.y(),ymax);
      zmax = (std::max)(p.z(),zmax);

      // update barycenter
      v = v + (p - CGAL::ORIGIN);
      norm += 1;
    }
    //
    Point p(xmin,ymin,zmin);
    Point q(xmax,ymax,zmax);
    m->bounding_box = Iso_cuboid(p,q);
    //
    m->barycenter = CGAL::ORIGIN + v / norm;

    // sphere of smallest volume enclosing the input points
    Min_sphere_d< CGAL::Optimisation_d_traits_3<Gt> > ms3(first, beyond);
    m->bounding_sphere = Sphere(ms3.center(), ms3.squared_radius()); // LS 05/2008: 1.2*was sqrt(ms3.squared_radius())

    // Compute standard deviation of the distance to barycenter
    typename Geom_traits::Compute_squared_distance_3 sqd;
    FT sq_radius = 0;
    for (InputIterator it = first; it != beyond; it++)
    {
        sq_radius += sqd(*it, m->barycenter);
    }
    sq_radius /= std::distance(first, beyond);
    m->diameter_standard_deviation = CGAL::sqrt(sq_radius);
  }

private:

  struct AlgebraicSphere {
    FT u0, u4;
    Vector u13;
    enum State {UNDETERMINED=0,PLANE=1,SPHERE=2};
    State state;
    Point center;
    FT radius;
    Vector normal;
    FT d;

    AlgebraicSphere() : state(UNDETERMINED) {}

    /** Converts the algebraix sphere to an explicit sphere or plane.
    */
    void finalize(void) {
      if (fabs(u4)>1e-9)
      {
        state = SPHERE;
        FT b = 1./u4;
        center = CGAL::ORIGIN + mul(-0.5*b,u13);
        radius = sqrt(dot(center - CGAL::ORIGIN,center - CGAL::ORIGIN) - b*u0);
      }
      else if (u4==0.)
      {
        state = PLANE;
        FT s = 1./length(u13);
        normal = mul(s,u13);
        d = u0*s;
      }
      else
      {
        state = UNDETERMINED;
      }
    }

    /** Compute the Euclidean distance between the algebraic surface and a point.
    */
    FT euclideanDistance(const Point& p) {
      if (state==SPHERE)
      {
        FT aux = length(sub(center,p)) - radius;
        if (u4<0.)
          aux = -aux;
        return aux;
      }

      if (state==PLANE)
        return dot(p - CGAL::ORIGIN,normal) + d;

      // else, tedious case, fall back to an iterative method:
      return iEuclideanDistance(p);
    }

    /** Euclidean distance via an iterative projection procedure.
    This is an optimized version of distance(x,iProject(x)).
    */
    inline FT iEuclideanDistance(const Point& x) const
    {
      FT d = 0.;
      Vector grad;
      Vector dir = add(mul(2.*u4,x-CGAL::ORIGIN),u13);
      FT ilg = 1./length(dir);
      dir = mul(ilg,dir);
      FT ad = u0 + dot(u13,x-CGAL::ORIGIN) + u4 * dot(x-CGAL::ORIGIN,x-CGAL::ORIGIN);
      FT delta = -ad*(ilg<1.?ilg:1.);
      Point p = add(x, mul(delta,dir));
      d += delta;
      for (int i=0 ; i<5 ; ++i)
      {
        grad = add(mul(2.*u4,p-CGAL::ORIGIN),u13);
        ilg = 1./length(grad);
        delta = -(u0 + dot(u13,p-CGAL::ORIGIN) + u4 * dot(p-CGAL::ORIGIN,p-CGAL::ORIGIN))*(ilg<1.?ilg:1.);
        p = add(p,mul(delta,dir));
        d += delta;
      }
      return -d;
    }

    /** Iterative projection.
    */
    inline Point iProject(const Point& x) const
    {
      Vector grad;
      Vector dir = add(mul(2.*u4,x-CGAL::ORIGIN),u13);
      FT ilg = 1./length(dir);
      dir = mul(ilg,dir);
      FT ad = u0 + dot(u13,x-CGAL::ORIGIN) + u4 * dot(x-CGAL::ORIGIN,x-CGAL::ORIGIN);
      FT delta = -ad*(ilg<1.?ilg:1.);
      Point p = add(x, mul(delta,dir));
      for (int i=0 ; i<5 ; ++i)
      {
        grad = add(mul(2.*u4,p-CGAL::ORIGIN),u13);
        ilg = 1./length(grad);
        delta = -(u0 + dot(u13,p-CGAL::ORIGIN) + u4 * dot(p-CGAL::ORIGIN,p-CGAL::ORIGIN))*(ilg<1.?ilg:1.);
        p = add(p,mul(delta,dir));
      }
      return p;
    }
  };

  /// Find a random point inside the surface.
  void find_inner_point()
  {
    m->inner_point = CGAL::ORIGIN;
    FT min_f = 1e38;

    // Try random points until we find a point / value < 0
    Point center = m->bounding_sphere.center();
    FT radius = sqrt(m->bounding_sphere.squared_radius());
    CGAL::Random_points_in_sphere_3<Point> rnd(radius);  
    while (min_f > 0)
    {
      // Create random point in bounding sphere
      Point p = center + (*rnd++ - CGAL::ORIGIN);
      FT value = (*this)(p);
      if(value < min_f)
      {
        m->inner_point = p;
        min_f = value;
      }
    }
  }

// Data members
private:

  struct Private {
    Private()
      : nofNeighbors(12), count(1)
    {}
    Tree* tree;
    std::vector<KdTreeElement> treeElements;
    std::vector<FT> radii;
    Iso_cuboid bounding_box; // Points' bounding box
    Sphere bounding_sphere; // Points' bounding sphere
    Point barycenter; // Points' barycenter
    FT diameter_standard_deviation; // Standard deviation of the distance to barycenter
    FT sqError;
    unsigned int nofNeighbors;
    mutable AlgebraicSphere as;
    Point inner_point; // Point inside the surface
    int count;
  };

  Private* m;

}; // end of APSS_implicit_function


CGAL_END_NAMESPACE

#endif // CGAL_APSS_IMPLICIT_FUNCTION_H
