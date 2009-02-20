// Copyright (c) 2007-09  ETH Zurich (France).
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

#ifndef CGAL_APSS_RECONSTRUCTION_FUNCTION_H
#define CGAL_APSS_RECONSTRUCTION_FUNCTION_H

#include <vector>
#include <algorithm>

#define MLS_USE_CUSTOM_KDTREE

#include <CGAL/Point_with_normal_3.h>
#include <CGAL/make_surface_mesh.h>
#ifdef MLS_USE_CUSTOM_KDTREE
#include "Neighborhood/knn_point_neighbor_search.h"
#else
#include <CGAL/Orthogonal_k_neighbor_search.h>
#endif
#include <CGAL/Search_traits_3.h>
#include <CGAL/Surface_mesher/Implicit_surface_oracle_3.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>
#include <CGAL/surface_reconstruction_assertions.h>

CGAL_BEGIN_NAMESPACE


/// APSS_reconstruction_function computes an implicit function
/// that defines a Point Set Surface (PSS) based on
/// moving least squares (MLS) fitting of algebraic spheres.
/// See "Algebraic Point Set Surfaces" by Guennebaud and Gross [Guennebaud07].
///
/// @heading Is Model for the Concepts:
/// Model of the ImplicitFunction concept.
///
/// @heading Design Pattern:
/// A model of ImplicitFunction is a
/// Strategy [GHJV95]: it implements a strategy of surface mesh reconstruction.
///
/// @heading Parameters:
/// @param Gt Geometric traits class.

template <class Gt>
class APSS_reconstruction_function
{
// Public types
public:

  typedef Gt Geom_traits; ///< Kernel's geometric traits

  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Vector_3 Vector;
  typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid;
  typedef typename Geom_traits::Sphere_3 Sphere;

  typedef Point_with_normal_3<Gt> Point_with_normal; ///< == Point_with_normal_3<Gt> (compact representation)
  typedef typename Point_with_normal::Normal Normal; ///< == Vector_3

// Private types
private:

  // Item in the Kd-tree: position (Point_3) + normal + index
  class KdTreeElement : public Point_with_normal
  {
  public:
    unsigned int index;

    KdTreeElement(const Origin& o = ORIGIN, unsigned int id=0)
      : Point_with_normal(o), index(id)
    {}
    template <class K, class N>
    KdTreeElement(const Point_with_normal_3<K,N>& pwn, unsigned int id=0)
      : Point_with_normal(pwn), index(id)
    {}
    KdTreeElement(const Point& p, unsigned int id=0)
      : Point_with_normal(p), index(id)
    {}
    KdTreeElement(const KdTreeElement& other)
      : Point_with_normal(other), index(other.index)
    {}
  };

  // Helper class for the Kd-tree
  class KdTreeGT : public Geom_traits
  {
  public:
    typedef KdTreeElement Point_3;
  };

  class TreeTraits : public Search_traits_3<KdTreeGT>
  {
    public:
      typedef Point PointType;
  };
  #ifdef MLS_USE_CUSTOM_KDTREE
  typedef knn_point_neighbor_search<TreeTraits> Neighbor_search;
  #else
  typedef Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  #endif
  typedef typename Neighbor_search::Tree Tree;

//   typedef Orthogonal_k_neighbor_search<TreeTraits> CNeighbor_search;
//   typedef typename CNeighbor_search::Tree CTree;

  typedef typename Neighbor_search::Point_with_transformed_distance
                                    Point_with_transformed_distance;

// Public methods
public:

  /// Create an APSS implicit function from a point set.
  ///
  /// Precondition: the value type of InputIterator must be convertible to Point_with_normal_3.
  ///
  /// @param first First point to add.
  /// @param beyond Past-the-end point to add.
  /// @param k Number of nearest neighbors.
  /// @param projection_error Dichotomy error when projecting point.
  template < class InputIterator >
  APSS_reconstruction_function(InputIterator first, InputIterator beyond,
                               unsigned int k,
                               FT projection_error = 3.16e-4) // sqrt(1e-7)
  {
    // Allocate smart pointer to data
    m = new Private;

    int nb_points = std::distance(first, beyond);

    // Number of nearest neighbors
    m->nofNeighbors = k;

    // Create kd-tree
    m->treeElements.reserve(nb_points);
    unsigned int i=0;
    for (InputIterator it=first ; it != beyond ; ++it,++i)
    {
      m->treeElements.push_back(KdTreeElement(*it,i));
    }
    #ifdef MLS_USE_CUSTOM_KDTREE
    m->tree = new Tree(ConstDataWrapper<Point>(static_cast<const Point*>(&(m->treeElements[0])), m->treeElements.size(), sizeof(KdTreeElement)), &(m->treeElements[0]));
    #else
    m->tree = new Tree(m->treeElements.begin(), m->treeElements.end());
    #endif
//     CTree* ctree = new CTree(m->treeElements.begin(), m->treeElements.end());

    // Compute the radius of each point = (distance max to KNN)/2.
    // The union of these balls defines the surface definition domain.
    m->radii.resize(nb_points);
//     for (int y=0; y<40; ++y)
    {
      int i=0;
      for (InputIterator it=first ; it != beyond ; ++it, ++i)
      {
        Neighbor_search search(*(m->tree), *it, 16); // why 16?
        FT maxdist2 = (--search.end())->second; // squared distance to furthest neighbor
        m->radii[i] = sqrt(maxdist2)/2.;

//         for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
//           std::cout << it->first.index << " ";
//         std::cout << "\n";
//         CNeighbor_search csearch(*ctree, *it, 8);
//         for (typename CNeighbor_search::iterator it = csearch.begin(); it != csearch.end(); ++it)
//           std::cout << it->first.index << " ";
//         std::cout << "\n\n";
      }
    }
//     exit(0);

    // Compute barycenter, bounding box, bounding sphere and standard deviation.
    update_bounds(first, beyond);

    // Find a point inside the surface.
    find_inner_point();

    // Dichotomy error when projecting point (squared)
    m->sqError = projection_error * projection_error * Gt().compute_squared_radius_3_object()(m->bounding_sphere);
  }

  /// Copy constructor
  APSS_reconstruction_function(const APSS_reconstruction_function& other) {
    m = other.m;
    m->count++;
  }

  /// operator =()
  APSS_reconstruction_function& operator = (const APSS_reconstruction_function& other) {
    m = other.m;
    m->count++;
  }

  /// Destructor
  ~APSS_reconstruction_function() {
    if (--(m->count)==0)
      delete m;
  }

  void set_numbers_of_neighbors(unsigned int k) {m->nofNeighbors = k;}

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

private:

  /** Fit an algebraic sphere on a set of neigbors in a Moving Least Square sense.
  The weight function is scaled such that the weight of the furthest neighbor is 0.
  */
  void fit(const Neighbor_search& search) const
  {
    FT r2 = (--search.end())->second; // squared distance to furthest neighbor

    Vector sumP(0,0,0);
    Vector sumN(0,0,0);
    FT sumDotPP = 0.;
    FT sumDotPN = 0.;
    FT sumW = 0.;

    r2 *= 1.001;
    FT invr2 = 1./r2;
    for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
      Vector p = it->first - CGAL::ORIGIN;
      const Vector& n = it->first.normal();
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

  /** Check whether the point 'p' is close to the input points or not.
      We assume that it's not if it is far from its nearest neighbor.
  */
  inline bool isValid(const Point_with_transformed_distance& nearest_neighbor, const Point& /* p */) const
  {
      FT r = 2. * m->radii[nearest_neighbor.first.index];
      return (r*r > nearest_neighbor.second);
  }

public:

  /// [ImplicitFunction interface]
  ///
  /// Evaluate implicit function for any 3D point.
  //
  // Implementation note: this function is called a large number of times,
  // thus us heavily optimized. The bottleneck is Neighbor_search's constructor,
  // which we try to avoid calling.
  FT operator()(const Point& p) const
  {
    // Is 'p' close to the surface?
    // Optimization: test first if 'p' is close to one of the neighbors
    //               computed during the previous call.
    typename Geom_traits::Compute_squared_distance_3 sqd;
    m->cached_nearest_neighbor.second = sqd(p, m->cached_nearest_neighbor.first);
    if (!isValid(m->cached_nearest_neighbor, p))
    {
      // Compute the nearest neighbor and cache it
      KdTreeElement query(p);
      Neighbor_search search_1nn(*(m->tree), query, 1);
      m->cached_nearest_neighbor = *(search_1nn.begin());

      // Is 'p' close to the surface?
      if (!isValid(m->cached_nearest_neighbor, p))
      {
        // If 'p' is far from the surface, project its nearest neighbor onto the surface...
        Vector n;
        Point pp = m->cached_nearest_neighbor.first;
        project(pp,n,1);
        // ...and return the (signed) distance to the surface
        Vector h = sub(p,pp);
        return length(h) * ( dot(n,h)>0. ? 1. : -1.);
      }
    }

    // Compute k nearest neighbors and cache the first one
    KdTreeElement query(p);
    Neighbor_search search_knn(*(m->tree), query, m->nofNeighbors);
    m->cached_nearest_neighbor = *(search_knn.begin());

    // If 'p' is close to the surface, fit an algebraic sphere
    // on a set of neigbors in a Moving Least Square sense.
    fit(search_knn);

    // return the distance to the sphere
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
      if (!isValid(*(search.begin()),p))
      {
        p = search.begin()->first;
        n = search.begin()->first.normal();
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
          if (!isValid(*(search.begin()),p))
          {
            // if the closest intersection is far away the input points,
            // then we take the other one.
            p = add(m->as.center,mul(-m->as.radius,dir));
            flipN = -flipN;
          }
          n = mul(flipN,dir);

          if (!isValid(*(search.begin()),p))
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
  void update_bounds(InputIterator first, InputIterator beyond)
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
    m->bounding_sphere = Sphere(ms3.center(), ms3.squared_radius());

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

    /** Converts the algebraic sphere to an explicit sphere or plane.
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

  struct Private
  {
    Private()
      : tree(NULL), count(1)
    {}

    ~Private()
    {
      delete tree; tree = NULL;
    }

    Tree* tree;
    std::vector<KdTreeElement> treeElements;
    std::vector<FT> radii;
    Iso_cuboid bounding_box; // Points' bounding box
    Sphere bounding_sphere; // Points' bounding sphere
    Point barycenter; // Points' barycenter
    FT diameter_standard_deviation; // Standard deviation of the distance to barycenter
    FT sqError; // Dichotomy error when projecting point (squared)
    unsigned int nofNeighbors; // Number of nearest neighbors
    Point inner_point; // Point inside the surface
    mutable AlgebraicSphere as;
    mutable Point_with_transformed_distance cached_nearest_neighbor;
    int count; // reference counter
  };

  Private* m; // smart pointer to data

}; // end of APSS_reconstruction_function


CGAL_END_NAMESPACE

#endif // CGAL_APSS_RECONSTRUCTION_FUNCTION_H
