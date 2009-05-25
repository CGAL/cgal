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
// Author(s)     : Gael Guennebaud, Laurent Saboret

#ifndef CGAL_APSS_RECONSTRUCTION_FUNCTION_H
#define CGAL_APSS_RECONSTRUCTION_FUNCTION_H

#include <vector>
#include <algorithm>

#include <CGAL/point_set_property_map.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Fast_orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>
#include <CGAL/surface_reconstruction_points_assertions.h>

CGAL_BEGIN_NAMESPACE


/// APSS_reconstruction_function computes an implicit function
/// that defines a Point Set Surface (PSS) based on
/// moving least squares (MLS) fitting of algebraic spheres.
///
/// This class implements a variant of the "Algebraic Point Set Surfaces" method
/// by Guennebaud and Gross [Guennebaud07].
///
/// Currently, the quality of the reconstruction highly depends on both the quality of input
/// normals and the smoothness parameter. Whereas the algorithm can tolerate
/// a little noise in the normal direction, the normals must be consistently oriented.
/// The smoothness parameter controls the width of the underlying low-pass filter as a factor
/// of the local point spacing. Larger value leads to smoother surfaces and longer computation
/// times. For clean datasets, this value should be set between 1.5 and 2.5. On the other hand,
/// as the amount of noise increases, this value should be increased as well. For these reasons,
/// we do not provide any default value for this parameter.
/// The radius property should correspond to the local point spacing which can be intuitively
/// defined as the average distance to its "natural" one ring neighbors. Currently, this
/// information is only used to define the "surface definition domain" as the union of
/// these balls. Outside this union of balls, the surface is not defined. Therefore, if the balls
/// do not overlap enough, then some holes might appear. If no radius is provided, then they are
/// automatically computed from a basic estimate of the local density based on the 16 nearest
/// neighbors. In the future, this information might be used as well to adjust the width
/// of the low pass filter.
///
/// Note that APSS reconstruction may create small "ghost" connected components
/// close to the reconstructed surface that you should delete.
/// For this purpose, you may call erase_small_polyhedron_connected_components()
/// after make_surface_mesh().
///
/// @heading Is Model for the Concepts:
/// Model of the ImplicitFunction concept.
///
/// @heading Parameters:
/// @param Gt Geometric traits class.

template <class Gt>
class APSS_reconstruction_function
{
// Public types
public:

  typedef Gt Geom_traits; ///< Geometric traits class

  // Geometric types
  typedef typename Geom_traits::FT FT; ///< == Geom_traits::FT
  typedef typename Geom_traits::Point_3 Point; ///< == Geom_traits::Point_3
  typedef typename Geom_traits::Vector_3 Vector; ///< == Geom_traits::Vector_3
  typedef typename Geom_traits::Sphere_3 Sphere; ///< == Geom_traits::Sphere_3

// Private types
private:

  // Item in the Kd-tree: position (Point_3) + normal + index
  class KdTreeElement : public Point_with_normal_3<Gt>
  {
    typedef Point_with_normal_3<Gt> Base; ///< base class

  public:
    unsigned int index;

    KdTreeElement(const Origin& o = ORIGIN, unsigned int id=0)
      : Base(o), index(id)
    {}
    KdTreeElement(const Point& p, const Vector& n = NULL_VECTOR, unsigned int id=0)
      : Base(p, n), index(id)
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

  typedef Fast_orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::Point_ptr_with_transformed_distance
                                    Point_ptr_with_transformed_distance;

// Private initialization method
private:

  /// Creates an APSS implicit function from the [first, beyond) range of points.
  ///
  /// @commentheading Template Parameters:
  /// @param InputIterator iterator over input points.
  /// @param PointPMap is a model of boost::ReadablePropertyMap with a value_type = Geom_traits::Point_3.
  ///        It can be omitted if InputIterator value_type is convertible to Geom_traits::Point_3.
  /// @param NormalPMap is a model of boost::ReadablePropertyMap with a value_type = Geom_traits::Vector_3.
  /// @param RadiusPMap is a model of boost::ReadablePropertyMap with a value_type = FT.
  ///        If it is omitted, a default radius is computed = (distance max to 16 nearest neighbors)/2.
  template <typename InputIterator,
            typename PointPMap,
            typename NormalPMap,
            typename RadiusPMap
  >
  void init(
    InputIterator first,  ///< iterator over the first input point.
    InputIterator beyond, ///< past-the-end iterator.
    PointPMap point_pmap, ///< property map InputIterator -> Point_3.
    NormalPMap normal_pmap, ///< property map InputIterator -> Vector_3.
    RadiusPMap radius_pmap, ///< property map InputIterator -> FT.
    FT smoothness) ///< smoothness factor.
  {
    // Allocate smart pointer to data
    m = new Private;
    m->cached_nearest_neighbor.first = 0;

    int nb_points = std::distance(first, beyond);

    set_smoothness_factor(smoothness);

    // Creates kd-tree
    m->treeElements.reserve(nb_points);
    unsigned int i=0;
    for (InputIterator it=first ; it != beyond ; ++it,++i)
    {
      m->treeElements.push_back(KdTreeElement(get(point_pmap,it), get(normal_pmap,it), i));
    }
    m->tree = new Tree(m->treeElements.begin(), m->treeElements.end());

    m->radii.resize(nb_points);
    if (boost::is_same<RadiusPMap,boost::dummy_property_map>::value)
    {
      // Compute the radius of each point = (distance max to 16 nearest neighbors)/2.
      // The union of these balls defines the surface definition domain.
      int i=0;
      for (InputIterator it=first ; it != beyond ; ++it, ++i)
      {
        Neighbor_search search(*(m->tree), get(point_pmap,it), 16);
        FT maxdist2 = search.begin()->second; // squared distance to furthest neighbor
        m->radii[i] = sqrt(maxdist2)/2.;
      }
    }
    else
    {
      // Copy the radii from given input data
      int i=0;
      for (InputIterator it=first ; it != beyond ; ++it, ++i)
      {
        m->radii[i] = boost::get(radius_pmap,*it);
      }
    }

    // Compute bounding sphere
    Min_sphere_d< CGAL::Optimisation_d_traits_3<Gt> > ms3(first, beyond);
    m->bounding_sphere = Sphere(ms3.center(), ms3.squared_radius());

    // Find a point inside the surface.
    find_inner_point();

    // Dichotomy error when projecting point (squared)
    m->sqError = 1e-7 * Gt().compute_squared_radius_3_object()(m->bounding_sphere);
  }

// Public methods
public:

  /// Creates an APSS implicit function from the [first, beyond) range of points.
  ///
  /// @commentheading Template Parameters:
  /// @param InputIterator iterator over input points.
  /// @param PointPMap is a model of boost::ReadablePropertyMap with a value_type = Geom_traits::Point_3.
  ///        It can be omitted if InputIterator value_type is convertible to Geom_traits::Point_3.
  /// @param NormalPMap is a model of boost::ReadablePropertyMap with a value_type = Geom_traits::Vector_3.
  /// @param RadiusPMap is a model of boost::ReadablePropertyMap with a value_type = FT.
  ///        If it is omitted, a default radius is computed = (distance max to 16 nearest neighbors)/2.

  // This variant requires all parameters.
  template <typename InputIterator,
            typename PointPMap,
            typename NormalPMap,
            typename RadiusPMap
  >
  APSS_reconstruction_function(
    InputIterator first,  ///< iterator over the first input point.
    InputIterator beyond, ///< past-the-end iterator.
    PointPMap point_pmap, ///< property map InputIterator -> Point_3 (access to the position of an input point).
    NormalPMap normal_pmap, ///< property map InputIterator -> Vector_3 (access to the *oriented* normal of an input point).
    RadiusPMap radius_pmap, ///< property map InputIterator -> FT (access to the local point spacing of an input point).
    FT smoothness) ///< smoothness factor. Typical choices are in the range 2 (clean datasets) and 8 (noisy datasets).
  {
    init(
      first, beyond,
      point_pmap,
      normal_pmap,
      radius_pmap,
      smoothness);
  }

  /// @cond SKIP_IN_MANUAL
  // This variant creates a default radius property map = boost::dummy_property_map.
  template <typename InputIterator,
            typename PointPMap,
            typename NormalPMap
  >
  APSS_reconstruction_function(
    InputIterator first,  ///< iterator over the first input point.
    InputIterator beyond, ///< past-the-end iterator over the input points.
    PointPMap point_pmap, ///< property map InputIterator -> Point_3.
    NormalPMap normal_pmap, ///< property map InputIterator -> Vector_3.
    FT smoothness) ///< smoothness factor.
  {
    init(
      first, beyond,
      point_pmap,
      normal_pmap,
      boost::dummy_property_map(),
      smoothness);
  }
  /// @endcond

  /// @cond SKIP_IN_MANUAL
  // This variant creates a default point property map = Dereference_property_map,
  // and a dummy radius property map
  template <typename InputIterator,
            typename NormalPMap
  >
  APSS_reconstruction_function(
    InputIterator first,  ///< iterator over the first input point.
    InputIterator beyond, ///< past-the-end iterator over the input points.
    NormalPMap normal_pmap, ///< property map InputIterator -> Vector_3.
    FT smoothness) ///< smoothness factor.
  {
    init(
      first, beyond,
      make_dereference_property_map(first),
      normal_pmap,
      boost::dummy_property_map(),
      smoothness);
  }
  /// @endcond

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

  /// Sets smoothness factor. Typical choices are in the range 2 (clean datasets) and 8 (noisy datasets).
  void set_smoothness_factor(FT smoothness) { m->nofNeighbors = 6*smoothness*smoothness; }

  /// Returns a sphere bounding the inferred surface.
  const Sphere& bounding_sphere() const
  {
    return m->bounding_sphere;
  }

private:

  /** Fit an algebraic sphere on a set of neigbors in a Moving Least Square sense.
  The weight function is scaled such that the weight of the furthest neighbor is 0.
  */
  void fit(const Neighbor_search& search) const
  {
    FT r2 = search.begin()->second; // squared distance to furthest neighbor

    Vector sumP(0,0,0);
    Vector sumN(0,0,0);
    FT sumDotPP = 0.;
    FT sumDotPN = 0.;
    FT sumW = 0.;

    r2 *= 1.001;
    FT invr2 = 1./r2;
    for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
      Vector p = *(it->first) - CGAL::ORIGIN;
      const Vector& n = it->first->normal();
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
  inline bool isValid(const Point_ptr_with_transformed_distance& nearest_neighbor, const Point& /* p */) const
  {
      FT r = FT(2) * m->radii[nearest_neighbor.first->index];
      return (r*r > nearest_neighbor.second);
  }

public:

  /// 'ImplicitFunction' interface: evaluates implicit function at 3D query point.
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
    if (m->cached_nearest_neighbor.first)
      m->cached_nearest_neighbor.second = sqd(p, *m->cached_nearest_neighbor.first);
    if (!(m->cached_nearest_neighbor.first && isValid(m->cached_nearest_neighbor, p)))
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
        Point pp = *m->cached_nearest_neighbor.first;
        project(pp,n,1);
        // ...and return the (signed) distance to the surface
        Vector h = sub(p,pp);
        return length(h) * ( dot(n,h)>0. ? 1. : -1.);
      }
    }

    // Compute k nearest neighbors and cache the nearest one
    KdTreeElement query(p);
    Neighbor_search search_knn(*(m->tree), query, m->nofNeighbors);
    m->cached_nearest_neighbor = search_nearest(search_knn);

    // If 'p' is close to the surface, fit an algebraic sphere
    // on a set of neigbors in a Moving Least Square sense.
    fit(search_knn);

    // return the distance to the sphere
    return m->as.euclideanDistance(p, *(m->cached_nearest_neighbor.first));
  }

  /// Returns a point located inside the inferred surface.
  Point get_inner_point() const
  {
    return m->inner_point;
  }

// Private methods:
private:

  const Point_ptr_with_transformed_distance& search_nearest(const Neighbor_search& search) const
  {
    typename Neighbor_search::iterator last=search.end(); --last;
    typename Neighbor_search::iterator nearest_it = last;
    for (typename Neighbor_search::iterator it = search.begin(); it != last; ++it)
      if (it->second < nearest_it->second)
        nearest_it = it;
    return *nearest_it;
  }

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

      // neighbors are not sorted anymore,
      // let's find the nearest
      Point_ptr_with_transformed_distance nearest = search_nearest(search);
      // if p is far away the input point cloud, start with the closest point.
      if (!isValid(nearest,p))
      {
        p = *nearest.first;
        n =  nearest.first->normal();
        delta2 = nearest.second;
      }
      else
      {
        fit(search);

        Point oldP = p;
        m->as.project(p,n,*(nearest.first));

        if (!isValid(nearest,p))
        {
          std::cout << "Invalid projection\n";
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

    /** Projects a point onto the surface of the sphere.
      * This function considers the projection which is the closest
      * to a given reference point.
      */
    void project(Point& p, Vector& n, const Point& reference_point) {
      if (state==AlgebraicSphere::SPHERE)
      {
        Vector dir = sub(p,center);
        FT l = length(dir);
        dir = dir * (1./l);
        p = add(center,mul( radius,dir));
        FT flip = u4<0. ? -1. : 1.;
        n = mul(flip,dir);

        Point other = add(center,mul(-radius,dir));
        typename Geom_traits::Compute_squared_distance_3 sqd;
        if (sqd(reference_point,other) < 0.9*sqd(reference_point,p))
        {
          p = other;
          n = -n;
        }
      }
      else if (state==AlgebraicSphere::PLANE)
      {
        // projection onto a plane.
        p = sub(p, mul(dot(normal,p-CGAL::ORIGIN) + d, normal));
        n = normal;
      }
      else
      {
        // iterative projection onto the algebraic sphere.
        p = iProject(p);
      }
    }

    /** Compute the Euclidean distance between the algebraic surface and a point 'p'.
      * This function uses the closest intersection to a given reference point.
      */
    FT euclideanDistance(const Point& p, const Point& reference_point) {
      if (state==SPHERE)
      {
        // tricky case because we have to pick to best intersection
        // that is the closest to the reference point
        Vector dir = sub(p,center);
        FT l = length(dir);
        FT inside = l < radius ? -1. : 1;
        dir = dir * (1./l);
        Point proj0 = add(center,mul( radius,dir));
        Point proj1 = add(center,mul(-radius,dir));
        FT flip = u4<0. ? -1. : 1.;

        typename Geom_traits::Compute_squared_distance_3 sqd;
        if (sqd(reference_point,proj1) < 0.9*sqd(reference_point,proj0))
        {
          proj0 = proj1;
          inside = -1;
        }

        return length(sub(p,proj0)) * flip * inside;
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
      // Creates random point in bounding sphere
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
    Sphere bounding_sphere; // Points' bounding sphere
    FT sqError; // Dichotomy error when projecting point (squared)
    unsigned int nofNeighbors; // Number of nearest neighbors
    Point inner_point; // Point inside the surface
    mutable AlgebraicSphere as;
    mutable Point_ptr_with_transformed_distance cached_nearest_neighbor;
    int count; // reference counter
  };

  Private* m; // smart pointer to data

}; // end of APSS_reconstruction_function


CGAL_END_NAMESPACE

#endif // CGAL_APSS_RECONSTRUCTION_FUNCTION_H
