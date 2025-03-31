
#pragma once

#include <array>

#include <CGAL/Filtered_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Dimension.h>

namespace CGAL {
namespace Test {

struct Vec3
{
  Vec3()
  {}

  Vec3(double x, double y, double z)
   : x(x), y(y), z(z)
    {}

  Vec3(double x, double y, double z, double w)
   : x(x/w), y(y/w), z(z/w)
    {}
  const double& hw() const
  {
      return CGAL::constant<double, 1>();
  }
  void foo() const
  {
    std::cout << "foo" << std::endl;
  }
    double x, y, z;
};

std::ostream& operator<<(std::ostream& os, const Vec3& v)
{
  return os << v.x << " " << v.y << " " << v.z << std::endl;
}

/// Functor that constructs a CGAL 3D bounding box given a 3D point (vector).
template<class BaseConstructBbox>
class ConstructBbox : public BaseConstructBbox
{
  public:
    using BaseConstructBbox::operator();

    /// Returns a bounding box that contains the given vector (point) \a v.
    ///
    /// \param v point to construct a bounding box for
    /// \return the bounding box that contains the given point \a v
    CGAL::Bbox_3 operator()(const Vec3 &v) const
    {
      return CGAL::Bbox_3(v.x, v.y, v.z, v.x, v.y, v.z);
    }
};

/// Functor that implements various interfaces for constructing a 3D point.
template<typename CustomKernel, typename CgalKernelBase>
class ConstructPoint
{
  private:
    typedef typename CustomKernel::RT RT;
    typedef typename CustomKernel::Point_3 Point_3;
    typedef typename CustomKernel::Weighted_point_3 Weighted_point_3;
    typedef typename Point_3::Rep Rep;

  public:

    template<typename>
    struct result {
      typedef Point_3 type;
    };

    template<typename F>
    struct result<F(Weighted_point_3)> {
      typedef const Point_3& type;
    };

    template<typename F>
    struct result<F(Point_3)> {
      typedef const Point_3& type;
    };

    // Note : the CGAL::Return_base_tag is really internal CGAL stuff.
    // Unfortunately it is needed for optimizing away copy-constructions,
    // due to current lack of delegating constructors in the C++ standard.
    /// \cond
    Rep operator()(CGAL::Return_base_tag, CGAL::Origin o) const { return Rep(0,0,0); }
    Rep operator()(CGAL::Return_base_tag, const RT &x, const RT &y, const RT &z) const { return Rep(x, y, z); }
    Rep operator()(CGAL::Return_base_tag, const RT &x, const RT &y, const RT &z, const RT &w) const { return Rep(x, y, z, w); }
    /// \endcond

    /// Returns a 3D point representing the origin for the parameterized kernel.
    ///
    /// \return a 3D point representing the origin for the parameterized kernel
    Point_3 operator()(const CGAL::Origin&) const { return Vec3(0, 0, 0); }

    const Point_3& operator()(const Point_3& p) const { return p;}

    const Point_3& operator()(const Weighted_point_3& wp) const {  return CgalKernelBase().construct_point_3_object()(wp); }
    /// Returns a 3D point representing the cartesian (\a x, \a y, \a z)
    /// coordinate.
    ///
    /// \param x X-coordinate
    /// \param y Y-coordinate
    /// \param z Z-coordinate
    /// \return a 3D point representing the cartesian (\a x, \a y, \a z)
    ///   coordinate
    Point_3 operator()(const RT &x, const RT &y, const RT &z) const { return Vec3(x, y, z); }

    Point_3 operator()(const RT &x, const RT &y, const RT &z, const RT &w) const { return Vec3(x/w, y/w, z/w); }

    //const Point_3 &operator()(const Point_3 &p) const { return p; }
};

/// Functor that returns an iterator to iterate over the x-, y- and z-coordinate
/// of a 3D point.
template <class  K>
class ConstructCoordIterator
{
  public:
    typedef const double* result_type;

    /// Returns the begin iterator for the given point \a v that iterates over
    /// the x-, y- and z-coordinates of the point.
    ///
    /// \param v point to iterate over
    /// \return begin iterator that iterates the coordinates of the given point
    ///   \a v
    const double* operator()(const Vec3 &v) const { return &v.x; }

    /// Returns the iterator for the \a n th coordinate of \a v.
    ///
    /// \param n offset from the x-coordinate to start iteration, for example,
    ///   for \a n = 1, this returns the  iterator pointing to the y-coordinate
    /// \param v point to iterate over
    /// \return iterator that iterates the coordinates of the given point
    ///   \a v starting at the given offset \a n
    const double* operator()(const Vec3 &v, int n) const
    {
      const double* pyptr = &v.x;
      return pyptr + n;
    }

     const double* operator()(const typename K::Vector_3& v)
     {
      return typename K::Construct_cartesian_const_iterator_3()(v);
     }
     const double* operator()(const typename K::Vector_3& v, int n)
     {
      return typename K::Construct_cartesian_const_iterator_3()(v, n);
     }
};

/// Functor that returns the x-coordinate of some 3D point.
template <typename CustomKernel, typename CgalKernelBase>
struct ComputeX
  : public CgalKernelBase::Compute_x_3
{
  using CgalKernelBase::Compute_x_3::operator();

  typedef const double& result_type;

  /// Returns the x-coordinate of the given point \a v.
  ///
  /// \param v point to return the x-coordinate for
  /// \return the x-coordinate of the given point \a v
  result_type operator()(const typename CustomKernel::Point_3&v) const { return v.rep().x; }

  // result_type operator()(const typename K::Vector_3 &v) const { return v.x(); } // I  hoped to just inherit
};

/// Functor that returns the y-coordinate of some 3D point.
template <typename CustomKernel, typename CgalKernelBase>
struct ComputeY
  : public CgalKernelBase::Compute_y_3
{
  using CgalKernelBase::Compute_y_3::operator();

  typedef const double& result_type;

  /// Returns the y-coordinate of the given point \a v.
  ///
  /// \param v point to return the y-coordinate for
  /// \return the y-coordinate of the given point \a v
  result_type operator()(const typename CustomKernel::Point_3&v) const { return v.rep().y; }
};

/// Functor that returns the z-coordinate of some 3D point.
template <typename CustomKernel, typename CgalKernelBase>
struct ComputeZ
  : public CgalKernelBase::Compute_z_3
{
  using CgalKernelBase::Compute_z_3::operator();

  typedef const double& result_type;

  /// Returns the z-coordinate of the given point \a v.
  ///
  /// \param v point to return the z-coordinate for
  /// \return the z-coordinate of the given point \a v
  result_type operator()(const typename CustomKernel::Point_3&v) const { return v.rep().z; }
};

template <typename CustomKernel, typename CgalKernelBase>
struct ComputeW : public CgalKernelBase::Compute_hw_3
{
  using CgalKernelBase::Compute_hw_3::operator();

  typedef const double& result_type;

  /// Returns the z-coordinate of the given point v.
  ///
  /// \param v point to return the homogeneous coordinate for
  /// \return the homogeneous coordinate of the given point v
    result_type operator()(const typename CustomKernel::Point_3&v) const { return v.rep().hw(); }
};

/// CGAL CRTP kernel type that enables derivation of some custom kernel type
/// from a predefined kernel type. This allows the usage of
/// `CGAL::Type_equality_wrapper` in the definition of the custom kernel type.
/// It defines a bunch of type aliases that point to the right types for
/// constructing 3D points, calculating coordinates for a point, etc.
template<typename CustomKernel, typename CgalKernel>
class CartesianCrtp
  : public CgalKernel::template Base<CustomKernel>::Type
{
  private:
    typedef typename CgalKernel::template Base<CustomKernel>::Type   CgalKernelBase;

public:
    typedef CustomKernel                                             Kernel;
    typedef Vec3                                                     Point_3;
    typedef ConstructPoint<Kernel, CgalKernelBase>                   Construct_point_3;

    typedef const double*                                            Cartesian_const_iterator_3;
    typedef ConstructCoordIterator<CgalKernelBase>                   Construct_cartesian_const_iterator_3;

    typedef ComputeX<CustomKernel, CgalKernelBase>                   Compute_x_3;
    typedef ComputeY<CustomKernel, CgalKernelBase>                   Compute_y_3;
    typedef ComputeZ<CustomKernel, CgalKernelBase>                   Compute_z_3;

    typedef Compute_x_3                                              Compute_hx_3;
    typedef Compute_y_3                                              Compute_hy_3;
    typedef Compute_z_3                                              Compute_hz_3;
    typedef ComputeW<CustomKernel, CgalKernelBase>                   Compute_hw_3;

    typedef ConstructBbox<typename CgalKernelBase::Construct_bbox_3> Construct_bbox_3;

    Construct_point_3 construct_point_3_object() const { return Construct_point_3(); }
    Construct_bbox_3 construct_bbox_3_object() const { return Construct_bbox_3(); }
    Construct_cartesian_const_iterator_3 construct_cartesian_const_iterator_3_object() const { return Construct_cartesian_const_iterator_3(); }

    Compute_x_3 compute_x_3_object() const { return Compute_x_3(); }
    Compute_y_3 compute_y_3_object() const { return Compute_y_3(); }
    Compute_z_3 compute_z_3_object() const { return Compute_z_3(); }
    Compute_hx_3 compute_hx_3_object() const { return Compute_hx_3(); }
    Compute_hy_3 compute_hy_3_object() const { return Compute_hy_3(); }
    Compute_hz_3 compute_hz_3_object() const { return Compute_hz_3(); }
    Compute_hw_3 compute_hw_3_object() const { return Compute_hw_3(); }

    template<typename KernelT>
    struct Base
    {
      typedef CartesianCrtp<KernelT, CgalKernel> Type;
    };
};

struct CartesianKernel
  : public CGAL::Type_equality_wrapper<CartesianCrtp<CartesianKernel,
                                                     CGAL::Simple_cartesian<double>>, CartesianKernel>
{
};

/// Defines a filtered kernel that allows for exact predicates. This is required
/// in some cases when performing mesh intersections for example.
using CustomKernel = CGAL::Filtered_kernel<CartesianKernel>;

} // namespace Test

template <>
struct Ambient_dimension<::CGAL::Test::Vec3, ::CGAL::Test::CustomKernel>{
  static const int value = 3;
  typedef Dimension_tag<3> type;
};
} // namespace CGAL
