// Author: Laurent Saboret

#ifndef UI_POINT_3_H
#define UI_POINT_3_H

#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Iterator_project.h>

#include <set>
#include <algorithm>


/// The UI_point_3 class represents a 3D point in Surface_reconstruction_points_3 demo.
/// It contains:
/// - a position,
/// - a normal,
/// - a radius,
/// - a selection flag.
///
/// @heading Parameters:
/// @param Gt   Geometric traits class.

template<class Gt>
class UI_point_3
  : public CGAL::Point_with_normal_3<Gt>
{
// Private types
private:

    // Base class
    typedef CGAL::Point_with_normal_3<Gt> Base;

// Public types
public:

    /// Base class
    typedef Base Point_with_normal; 

    // Repeat base class public types
    typedef Gt Geom_traits; ///< Geometric traits class.
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::RT RT;
    typedef typename Geom_traits::Point_2  Point_2;  ///< typedef to Geom_traits::Point_2
    typedef typename Geom_traits::Point_3  Point_3;  ///< typedef to Geom_traits::Point_3
    typedef typename Geom_traits::Vector_3 Vector_3; ///< typedef to Geom_traits::Vector_3

// Public methods
public:

    /// Point is (0,0,0) by default.
    /// Normal is (0,0,0) by default.
    UI_point_3(const CGAL::Origin& o = CGAL::ORIGIN)
    : Base(o)
    {
      m_is_selected = false;
      m_radius = FT(0);
    }
    UI_point_3(FT x, FT y, FT z,
               const Vector_3& normal = CGAL::NULL_VECTOR)
    : Base(x,y,z,normal)
    {
      m_is_selected = false;
      m_radius = FT(0);
    }
    UI_point_3(RT hx, RT hy, RT hz, RT hw,
               const Vector_3& normal = CGAL::NULL_VECTOR)
    : Base(hx,hy,hz,hw,normal)
    {
      m_is_selected = false;
      m_radius = FT(0);
    }
    UI_point_3(const Point_3& point,
               const Vector_3& normal = CGAL::NULL_VECTOR)
    : Base(point, normal)
    {
      m_is_selected = false;
      m_radius = FT(0);
    }
    template <class K>
    UI_point_3(const CGAL::Point_with_normal_3<K>& pwn)
    : Base(pwn)
    {
      m_is_selected = false;
      m_radius = FT(0);
    }

    /// Copy constructor
    UI_point_3(const UI_point_3& upt)
    : Base(upt)
    {
      m_is_selected = upt.m_is_selected;
      m_radius = upt.m_radius;
    }
    template<class K>
    UI_point_3(const UI_point_3<K>& upt)
    : Base(upt)
    {
      m_is_selected = upt.is_selected();
      m_radius = upt.radius();
    }
    /// Operator =()
    UI_point_3& operator=(const UI_point_3& upt)
    {
      Base::operator=(upt);
      m_is_selected = upt.m_is_selected;
      m_radius = upt.m_radius;
      return *this;
    }

    // Inherited operators ==() and !=() are fine.

    /// Selection flag.
    bool is_selected() const { return m_is_selected; }
    void select(bool is_selected=true) { m_is_selected = is_selected; }

    /// Gets/sets radius.
    FT radius() const { return m_radius; }
    FT& radius() { return m_radius; }

// Data
private:

    // Selection flag.
    bool m_is_selected;

    /// radius.
    FT m_radius;
};

#include <CGAL/Kernel_traits.h>

namespace CGAL{
template<class Gt>
struct Kernel_traits< ::UI_point_3<Gt> >{
  typedef Gt Kernel;
};
} //end of CGAL namespace

#endif //UI_POINT_3_H

