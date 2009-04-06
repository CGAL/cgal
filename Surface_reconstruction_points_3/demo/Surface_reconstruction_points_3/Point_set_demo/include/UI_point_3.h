// Author: Laurent Saboret

#ifndef UI_POINT_3_H
#define UI_POINT_3_H

#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Orientable_normal_3.h>
#include <CGAL/Iterator_project.h>

#include <set>
#include <algorithm>


/// The UI_point_3 class represents a 3D point in Surface_reconstruction_points_3 demo.
/// It contains:
/// - a position,
/// - a normal (oriented or not),
/// - a selection flag.
///
/// @heading Is Model for the Concepts:
/// Model of the PointWithOrientableNormal_3 concept.
///
/// @heading Parameters:
/// @param Gt   Kernel's geometric traits.

template<class Gt>
class UI_point_3
  : public CGAL::Point_with_normal_3<Gt, CGAL::Orientable_normal_3<Gt> >
{
// Private types
private:

    // Base class
    typedef CGAL::Point_with_normal_3<Gt, CGAL::Orientable_normal_3<Gt> > Base;

// Public types
public:

    // Base class
    typedef Base Point_with_normal; ///< Model of the PointWithOrientableNormal_3 concept.

    // Repeat base class public types
    typedef Gt Geom_traits; ///< Kernel's geometric traits.
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::RT RT;
    typedef typename Geom_traits::Point_2 Point_2; ///< Kernel's Point_2 class.
    typedef typename Geom_traits::Point_3 Point_3; ///< Kernel's Point_3 class.
    typedef typename Geom_traits::Vector_3 Vector_3; ///< Kernel's Vector_3 class.
    typedef typename Point_with_normal::Normal Normal; ///< Model of OrientableNormal_3 concept.

// Public methods
public:

    /// Point is (0,0,0) by default.
    /// Normal is (0,0,0) and is oriented by default.
    UI_point_3(const CGAL::Origin& o = CGAL::ORIGIN)
    : Base(o)
    {
      m_is_selected = false;
    }
    UI_point_3(FT x, FT y, FT z,
               const Normal& normal = CGAL::NULL_VECTOR)
    : Base(x,y,z,normal)
    {
      m_is_selected = false;
    }
    UI_point_3(RT hx, RT hy, RT hz, RT hw,
               const Normal& normal = CGAL::NULL_VECTOR)
    : Base(hx,hy,hz,hw,normal)
    {
      m_is_selected = false;
    }
    UI_point_3(const Point_3& point,
               const Normal& normal = CGAL::NULL_VECTOR)
    : Base(point, normal)
    {
      m_is_selected = false;
    }
    template <class K, class N>
    UI_point_3(const CGAL::Point_with_normal_3<K,N>& pwn)
    : Base(pwn)
    {
      m_is_selected = false;
    }

    /// Copy constructor
    UI_point_3(const UI_point_3& upt)
    : Base(upt)
    {
      m_is_selected = upt.m_is_selected;
    }
    template<class K>
    UI_point_3(const UI_point_3<K>& upt)
    : Base(upt)
    {
      m_is_selected = upt.is_selected();
    }
    /// Operator =()
    UI_point_3& operator=(const UI_point_3& upt)
    {
      Base::operator=(upt);
      m_is_selected = upt.m_is_selected;
      return *this;
    }

    // Inherited operators ==() and !=() are fine.

    /// Selection flag.
    bool is_selected() const { return m_is_selected; }
    void select(bool is_selected=true) { m_is_selected = is_selected; }

// Data
private:

    // Selection flag.
    bool m_is_selected;
};


#endif //UI_POINT_3_H

