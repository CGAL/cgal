// Author: Laurent Saboret

#ifndef UI_POINT_3_H
#define UI_POINT_3_H

#include <Gyroviz_point_3.h>
#include <CGAL/Orientable_normal_3.h>
#include <CGAL/Iterator_project.h>

#include <set>
#include <algorithm>


/// The UI_point_3 class represents a 3D point in Surface_reconstruction_points_3 demo. 
/// It contains:
/// - a position,
/// - a normal (oriented or not),
/// - an original normal (optional, always oriented),
/// - a list of camera/2D point pairs used to reconstruct the point from images,
/// - a selection flag.
///
/// @heading Is Model for the Concepts: 
/// Model of the PointWithOrientableNormal_3 concept.
///
/// @heading Parameters:
/// @param Gt   Kernel's geometric traits.

template<class Gt>
class UI_point_3 
  : public Gyroviz_point_3<Gt>
{
// Private types
private:

    // Base class
    typedef Gyroviz_point_3<Gt> Base;

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
    /// Camera list is empty by default.
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
    template <class K>
    UI_point_3(const Gyroviz_point_3<K>& gpt)
    : Base(gpt)
    {
      m_is_selected = false;
    }
    template < class InputIterator >
    UI_point_3(const Point_3& point,
               InputIterator first_camera_point2_pair, 
               InputIterator beyond_camera_point2_pair)
    : Base(point, first_camera_point2_pair, beyond_camera_point2_pair)
    {
      m_is_selected = false;
    }
    template < class InputIterator >
    UI_point_3(const Point_3& point,
               const Normal& normal,
               InputIterator first_camera_point2_pair, 
               InputIterator beyond_camera_point2_pair)
    : Base(point, normal, first_camera_point2_pair, beyond_camera_point2_pair)
    {
      m_is_selected = false;
    }

    /// Copy constructor
    UI_point_3(const UI_point_3& upt)
    : Base(upt)
    {
      m_is_selected = upt.m_is_selected;
      m_original_normal = upt.m_original_normal;
    }
    template<class K>
    UI_point_3(const UI_point_3<K>& upt)
    : Base(upt)
    {
      m_is_selected = upt.is_selected();
      m_original_normal = upt.m_original_normal;
    }
    /// Operator =()
    UI_point_3& operator=(const UI_point_3& upt)
    {
      Base::operator=(upt);
      m_is_selected = upt.m_is_selected;
      m_original_normal = upt.m_original_normal;
      return *this;
    }

    /// Merge points, including lists of camera/2D point pairs.
    void merge(const UI_point_3& upt)
    { 
      Base::merge(upt); 
    }

    // Inherited operators ==() and !=() are fine.

    /// Selection flag.
    bool is_selected() const { return m_is_selected; }
    void select(bool is_selected=true) { m_is_selected = is_selected; }

    /// Get/set *original* normal.
    const Vector_3& original_normal() const { return m_original_normal; }
    Vector_3&       original_normal()       { return m_original_normal; }

// Data
private:

    // Selection flag.
    bool m_is_selected;

    /// *Original* normal.
    Vector_3  m_original_normal;
};


#endif //UI_POINT_3_H

