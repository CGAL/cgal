#ifndef TRIANGULAR_SURFACE_3_H
#define TRIANGULAR_SURFACE_3_H

#include <algorithm>
#include <vector>
#include <deque>
#include <GL/gl.h>


/// The Triangular_surface_3 class is array of triangles + normals 
/// of type std::pair<Triangle_3,Vector_3>. 
/// It provides:
/// - normals computation
/// - OpenGL rendering
/// 
/// @heading Parameters:
/// @param Gt       Geometric traits class.

template <class Gt>
class Triangular_surface_3 
  : public std::vector<std::pair<typename Gt::Triangle_3, typename Gt::Vector_3> >
{
// Private types
private:

  // Base class 
  typedef std::vector<std::pair<typename Gt::Triangle_3, typename Gt::Vector_3> > 
                                  Base;

// Public types
public:

  // Repeat base class' types
  /// @cond SKIP_IN_MANUAL
  using Base::iterator;
  using Base::const_iterator;
  /// @endcond

  // Geometric types
  typedef Gt  Geom_traits; ///<Geometric traits class.
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;
  typedef typename Geom_traits::Vector_3 Vector;
  typedef typename Geom_traits::Triangle_3 Triangle;

  // Type of the container's elements
  typedef typename std::pair<Triangle,Vector> Facet_with_normal;

// Public methods
public:

  // Default constructor, copy constructor and operator =() are fine

  // Repeat base class' public methods used below
  /// @cond SKIP_IN_MANUAL
  Base::begin;
  Base::end;
  Base::size;
  Base::push_back;
  /// @endcond
  
  /// Compute triangle normal and append it.
  void push_back(const Triangle& t)
  {
      insert(end(), t);
  }

  /// Compute triangle normal and insert it.
  iterator insert(iterator loc, const Triangle& t)
  {
      // Compute normal
      Vector u = t[1] - t[0];
      Vector v = t[2] - t[0];
      Vector n = CGAL::cross_product(u,v);
      n = n / std::sqrt(n*n);
      
      return Base::insert(loc, Facet_with_normal(t,n));
  }

  /// Insert range of triangles.
  ///
  /// Precondition: the value type of InputIterator must 'Triangle'.
  ///
  /// @param first First triangle to add.
  /// @param beyond Past-the-end triangle to add.
  template < class InputIterator >
  void insert(iterator loc, 
              InputIterator first, InputIterator beyond)
  {
    // Compute normals in temporary container
    std::deque<Facet_with_normal> tmp;
    for (InputIterator it = first; it != beyond; ++it)
    {
      const Triangle& t = *it;
      Vector u = t[1] - t[0];
      Vector v = t[2] - t[0];
      Vector n = CGAL::cross_product(u,v);
      n = n / std::sqrt(n*n);
      
      tmp.push_back(Facet_with_normal(t,n));
    }

    Base::insert(loc, tmp.begin(), tmp.end());
  }

  // Draw triangles using OpenGL calls.
  void gl_draw_surface() const
  {
    ::glBegin(GL_TRIANGLES);

    for(const_iterator it = begin(); it != end(); it++)
    {
      const Facet_with_normal& f = *it;

      const Vector& n = f.second;
      ::glNormal3d(n.x(),n.y(),n.z());

      const Triangle& t = f.first;
      const Point& a = t[0];
      const Point& b = t[1];
      const Point& c = t[2];
      ::glVertex3d(a.x(),a.y(),a.z());
      ::glVertex3d(b.x(),b.y(),b.z());
      ::glVertex3d(c.x(),c.y(),c.z());
    }
    
    ::glEnd();
  }

}; // end of class Triangular_surface_3


#endif // TRIANGULAR_SURFACE_3_H
