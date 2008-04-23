#ifndef Gyroviz_point_dt3_H_INCLUDED
#define Gyroviz_point_dt3_H_INCLUDED

#include <CGAL/Point_3.h>

#include <vector>
#include <algorithm>


/// The Gyroviz_point_dt3 class represents a 3D point with:
/// - a position,
///  - a list of cameras used to reconstruct the point from an image sequence.
///
/// @heading Is Model for the Concepts: Model of the Point_3 concept.
///
/// @heading Parameters:
/// @param Gt   Kernel's geometric traits.

template<class Gt>
class Gyroviz_point_dt3 : public Gt::Point_3 
{
	// Private types
private:

	typedef typename Gt::Point_3 Base;

	// Public types
public:

	// Repeat Point_3 public types
	typedef Gt Geom_traits; ///< Kernel's geometric traits
	typedef typename Geom_traits::FT FT;
	typedef typename Geom_traits::Point_3  Point;  ///< Kernel's Point_3 class.

	// Public methods
public:

	/// Point is (0,0,0) by default.
	/// Camera list is empty by default.
	Gyroviz_point_dt3(const CGAL::Origin& o = CGAL::ORIGIN)
		: Base(o)
	{
	}
	Gyroviz_point_dt3(FT x, FT y, FT z)
		: Base(x,y,z)
	{
	}
	Gyroviz_point_dt3(const Point& point)
		: Base(point)
	{
	}
	template < class InputIterator >
	Gyroviz_point_dt3(const Point& point,
		InputIterator first_camera, InputIterator beyond_camera)
		: Base(point)
	{
		std::copy(first_camera, beyond_camera, std::back_inserter(list_of_cameras));
	}

	// Default copy constructor and operator =() are fine

	/// Compare positions
	bool operator==(const Gyroviz_point_dt3& that)
	{ 
		return Base::operator==(that); 
	}
	bool operator!=(const Gyroviz_point_dt3& that)
	{ 
		return ! (*this == that); 
	}

	// Get/set cameras
	template < class InputIterator >
	void set_cameras(InputIterator first_camera, InputIterator beyond_camera)
	{
		list_of_cameras.clear();
		std::copy(first_camera, beyond_camera, std::back_inserter(list_of_cameras));
	}
	const std::vector<Point>& cameras() const { return	list_of_cameras; }

	// Data
private:

	// List of cameras
	std::vector<Point> list_of_cameras;
};


#endif //Gyroviz_point_dt3_H_INCLUDED

