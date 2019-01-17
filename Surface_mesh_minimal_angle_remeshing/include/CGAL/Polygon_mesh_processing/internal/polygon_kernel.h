#ifndef POLYGON_KERNEL_H
#define POLYGON_KERNEL_H

#include <CGAL/convex_hull_2.h>

#include "console_color.h"

template <class Kernel>
class Polygon_kernel
{
private:
	typedef typename Kernel::FT FT;
	typedef typename Kernel::Line_2     Line;
	typedef typename Kernel::Point_2    Point;
	typedef typename Kernel::Vector_2   Vector;
	typedef typename Kernel::Segment_2  Segment;

public:
	Polygon_kernel() {}
	~Polygon_kernel() {}

public:
	template <class InputIterator, // 2D segments
	          class OutputIterator> // 2D points
	bool run(InputIterator begin,
	InputIterator end,
	const Point& inside_point,
	OutputIterator out)
	{
		// compute dual of translated polygon w.r.t. inside point.
		Vector translate = inside_point - CGAL::ORIGIN;
		std::list<Point> dual_points;
		InputIterator it;
		for(it = begin; it != end; it++)
		{
		  const Segment& segment = *it; // no care about orientation?
		  Point dual_point = dual(segment, -translate);
		  dual_points.push_back(dual_point);
		}

		kernel_from_dual_points(dual_points, translate, out);
		return true;
	}

	template <class OutputIterator> 
	void kernel_from_dual_points(std::list<Point>& dual_points,
		                           const Vector& translate,
															 OutputIterator out)
	{
    // compute convex hull in dual space (result is ccw)
    std::vector<Point> hull;
    CGAL::convex_hull_2(dual_points.begin(), 
			                  dual_points.end(), 
												std::back_inserter(hull));

    // dualize and translate back 
    for(unsigned int i=0;i<hull.size();i++)
    {
      const Point& a = hull[i];
      const Point& b = hull[(i+1)%hull.size()];
      Segment ab(a,b);
			Point point = dual(ab, CGAL::NULL_VECTOR) + translate;
			*out++ = point;
    }
	}

  // translate segment, then compute dual
  Point dual(const Segment& segment,
             Vector translate)
  {
    const Point a = segment[0] + translate;
    const Point b = segment[1] + translate;
    Vector normal = unit_normal(b-a);
    Line line(a,b); 

		// flip normal if required
		// project origin
		Point projection = line.projection(CGAL::ORIGIN);
		if(normal * (projection - CGAL::ORIGIN) < 0.0)
			normal = -normal;

    FT d_origin = distance_to_origin(line);
		if(d_origin == 0.0)
		{
			// FIXME
			// std::cerr << red << "dual: null distance to origin" << white << std::endl;
			const FT a_lot = 1e10; // std::numeric_limits<FT>::max();
			CGAL::ORIGIN + normal * a_lot;
			return CGAL::ORIGIN; 
		}
		else
			return CGAL::ORIGIN + normal / d_origin;
  }

  Vector normalize(const Vector& v)
  {
		if(v*v != 0.0)
			return v / std::sqrt(v*v);
		else
		{
			std::cerr << "zero norm in normalize(vector)" << std::endl;
			return v;
		}
  }

  Vector unit_normal(const Vector& v)
  {
    Vector normal(v.y(), -v.x()); // 90 deg CW rotation 
    return normalize(normal);
  }

  double distance_to_origin(const Line& line)
  {
    Point origin = CGAL::ORIGIN;
    return std::sqrt(CGAL::squared_distance(origin, line));
  }
};

#endif // POLYGON_KERNEL_H
