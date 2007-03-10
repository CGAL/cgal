#ifndef _DT3_
#define _DT3_

#include <CGAL/Delaunay_triangulation_3.h>

template < class kernel, class Tds >
class DT3 : public CGAL::Delaunay_triangulation_3<kernel, Tds>
{
public:
	typedef typename kernel::FT FT;
	typedef typename kernel::Point_3 Point;
	typedef typename kernel::Segment_3 Segment;

public:

public:
	DT3()
	{
	}
	~DT3()
	{
	}

public:

	// draw Delaunay vertices
	void gl_draw_delaunay_vertices(const unsigned char r,
		                  const unsigned char g,
											const unsigned char b,
											const float size)
	{
		::glPointSize(size);
		::glColor3ub(r,g,b);
		::glBegin(GL_POINTS);

		// Voronoi vertices are dual of Voronoi cells
		Point_iterator it;
		for(it = points_begin();
			  it != points_end();
				it++)
		{
			const Point& p = *it;
			::glVertex3d(p.x(),p.y(),p.z());
		}
		::glEnd();
	}


	// draw Delaunay edges
	void gl_draw_voronoi_vertices(const unsigned char r,
		                            const unsigned char g,
																const unsigned char b,
																const float size)
	{
		::glPointSize(size);
		::glColor3ub(r,g,b);
		::glBegin(GL_POINTS);

		// Voronoi vertices are dual of Voronoi cells
		Finite_cells_iterator it;
		for(it = finite_cells_begin();
			  it != finite_cells_end();
				it++)
		{
			Point p = dual(it);
			::glVertex3d(p.x(),p.y(),p.z());
		}
		::glEnd();
	}

	// draw Delaunay edges
	void gl_draw_delaunay_edges(const unsigned char r,
		                          const unsigned char g,
															const unsigned char b,
															const float width)
	{
		::glLineWidth(width);
		::glColor3ub(r,g,b);
		::glBegin(GL_LINES);
		Finite_edges_iterator it;
		for(it = finite_edges_begin();
			  it != finite_edges_end();
				it++)
		{
			Segment s = segment(*it);
			Point p1 = s.source();
			Point p2 = s.target();
			::glVertex3d(p1.x(),p1.y(),p1.z());
			::glVertex3d(p2.x(),p2.y(),p2.z());
		}
		::glEnd();
	}

	// draw Voronoi edges
	void gl_draw_voronoi_edges(const unsigned char r,
			                        const unsigned char g,
															const unsigned char b,
															const float width)
	{
		::glLineWidth(width);
		::glColor3ub(r,g,b);

		::glBegin(GL_LINES);
		Finite_facets_iterator it;
		for(it = finite_facets_begin();
				it != finite_facets_end();
				it++)
		{
			Object object = dual(*it);
			Segment segment;
			if(assign(segment,object))
			{
				Point p1 = segment.source();
				Point p2 = segment.target();
				::glVertex3d(CGAL_NTS to_double(p1.x()),
					            CGAL_NTS to_double(p1.y()),
					            CGAL_NTS to_double(p1.z()));
				::glVertex3d(CGAL_NTS to_double(p2.x()),
						          CGAL_NTS to_double(p2.y()),
											CGAL_NTS to_double(p2.z()));
			}
			Ray ray;
			if(assign(ray,object))
			{
				Point p1 = ray.source();
				Point p2 = ray.point(1);
				::glVertex3d(CGAL_NTS to_double(p1.x()),
						          CGAL_NTS to_double(p1.y()),
						          CGAL_NTS to_double(p1.z()));
				::glVertex3d(CGAL_NTS to_double(p2.x()),
						          CGAL_NTS to_double(p2.y()),
											CGAL_NTS to_double(p2.z()));
			}
		}
		::glEnd();
	}
};

#endif // _DT3_
