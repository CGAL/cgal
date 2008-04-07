// Author     : Nader Salman

#ifndef _Gyroviz_border_points_dt2_
#define _Gyroviz_border_points_dt2_

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <CGAL/Delaunay_triangulation_2.h>

template < class Gt, class Tds >
class Gyroviz_border_points_dt2 : public CGAL::Delaunay_triangulation_2<Gt, Tds>
{
	// Private types
private:

	typedef CGAL::Delaunay_triangulation_2<Gt, Tds>  Base;

	// Public types
public:

	// Repeat Delaunay_triangulation_2 public types
	typedef Tds Triangulation_data_structure;
	typedef Gt  Geom_traits;
	typedef typename Geom_traits::FT FT;
	typedef typename Geom_traits::Point_2        Point_2;
	typedef typename Geom_traits::Segment_2      Segment_2;
	typedef typename Base::Face_handle           Face_handle;
	typedef typename Base::Vertex_handle         Vertex_handle;
	typedef typename Base::Edge                  Edge;
	typedef typename Base::Finite_edges_iterator Finite_edges_iterator;
	typedef typename Base::Finite_faces_iterator Finite_faces_iterator;
	typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;

	// Data members
private:
	std::vector<Vertex_handle> input;

	// Public methods
public:


	// Constructors
	Gyroviz_border_points_dt2(){}
	Gyroviz_border_points_dt2(std::vector<Vertex_handle> in)
	{
		for(int i=0; i<in.size(); ++i)
		{
			this->insert(in[i]->point());
		}
	}


	// Functions

	// this function will store as a vector of Gyroviz_vertex_segment_2 
	// the entire 2D Delaunay triangulation
	std::vector<Segment_2> segments_out_of_dt2()
	{
		std::vector<Segment_2> result;

		Finite_edges_iterator fe = this->finite_edges_begin();
		for(; fe != this->finite_edges_end(); ++fe)
		{
			Segment_2 curr_segment(fe->first->vertex(ccw(fe->second))->point(),
								   fe->first->vertex(cw(fe->second))->point());
			result.push_back(curr_segment);
		}

		return result;
	}

};

#endif // _Gyroviz_border_points_dt2_
