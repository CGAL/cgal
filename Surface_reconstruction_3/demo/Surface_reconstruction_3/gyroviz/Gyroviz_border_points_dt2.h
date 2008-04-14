// Author     : Nader Salman

#ifndef _Gyroviz_border_points_dt2_
#define _Gyroviz_border_points_dt2_

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Gyroviz_vertex_segment_2.h"
#include <CGAL/Delaunay_triangulation_2.h>

template < class Gt, class Tds, class Cdt2VertexHandle >
class Gyroviz_border_points_dt2 :  public CGAL::Delaunay_triangulation_2<Gt, Tds>
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
	typedef typename Geom_traits::Point_2			   Point_2;
	typedef typename Geom_traits::Segment_2			   Segment_2;
	typedef typename Base::Face_handle				   Face_handle;
	typedef typename Base::Vertex_handle			   Vertex_handle;
	typedef typename Base::Vertex					   Vertex;
	typedef typename Base::Edge						   Edge;
	typedef typename Base::Finite_edges_iterator	   Finite_edges_iterator;
	typedef typename Base::Finite_faces_iterator       Finite_faces_iterator;
	typedef typename Base::Finite_vertices_iterator    Finite_vertices_iterator;
	typedef typename Gyroviz_vertex_segment_2<Cdt2VertexHandle> Gyroviz_vertex_segment_2;
	
//	// Data members
//private:
//	std::vector<Vertex_handle> input;

	// Public methods
public:


	// Constructors
	Gyroviz_border_points_dt2(){}
	Gyroviz_border_points_dt2(std::vector<Cdt2VertexHandle> original_vertices)
	{
		for(int i=0; i<original_vertices.size(); ++i)
		{
			Vertex_handle vh = this->insert(original_vertices[i]->point());
			//vh->info().set_ptr(&*original_vertices[i]);
			vh->info() = original_vertices[i];
		}
	}


	// Functions

	// this function will store as a vector of Gyroviz_vertex_segment_2 
	// the entire 2D Delaunay triangulation
	std::vector<Gyroviz_vertex_segment_2> segments_out_of_dt2()
	{
		std::vector<Gyroviz_vertex_segment_2> result;

		Finite_edges_iterator fe = this->finite_edges_begin();
		for(; fe != this->finite_edges_end(); ++fe)
		{
			Vertex_handle first_vertex(fe->first->vertex(ccw(fe->second)));
			//result.push_back(first_vertex);
			Vertex_handle second_vertex(fe->first->vertex(cw(fe->second)));
			//result.push_back(second_vertex);
			//Vertex_handle original_first_vertex  = (Vertex*) first_vertex->info().get_ptr();
			//Vertex_handle original_second_vertex = (Vertex*) second_vertex->info().get_ptr();
			Cdt2VertexHandle original_first_vertex  = first_vertex->info();
			Cdt2VertexHandle original_second_vertex = second_vertex->info();
			Gyroviz_vertex_segment_2 curr_segment(original_first_vertex, original_second_vertex);
			result.push_back(curr_segment);
		}

		return result;
	}

};

#endif // _Gyroviz_border_points_dt2_
