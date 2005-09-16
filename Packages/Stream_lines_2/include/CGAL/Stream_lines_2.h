// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Abdelkrim Mebarki <Abdelkrim.Mebarki@sophia.inria.fr>

#ifndef CGAL_STREAM_LINES_2_H_ 
#define CGAL_STREAM_LINES_2_H_

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <queue>
#include <math.h>

#include <fstream>
#include <iostream>

#include "streamlines_assertions.h"

#define forward 1;
#define backward -1;

CGAL_BEGIN_NAMESPACE

template <class VectorField_2, class Integrator_2>
class Stream_lines_2{
public:
	typedef typename VectorField_2::Vector_field_2                                          Vector_field_2;
	typedef typename VectorField_2::Geom_traits                                             Geom_traits;
	typedef typename VectorField_2::FT                                                      FT;
	typedef typename VectorField_2::Point_2                                                 Point_2;
	typedef typename VectorField_2::Vector_2                                                Vector_2;
protected:
	typedef CGAL::Cartesian<FT>                                                             K1;
	typedef struct Kernel : public K1 {};
	typedef CGAL::Triangulation_vertex_base_2<Kernel>                                       Vb;
	typedef CGAL::Triangulation_face_base_2<Kernel>                                         Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                                     TDS;
	typedef CGAL::Delaunay_triangulation_2<Kernel,TDS>                                      DT;
	typedef typename DT::Vertex_handle                                                      Vertex_handle;
	typedef typename DT::Face_handle                                                        Face_handle;
	typedef typename DT::Face_circulator                                                    Face_circulator;
	typedef typename DT::Edge_iterator                                                      Edge_iterator;
	typedef std::pair<Point_2,FT>                                                           Circle;
	typedef
	CGAL::Quadruple<Vertex_handle,Vertex_handle,Vertex_handle,Circle>                       Pq_element;
	Pq_element                                                                              Biggest_circle;
	FT distance(Point_2 p, Point_2 q){
		return sqrt(((p.x() - q.x())*(p.x() - q.x()))+((p.y() - q.y())*(p.y() - q.y())));};
	int          ir;
	int          il;
	Pq_element   Pq_element_max_r;
	Pq_element   Pq_previous_r,Pq_current_r,Pq_next_r;
	Pq_element   Pq_element_max_l;
	Pq_element   Pq_previous_l,Pq_current_l,Pq_next_l;
public:
	DT           m_DT;
	typedef std::list<Point_2>                                                              Point_container_2;
	typedef typename Point_container_2::iterator                                            Point_iterator_2;
	typedef std::list<std::pair<Point_iterator_2, Point_iterator_2> >                       Iterator_container_2;
	Iterator_container_2                                                                    iterator_container;
	typedef std::list<Vertex_handle>                                                        Vertex_container_2;
	typedef typename  Iterator_container_2::iterator                                        Stream_line_iterator_2;
	typedef std::list<Point_container_2>                                                    Stream_line_container_2;
protected:
	Stream_line_container_2                                        stl_container;
	class C{
	public:
		bool operator()(const Pq_element &a1, const Pq_element &a2){
			return a1.fourth.second < a2.fourth.second ;}};
	std::priority_queue<Pq_element, std::vector<Pq_element>, C>    pq;
	int                                                            iOrder_insertion;
	FT                                                             fSepStl_seed;
	FT                                                             separating_distance;
	FT                                                             saturation_ratio;
	Point_2                                                        seed_point;
	unsigned int                                                   _number_of_lines;
protected:
	void place_stream_lines(const Vector_field_2 & vector_field_2, const Integrator_2 & integrator,
			const int & sampling_step, const bool & step_by_step = false);
	bool get_next_seed_point(FT & distanceg, Point_2 & seed_point);
	FT find_smallest_circle(const Vertex_handle & pVertex_handle);
	Vertex_handle insert_point(const Point_2 & pPoint, FT& fDist,bool bDistanceCalculation);
	Vertex_handle insert_point(const Point_2 & pPoint, const Face_handle & m_Face_handle, FT& fDist,
			bool bDistanceCalculation);
	void integrate_streamline(const Vector_field_2 & vector_field_2,	const Integrator_2 & integrator,
		Point_container_2& stl, Point_2& seed_point, Vertex_container_2& stl_vertices, const int & sampling_step);
	void integrate_forward(const Vector_field_2 & vector_field_2,	const Integrator_2 & integrator,
			Point_container_2& stl,Point_2&	seed_point,
			Vertex_container_2& stl_vertices, const int & sampling_step);
	void integrate_backward(const Vector_field_2 & vector_field_2,	const Integrator_2 & integrator, Point_container_2& 
			stl, Vertex_container_2& stl_vertices, const int & sampling_step);
	void insert_streamline(const Vector_field_2 & vector_field_2, Point_container_2 stl,
			 Vertex_container_2 stl_vertices);
	void pq_elements(const Vector_field_2 & vector_field_2, Vertex_container_2 stl_vertices, int i,
			const Vertex_handle & m_Vertex_handle, int before_end);
	void make_iterator();
public:
	Stream_lines_2(const Vector_field_2 & m_vector_field_2, const Integrator_2 & m_integrator, const FT
			& m_separating_distance, const FT & m_saturation_ratio, const int & sampling_insertion = 0, const bool & 
			step_by_step = false);
	bool continue_next(const Vector_field_2 & vector_field_2, const Integrator_2 & integrator, const int & sampling_step);
	Stream_line_iterator_2 begin();
	Stream_line_iterator_2 end();
	// for visualizing streamlines
	void print_stream_lines(std::ofstream & fw);
	std::list<Point_2> get_pq();
// 	void print_stream_lines_eps(std::ofstream & fw);
	unsigned int number_of_lines(){return _number_of_lines;};
	std::list< std::pair<Point_2, Point_2> > get_tr()
	{
		std::list< std::pair<Point_2, Point_2> > _list;
		Edge_iterator eit = m_DT.edges_begin();
		for (;eit != m_DT.edges_end();eit++)
		{
			Point_2 p1 = (*eit).first->vertex(m_DT.ccw((*eit).second))->point();
			Point_2 p2 = (*eit).first->vertex(m_DT.cw((*eit).second))->point();
			_list.push_front(std::pair<Point_2, Point_2>(p1, p2));
		}
		return _list;
	}	
	std::pair<Point_2, FT> get_biggest_circle()
	{
		Pq_element m_Pq = Biggest_circle;
		std::pair<Point_2, FT> circle(m_Pq.fourth.first, m_Pq.fourth.second);
		return circle;
	}
protected:
	FT max_x;
	FT min_x;
	FT max_y;
	FT min_y;
protected:
	Stream_line_iterator_2 begin_iterator;
	Stream_line_iterator_2 end_iterator;	
private:
	int number_of_points;
};

template <class VectorField_2, class Integrator_2>
Stream_lines_2<VectorField_2, Integrator_2>::Stream_lines_2(const Vector_field_2 &
vector_field_2, const Integrator_2 & m_integrator, const FT & m_separating_distance, const FT
& m_saturation_ratio, const int & sampling_step, const bool & step_by_step){
	separating_distance = m_separating_distance;
	saturation_ratio = m_saturation_ratio;
	ir = il = 0; // initialization
	fSepStl_seed = separating_distance*saturation_ratio;
	max_x = vector_field_2.iso_rectangle().xmax();
	min_x = vector_field_2.iso_rectangle().xmin();
	max_y = vector_field_2.iso_rectangle().ymax();
	min_y = vector_field_2.iso_rectangle().ymin();
	m_DT.clear();
	Point_2 pPoint;
	pPoint = Point_2(min_x-separating_distance,min_y-separating_distance);
	m_DT.insert(pPoint);
	pPoint = Point_2(min_x-separating_distance,max_y+separating_distance);
	m_DT.insert(pPoint);
	pPoint = Point_2(max_x+separating_distance,min_y-separating_distance);
	m_DT.insert(pPoint);
	pPoint = Point_2(max_x+separating_distance,max_y+separating_distance);
	m_DT.insert(pPoint);
	for (int i=(int) (min_x-separating_distance);i<max_x+(int) separating_distance;i=i+(int) (fSepStl_seed)){
		pPoint = Point_2((FT)i,(FT)max_y+separating_distance);
		m_DT.insert(pPoint);
		pPoint = Point_2((FT)max_x+separating_distance,(FT)i);
		m_DT.insert(pPoint);
		pPoint = Point_2((FT)i,min_y-separating_distance);
		m_DT.insert(pPoint);
		pPoint = Point_2(min_x-separating_distance,(FT)i);
		m_DT.insert(pPoint);}
		_number_of_lines = 0;
	place_stream_lines(vector_field_2, m_integrator, sampling_step, step_by_step);}

template <class VectorField_2, class Integrator_2>
void Stream_lines_2<VectorField_2, Integrator_2>::place_stream_lines(const Vector_field_2 & vector_field_2,	const Integrator_2 & integrator, const int & sampling_step, const bool & step_by_step)
{
	seed_point = Point_2((max_x+min_x)/2.0,(max_y+min_y)/2.0);
// the first chosen point can be not valid
	FT xrange = max_x - min_x;
	FT yrange = max_y - min_y;
	while(!vector_field_2.is_in_domain(seed_point))
	{
		std::cout << "searching valid seed point..\n";
		FT x = min_x + (FT) (((FT) rand() * xrange)/((FT) RAND_MAX));
		FT y = min_y + (FT) (((FT) rand() * yrange)/((FT) RAND_MAX));
		seed_point = Point_2(x, y);
	}
	std::cout << seed_point << " first seed point\n";
	std::cout << "creating the placement..\n";
	FT distance = (FT) max_x * (1.0/2.0);
	bool b = (distance>fSepStl_seed);
// 	int i=0;
	if (!step_by_step)
		while(b)
		{
			Point_container_2 stl;
			Vertex_container_2 stl_vertices;
			integrate_streamline(vector_field_2, integrator, stl, seed_point, stl_vertices, sampling_step);
			insert_streamline(vector_field_2, stl, stl_vertices);
			_number_of_lines++;
			b = get_next_seed_point(distance,seed_point);
		}
	else
	{
		Point_container_2 stl;
		Vertex_container_2 stl_vertices;
		integrate_streamline(vector_field_2, integrator, stl, seed_point, stl_vertices, sampling_step);
		insert_streamline(vector_field_2, stl, stl_vertices);
		_number_of_lines++;
		b = get_next_seed_point(distance,seed_point);
	}
	make_iterator();
}
	
template <class VectorField_2, class Integrator_2>
bool Stream_lines_2<VectorField_2, Integrator_2>::continue_next(const Vector_field_2 & vector_field_2, const Integrator_2 & integrator, const int & sampling_step)
{
	FT distance;
	Point_container_2 stl;
	Vertex_container_2 stl_vertices;
	integrate_streamline(vector_field_2, integrator, stl, seed_point, stl_vertices, sampling_step);
	insert_streamline(vector_field_2, stl, stl_vertices);
	_number_of_lines++;
	make_iterator();
	return get_next_seed_point(distance,seed_point);
}

// get the next seed point
template <class VectorField_2, class Integrator_2>
void 
Stream_lines_2<VectorField_2,Integrator_2>::integrate_streamline(const Vector_field_2 & vector_field_2, const Integrator_2 & integrator, Point_container_2& stl, Point_2&
seed_point, Vertex_container_2& stl_vertices, const int & sampling_step)
{
	integrate_forward(vector_field_2, integrator, stl, seed_point, stl_vertices, sampling_step);
	integrate_backward(vector_field_2, integrator, stl, stl_vertices, sampling_step);
}

template <class VectorField_2, class Integrator_2>
inline
typename Stream_lines_2<VectorField_2, Integrator_2>::FT 
Stream_lines_2<VectorField_2, Integrator_2>::find_smallest_circle(const Vertex_handle & pVertex_handle)
{
	Face_circulator pFace_handle = m_DT.incident_faces(pVertex_handle);
	Face_circulator pEnd = pFace_handle;
	FT fMin = max_x;
	CGAL_For_all(pFace_handle,pEnd)
	{
		FT fDist =
			CGAL::squared_radius(
				pFace_handle->vertex(0)->point(),
				pFace_handle->vertex(1)->point(),
				pFace_handle->vertex(2)->point()) * 4.0;
		fDist = sqrt(fDist);
		if (fDist < fMin)
		{
			fMin = fDist;
		}
	}
	return fMin;
}

template <class VectorField_2, class Integrator_2>
inline 
typename Stream_lines_2<VectorField_2, Integrator_2>::Vertex_handle 
Stream_lines_2<VectorField_2, Integrator_2>::insert_point(const Point_2 & pPoint, FT& fDist,bool bDistanceCalculation)
{
	Vertex_handle pVertex_handle = m_DT.insert(pPoint);
	if (bDistanceCalculation)
		fDist = find_smallest_circle(pVertex_handle);
	else
		fDist = 0.0;
	return (pVertex_handle);
}

template <class VectorField_2, class Integrator_2>
inline
typename Stream_lines_2<VectorField_2, Integrator_2>::Vertex_handle
Stream_lines_2<VectorField_2,Integrator_2>::insert_point(const Point_2 & pPoint, const Face_handle & m_Face_handle, FT& fDist,bool bDistanceCalculation)
{
	Vertex_handle pVertex_handle = m_DT.insert(pPoint,m_Face_handle);
	if (bDistanceCalculation)
		fDist = find_smallest_circle(pVertex_handle);
	else
		fDist = 0.0;
	return (pVertex_handle);
}

template <class VectorField_2, class Integrator_2>
void
Stream_lines_2<VectorField_2, Integrator_2>::integrate_forward(const Vector_field_2 & vector_field_2, const Integrator_2 & integrator, Point_container_2& stl, Point_2& seed_point, Vertex_container_2& stl_vertices, const int & sampling_step)
{
	int sampling = 0; // sampling step;
	int insertion = 0; // insertion order;
	int insertion_step = 0;
	Point_container_2 list_of_point;
	Vertex_container_2 list_of_vertex;
	number_of_points = 0;
	Point_2 pPoint1;
	bool bEnd = false;
	FT dist;
	Point_2 new_point = Point_2 (seed_point.x(), seed_point.y());
	Vertex_handle m_Vertex_handle = insert_point(new_point, dist, true);
	stl_vertices.push_front(m_Vertex_handle);
	stl.push_front(new_point);
	number_of_points++;
	Point_2 old_point;
	insertion_step = (int) (((dist)-fSepStl_seed) / std::max((FT) sampling_step,vector_field_2.get_integration_step()));
	if (insertion_step < 0) insertion_step = 0;
	while (!bEnd)
	{
		Point_2 ex_old_point = old_point;
		old_point = new_point;
		CGAL_streamlines_precondition(vector_field_2.is_in_domain(old_point));
		new_point = integrator(old_point,vector_field_2,true);
		bEnd = !vector_field_2.is_in_domain(new_point);
		bEnd = bEnd || (new_point == old_point);/* to review */
		if(number_of_points > 30)
			bEnd = bEnd || ((distance(stl.front(), stl.back()))<vector_field_2.get_integration_step());
		FT dist_ = distance(ex_old_point,new_point);
		bEnd = bEnd || dist_ <= vector_field_2.get_integration_step();
		if (!bEnd)
		{
			if(sampling != sampling_step)
			{
				stl.push_front(new_point);
				number_of_points ++;
				sampling++;
			}
			else
			{
				if (insertion != insertion_step)
				{
					stl.push_front(new_point); 
					number_of_points++;
					insertion++;
					list_of_point.push_front(new_point);
				}
				else
				{
					stl.push_front(new_point);
					number_of_points++;
					list_of_point.push_front(new_point);
					list_of_point.pop_front();
					m_Vertex_handle = insert_point(new_point, stl_vertices.front()->face(), dist, true);
					while ((dist <= separating_distance)&&(!list_of_point.empty()))
					{
						m_DT.remove(m_Vertex_handle);
						for (int i=0;i<=sampling_step;i++){
							stl.pop_front();
							number_of_points--;}
						new_point = list_of_point.front();
						list_of_point.pop_front();
						m_Vertex_handle = insert_point(new_point, stl_vertices.front()->face(), dist,true);
					}
						// adaptive insertion order coefficient
					insertion_step = (int) (((dist)-fSepStl_seed) /
							std::max((FT) sampling_step,vector_field_2.get_integration_step()));
					if (insertion_step < 0) insertion_step = 0;
 					list_of_vertex.push_front(m_Vertex_handle);
					(bEnd) = (((bEnd))||(dist<separating_distance));
					while (!list_of_point.empty())
					{
						Point_2 p = list_of_point.front();
						m_Vertex_handle = insert_point(p, stl_vertices.front()->face(), dist, false);
						list_of_vertex.push_front(m_Vertex_handle);
						list_of_point.pop_front();
					}
					while(!list_of_vertex.empty())
					{
						stl_vertices.push_front(list_of_vertex.front());
						list_of_vertex.pop_front();
					}
					insertion = 0;
				}
				sampling = 0;
			}
		}
		else
		{
			if (!list_of_point.empty())
			{
				new_point = list_of_point.front();
				list_of_point.pop_front();
				Vertex_handle m_Vertex_handle = insert_point(new_point, dist, true);
				while ((dist <= separating_distance)&&(!list_of_point.empty()))
				{
					m_DT.remove(m_Vertex_handle);
					for (int i=0;i<=sampling_step;i++)
					{
						stl.pop_front();
						number_of_points--;
					}
					new_point = list_of_point.front();
					list_of_point.pop_front();
					m_Vertex_handle = insert_point(new_point, stl_vertices.front()->face(), dist, true);
				}
				insertion_step =	(int) (((dist)-fSepStl_seed) / 
						std::max((FT) sampling_step,vector_field_2.get_integration_step()));
				if (insertion_step < 0) insertion_step = 0;
				(bEnd) = (((bEnd))||(dist<separating_distance));
			}
			while (!list_of_point.empty())
			{
				Point_2 p = list_of_point.front();
				m_Vertex_handle = insert_point(p, stl_vertices.front()->face(), dist, false);
 				list_of_vertex.push_front(m_Vertex_handle);
				list_of_point.pop_front();
			}
			while(!list_of_vertex.empty())
			{
				stl_vertices.push_front(list_of_vertex.front());
				list_of_vertex.pop_front();
			}
		}
	}
}

template <class VectorField_2, class Integrator_2>
void Stream_lines_2<VectorField_2, Integrator_2>::integrate_backward(const Vector_field_2 & vector_field_2, const Integrator_2 & integrator, Point_container_2& stl, Vertex_container_2& stl_vertices, const int & sampling_step)
{
	int sampling = 0; // sampling step;
	int insertion = 0; // insertion order;
	int insertion_step = 0;
	Point_container_2 list_of_point;
	Vertex_container_2 list_of_vertex;
	Point_2 pPoint1;
	bool bEnd = false;
	FT dist;
	Point_2 new_point = Point_2 (stl.back().x(),stl.back().y());
	Vertex_handle m_Vertex_handle = insert_point(new_point, stl_vertices.back()->face(), dist,true);
	stl_vertices.push_back(m_Vertex_handle);
	stl.push_back(new_point);
	number_of_points++;
	Point_2 old_point;
	while (!bEnd)
	{
		Point_2 ex_old_point = old_point;
		old_point = new_point;
		std::pair<Vector_2, FT> field_vector;
		CGAL_streamlines_precondition(vector_field_2.is_in_domain(old_point));
		new_point = integrator(old_point,vector_field_2,false);
		bEnd = !vector_field_2.is_in_domain(new_point);
		FT dist_ = distance(ex_old_point,new_point);
		bEnd = bEnd || dist_ <= vector_field_2.get_integration_step() || (new_point == old_point);/* to review */	
		if(number_of_points > 30)
			bEnd = bEnd || ((distance(stl.front(), stl.back()))<vector_field_2.get_integration_step());
// 		bEnd = bEnd || (number_of_points > 3000);
		if (!bEnd)
		{
			if(sampling != sampling_step)
			{
				stl.push_back(new_point);
				number_of_points ++;
				sampling++;
			}
			else
			{
				if (insertion != insertion_step)
				{
					stl.push_back(new_point); 
					number_of_points++;
					insertion++;
					list_of_point.push_back(new_point);
				}
				else
				{
					stl.push_back(new_point);
					number_of_points++;
					list_of_point.push_back(new_point);
					list_of_point.pop_back();
					m_Vertex_handle = insert_point(new_point, stl_vertices.back()->face(), dist,true);
					while ((dist <= separating_distance)&&(!list_of_point.empty()))
					{
						m_DT.remove(m_Vertex_handle);
						for (int i=0;i<=sampling_step;i++)
						{
							stl.pop_back();
							number_of_points--;
						}
						new_point = list_of_point.back();
						list_of_point.pop_back();
						m_Vertex_handle = insert_point(new_point, stl_vertices.back()->face(), dist,true);
					}
					// adaptive insertion order coefficient
					insertion_step =	(int) (((dist)-fSepStl_seed) / 
							std::max((FT) sampling_step,vector_field_2.get_integration_step()));
					if (insertion_step < 0) insertion_step = 0;
 					list_of_vertex.push_back(m_Vertex_handle);
					(bEnd) = (((bEnd))||(dist<separating_distance));
					while (!list_of_point.empty())
					{
						Point_2 p = list_of_point.back();
						m_Vertex_handle = insert_point(p, stl_vertices.back()->face(), dist, false);
						list_of_vertex.push_front(m_Vertex_handle);
						list_of_point.pop_back();
					}
					while(!list_of_vertex.empty())
					{
						stl_vertices.push_back(list_of_vertex.back());
						list_of_vertex.pop_back();
					}
					insertion = 0;
				}
			sampling = 0;
			}
		}
		else
		{
			if (!list_of_point.empty())
			{
				new_point = list_of_point.back();
				list_of_point.pop_back();
				Vertex_handle m_Vertex_handle = insert_point(new_point, stl_vertices.back()->face(), dist, true);
				while ((dist <= separating_distance)&&(!list_of_point.empty()))
				{
					m_DT.remove(m_Vertex_handle);
					for (int i=0;i<=sampling_step;i++)
					{
						stl.pop_back();
						number_of_points--;
					}
					new_point = list_of_point.back();
					list_of_point.pop_back();
					m_Vertex_handle = insert_point(new_point, stl_vertices.back()->face(), dist, true);
				}
				// adaptive insertion order coefficient
				insertion_step =	(int) (((dist)-fSepStl_seed) / 
						std::max((FT) sampling_step,vector_field_2.get_integration_step()));
				if (insertion_step < 0) insertion_step = 0;
// 				list_of_vertex.push_front(m_Vertex_handle);
				(bEnd) = (((bEnd))||(dist<separating_distance));
			}
			while (!list_of_point.empty())
			{
				Point_2 p = list_of_point.back();
				m_Vertex_handle = insert_point(p, stl_vertices.back()->face(), dist, false);
 				list_of_vertex.push_back(m_Vertex_handle);
				list_of_point.pop_back();
			}
			while(!list_of_vertex.empty())
			{
				stl_vertices.push_back(list_of_vertex.back());
				list_of_vertex.pop_back();
			}
		}
	}
}

template <class VectorField_2, class Integrator_2>
inline void
Stream_lines_2<VectorField_2, Integrator_2>::
insert_streamline(const Vector_field_2 & vector_field_2, Point_container_2 stl, Vertex_container_2 stl_vertices){
	stl_container.push_back(stl);
	Vertex_handle m_Vertex_handle = NULL;
	int i = 1;
	unsigned int size_ = (int) (stl_vertices.size());
	ir = il = 0;
	while (!stl_vertices.empty()){
		pq_elements(vector_field_2, stl_vertices, i, m_Vertex_handle, size_);
		m_Vertex_handle = stl_vertices.front();
		stl_vertices.pop_front();
		i++;}}

template <class VectorField_2, class Integrator_2>
void Stream_lines_2<VectorField_2, Integrator_2>::
pq_elements(const Vector_field_2 & vector_field_2, Vertex_container_2 stl_vertices, int i,
		const Vertex_handle & m_Vertex_handle, int size_){
	if ((i!=0) && (i!=(size_)) && (std::div(i,10).rem!=0)){
			Vertex_handle pVertex_handle = stl_vertices.front();
			Face_handle pFace_handle;
			int iIndex;
			if (m_DT.is_edge(pVertex_handle,m_Vertex_handle,pFace_handle,iIndex)){
				Point_2 p0 = pVertex_handle->point();
				Point_2 c = m_DT.circumcenter(pFace_handle);
				FT fDist = distance(p0,c);
				bool b = vector_field_2.is_in_domain(c) && (fDist >= fSepStl_seed);
				if (b){
					Circle pCircle(c,fDist);
					Pq_element m_Pq_element = Pq_element(
						pFace_handle->vertex(0),
						pFace_handle->vertex(1),
						pFace_handle->vertex(2),
						pCircle);
					if (ir == 0){
						Pq_previous_r = m_Pq_element;
						Pq_element_max_r = m_Pq_element;
						ir++;}
					else if (ir == 1){
						Pq_current_r = m_Pq_element;ir++;}
					else{
						Pq_next_r = m_Pq_element;
						if (Pq_element_max_r.fourth.second <= Pq_next_r.fourth.second)
							Pq_element_max_r = Pq_next_r;
						if ((Pq_current_r.fourth.second>=Pq_previous_r.fourth.second)
								&&(Pq_current_r.fourth.second>=Pq_next_r.fourth.second)){
							pq.push(Pq_current_r);}
						Pq_previous_r = Pq_current_r;
						Pq_current_r = Pq_next_r;
						ir++;}}
				p0 = pFace_handle->neighbor(iIndex)->vertex(0)->point();
				c = m_DT.circumcenter(pFace_handle->neighbor(iIndex));
				fDist = distance(p0,c);
				b = vector_field_2.is_in_domain(c) && (fDist >= fSepStl_seed);
				if (b){
					Circle pCircle(c,fDist);
					Pq_element m_Pq_element = Pq_element(
							pFace_handle->neighbor(iIndex)->vertex(0),
							pFace_handle->neighbor(iIndex)->vertex(1),
							pFace_handle->neighbor(iIndex)->vertex(2),
							pCircle);
				if (il == 0){
					Pq_previous_l = m_Pq_element;
					Pq_element_max_l = m_Pq_element;
					il++;}
				else if (il == 1){
					Pq_current_l = m_Pq_element;il++;}
				else{
					Pq_next_l = m_Pq_element;
					if (Pq_element_max_l.fourth.second <= Pq_next_l.fourth.second) 
						Pq_element_max_l = Pq_next_l;
					if ((Pq_current_l.fourth.second>=Pq_previous_l.fourth.second)
							&&(Pq_current_l.fourth.second>=Pq_next_l.fourth.second)){
						pq.push(Pq_current_l);}
					Pq_previous_l = Pq_current_l;
					Pq_current_l = Pq_next_l;
					il++;}}}
		if ((ir+il == (int) size_-2)&&(size_>2)){
			pq.push(Pq_element_max_l);
			pq.push(Pq_element_max_r);}}
	else{
		Vertex_handle pVertex_handle = stl_vertices.front();
		Face_circulator pFace_handle = m_DT.incident_faces(pVertex_handle);
		Face_circulator pEnd = pFace_handle;
		CGAL_For_all(pFace_handle,pEnd){
			Point_2 p0 = pFace_handle->vertex(0)->point();
			Point_2 c = m_DT.circumcenter(pFace_handle);
			bool b = vector_field_2.is_in_domain(c);
			if (b){
				FT fDist = distance(p0,c);
				if (fDist >= fSepStl_seed){
					Circle pCircle(c,fDist);
					Pq_element m_Pq_element = Pq_element(
							pFace_handle->vertex(0),
							pFace_handle->vertex(1),
							pFace_handle->vertex(2),
							pCircle);
					pq.push(m_Pq_element);}}}}}

// get the next seed point
template <class VectorField_2, class Integrator_2>
inline 
bool
Stream_lines_2<VectorField_2, Integrator_2>::get_next_seed_point(FT & distance, Point_2 & seed_point){
	Vertex_handle v0, v1, v2; 
	Face_handle fr;
	bool b0,b;
	Pq_element m_Pq_element;
	do{
		m_Pq_element = pq.top();
		v0 = m_Pq_element.first;
		v1 = m_Pq_element.second;
		v2 = m_Pq_element.third;
		distance = m_Pq_element.fourth.second;
		pq.pop();
		b0 = m_DT.is_face(v0,v1,v2,fr);
		if (b0){
			seed_point = m_Pq_element.fourth.first;}
		b = (!pq.empty());
	}while ((b)&&(!b0));
	Biggest_circle = m_Pq_element;
	return b;}
	
template <class VectorField_2, class Integrator_2>
typename Stream_lines_2<VectorField_2, Integrator_2>::Stream_line_iterator_2
Stream_lines_2<VectorField_2, Integrator_2>::begin(){
	return begin_iterator;}

template <class VectorField_2, class Integrator_2>
typename Stream_lines_2<VectorField_2, Integrator_2>::Stream_line_iterator_2
Stream_lines_2<VectorField_2, Integrator_2>::end(){
	return end_iterator;}

template <class VectorField_2, class Integrator_2>
inline
void Stream_lines_2<VectorField_2, Integrator_2>::make_iterator(){ 
	for(typename Stream_line_container_2::iterator begin=stl_container.begin(); begin!=stl_container.end();begin++){
		std::pair<Point_iterator_2, Point_iterator_2>
		iterator_pair((*begin).begin(), (*begin).end());
		iterator_container.push_front(iterator_pair);}
	begin_iterator = iterator_container.begin();
	end_iterator = iterator_container.end();}

// output an stl file
template <class VectorField_2, class Integrator_2>
void Stream_lines_2<VectorField_2, Integrator_2>::print_stream_lines(std::ofstream & fw)
{
	typename Stream_line_container_2::iterator begin_iterator;
	Stream_line_container_2 stl_container_temp = stl_container;
	Point_container_2 stl;
	fw << max_x - min_x << " " << max_y - min_y << "\n";
	fw << stl_container.size() << "\n";
	for(begin_iterator=stl_container_temp.begin();begin_iterator!=stl_container_temp.end();begin_iterator++){
		fw << (*begin_iterator).size() << "\n";
		typename Point_container_2::iterator begin_point_iterator = (*begin_iterator).begin();
		typename Point_container_2::iterator end_point_iterator = (*begin_iterator).end();
		FT i , j;
		for(;begin_point_iterator!=end_point_iterator;begin_point_iterator++){
			Point_2 p = *begin_point_iterator;
			i = p.x() - min_x;
			j = p.y() - min_y;
				fw << i << " " << j << "\n";}};
	fw.close();}

template <class VectorField_2, class Integrator_2>
std::list<typename Stream_lines_2<VectorField_2, Integrator_2>::Point_2> 
Stream_lines_2<VectorField_2, Integrator_2>::get_pq()
{
	std::list<Point_2> _list;
	std::priority_queue<Pq_element, std::vector<Pq_element>, C> pq_temp;
	pq_temp = pq;
	while (!pq_temp.empty())
	{
		Pq_element m_Pq_element = pq_temp.top();
		Vertex_handle v0 = m_Pq_element.first;
		Vertex_handle v1 = m_Pq_element.second;
		Vertex_handle v2 = m_Pq_element.third;
		pq_temp.pop();
		Face_handle fr;
		bool b0 = m_DT.is_face(v0,v1,v2,fr);
		Point_2 sdPoint = m_Pq_element.fourth.first;
		if (b0)
			_list.push_front(sdPoint);
	}
	return _list;
}

CGAL_END_NAMESPACE

#endif
