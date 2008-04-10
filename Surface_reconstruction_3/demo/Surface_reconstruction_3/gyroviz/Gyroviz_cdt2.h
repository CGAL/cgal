// Author     : Nader Salman

#ifndef _Gyroviz_cdt2_
#define _Gyroviz_cdt2_


// the idea in here is to read a pnt file
// segment the corresponding image
// extract the vertices on the border of the segmentation (3x3 mask)
// link all these vertices using segments
// Adapt Smith-Waterman/Needleman-Wunsch algorithms 
// for performing local/global sequence alignement to filter these segments
// the remaining segments will be used as constraints in the 2D Delaunay triangulation

#include <algorithm>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/origin.h>

#include "Gyroviz_info_for_cdt2.h"
//#include "Gyroviz_vertex_segment_2.h"
#include "Gyroviz_border_points_dt2.h"
#include "Gyroviz_triangle_with_cam.h"
#include "Gyroviz_constrained_triangle_soup.h"

#include <list>
#include <vector>
#include <sstream>
#include <cmath>

#include <CImg.h>
using namespace cimg_library;

// Tds vertex class must inherit from Triangulation_vertex_base_with_info_2<CGAL::Point_3,K>
template < class Gt, class Tds, class Itag >
class Gyroviz_cdt2 : public CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>
{
	// Private types
private:

	typedef CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>  Base;

	// Public types
public:

	// Repeat Constrained_Delaunay_triangulation_2 public types
	typedef Tds Triangulation_data_structure;
	typedef CGAL::Exact_predicates_tag                 Itag;
	typedef Gt  Geom_traits;
	typedef typename Geom_traits::FT FT;
	typedef typename Geom_traits::Point_2              Point_2;
	typedef typename Geom_traits::Point_3              Point_3;
	typedef typename Geom_traits::Vector_2             Vector_2;
	typedef typename Geom_traits::Vector_3             Vector_3;
	typedef typename Geom_traits::Segment_2            Segment_2;
	typedef typename Geom_traits::Segment_3            Segment_3;
	typedef typename Geom_traits::Triangle_2		   Triangle_2;
	typedef typename Geom_traits::Triangle_3		   Triangle_3;
	typedef typename Geom_traits::Sphere_3             Sphere;
	typedef typename Geom_traits::Iso_cuboid_3         Iso_cuboid_3;
	typedef typename Base::Face_handle                 Face_handle;
	typedef typename Base::Vertex_handle               Vertex_handle;
	typedef typename Base::Edge                        Edge;
	typedef typename Base::Edge_circulator             Edge_circulator;
	typedef typename Base::Finite_edges_iterator       Finite_edges_iterator;
	typedef typename Base::Finite_faces_iterator       Finite_faces_iterator;
	typedef typename Base::Finite_vertices_iterator    Finite_vertices_iterator;
	//typedef typename Gyroviz_vertex_segment_2<Gyroviz_cdt2>    Gyroviz_vertex_segment_2;
	typedef typename Gyroviz_triangle_with_cam<Gyroviz_cdt2>   Gyroviz_triangle_with_cam;
	typedef typename Gyroviz_constrained_triangle_soup<Gyroviz_cdt2>   Gyroviz_constrained_triangle_soup;

	// Data members
private:

	Iso_cuboid_3 m_bounding_box; // Triangulation's bounding box
	Point_3 m_barycenter; // Triangulation's barycenter
	FT m_standard_deviation; // Triangulation's standard deviation

	// Public methods
public:

	// Default constructor, copy constructor and operator =() are fine

	bool save_pnt(char *pFilename)
	{
		// TODO
	}

	bool read_pnt(const char *pFilename)
	{

		//DEBUG
		int number_of_vertices_pnt = 0;

		//extract from pFilename the image number
		std::string temp ( pFilename );
		std::string filename_without_path = temp.substr(temp.size()-15);
		std::string extract_number = filename_without_path.substr(7,filename_without_path.size()-4);

		int image_number = atoi(extract_number.c_str());

		FILE *pFile = fopen(pFilename,"r");
		if(pFile == NULL)
			return false;

		//scan vertices and add them to triangulation with corresponding info
		int lineNumber = 0;
		char pLine[512];

		while( fgets(pLine, 512, pFile))
		{
			lineNumber++;

			// read 2D/3D coordinates 
			//(on suppose avoir fait du traitement 
			// des fichiers pnt a l'etape precedente)
			int  unused1, is_reconstructed;
			double p,q,x,y,z;

			if (sscanf(pLine,"%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf", &p,&q,&unused1,&is_reconstructed,&x,&y,&z) == 7)
			{
				if (is_reconstructed == 2)
				{
					Point_2 point_2(p,q);
					Point_3 point_3(x,y,z);

					Vertex_handle vh = this->insert(point_2);
					vh->info() = Gyroviz_info_for_cdt2(point_3,image_number,false);
					number_of_vertices_pnt++;
				}
			}
		}
		fclose(pFile);
		update_bounding_box();
		return (this->number_of_vertices() > 0);
	}


	// set flag to true for the vertices positionned on the border and return a vector of these
	std::vector<Vertex_handle> set_on_border_2D_vertices(const CImg <unsigned char> & image)
	{

		//DEBUG
		int number_of_border_vertices = 0;

		std::vector<Vertex_handle> vector_of_border_points;

		Finite_vertices_iterator fv = this->finite_vertices_begin();
		for(; fv != this->finite_vertices_end(); ++fv)
		{
			// if pixel any of 3x3 surrounding pixels is on border 
			// keep vertex(white color is used in frei-chen gradient operator)

			if (image((unsigned int)fv->point().x(),(unsigned int)fv->point().y(),0,0)==255)
			{
				fv->info().set_flag(true);
				vector_of_border_points.push_back(fv);
			}

			else if((unsigned int)fv->point().x() == 0 && (unsigned int)fv->point().y() == 0) //upper left pixel
			{
				if (image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x(),(unsigned int)fv->point().y()+1,0,0)  == 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y()+1,0,0)== 255)
				{  
					fv->info().set_flag(true);
					vector_of_border_points.push_back(fv);
				}

			}
			else if((unsigned int)fv->point().x() == image.dimx() && (unsigned int)fv->point().y() == 0) //upper right pixel
			{
				if (image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x(),(unsigned int)fv->point().y()+1,0,0)  == 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y()+1,0,0)== 255) 
				{
					fv->info().set_flag(true);
					vector_of_border_points.push_back(fv);
				}
			}

			else if((unsigned int)fv->point().x() == 0 && (unsigned int)fv->point().y() == image.dimy()) //lower left pixel
			{
				if (image((unsigned int)fv->point().x(),(unsigned int)fv->point().y()-1,0,0)  == 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y()-1,0,0)== 255) 
				{
					fv->info().set_flag(true);
					vector_of_border_points.push_back(fv);
				}

			}

			else if((unsigned int)fv->point().x() == image.dimx() && (unsigned int)fv->point().y() == image.dimy()) //lower right pixel
			{
				if (image((unsigned int)fv->point().x(),(unsigned int)fv->point().y()-1,0,0)  == 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y()-1,0,0)== 255) 
				{
					fv->info().set_flag(true);
					vector_of_border_points.push_back(fv);
				}
			}

			else if((unsigned int)fv->point().x() > 0 && (unsigned int)fv->point().x() < image.dimx() && (unsigned int)fv->point().y() == 0) // upper band
			{
				if (image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y()+1,0,0)== 255 ||
					image((unsigned int)fv->point().x(),  (unsigned int)fv->point().y()+1,0,0)== 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y()+1,0,0)== 255) 
				{
					fv->info().set_flag(true);
					vector_of_border_points.push_back(fv);
				}      
			}

			else if((unsigned int)fv->point().x() > 0 && (unsigned int)fv->point().x() < image.dimx() && (unsigned int)fv->point().y() == image.dimy()) // lower band
			{
				if (image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x(),  (unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y(),0,0)  == 255) 
				{
					fv->info().set_flag(true);
					vector_of_border_points.push_back(fv);
				}
			}

			else if((unsigned int)fv->point().x() == 0 && (unsigned int)fv->point().y() > 0  && (unsigned int)fv->point().y() < image.dimy()) // left band
			{
				if (image((unsigned int)fv->point().x(),  (unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y()+1,0,0)== 255 ||
					image((unsigned int)fv->point().x(),  (unsigned int)fv->point().y()+1,0,0)== 255) 
				{
					fv->info().set_flag(true);
					vector_of_border_points.push_back(fv);
				}
			}

			else if((unsigned int)fv->point().x() == image.dimx() && (unsigned int)fv->point().y() > 0  && (unsigned int)fv->point().y() < image.dimy()) // right band
			{
				if (image((unsigned int)fv->point().x(),  (unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y()+1,0,0)== 255 ||
					image((unsigned int)fv->point().x(),  (unsigned int)fv->point().y()+1,0,0)== 255) 
				{
					fv->info().set_flag(true);
					vector_of_border_points.push_back(fv);
				}
			}

			else // middle of the image corner and bands excluded
			{
				if (image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x()-1,(unsigned int)fv->point().y()+1,0,0)== 255 ||
					image((unsigned int)fv->point().x(),  (unsigned int)fv->point().y()+1,0,0)== 255 ||
					image((unsigned int)fv->point().x(),  (unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y()-1,0,0)== 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y(),0,0)  == 255 ||
					image((unsigned int)fv->point().x()+1,(unsigned int)fv->point().y()+1,0,0)== 255)
				{
					fv->info().set_flag(true);
					vector_of_border_points.push_back(fv);
				}
			}
		}

		number_of_border_vertices = vector_of_border_points.size();

		return vector_of_border_points;
	}


	// (DEPRECATED) this function will flag on the input image the vertices on the border 
	CImg <unsigned char> image_with_vertex_on_border(const CImg <unsigned char>& image) 
	{
		CImg <unsigned char> result = image;

		Finite_vertices_iterator fv = this->finite_vertices_begin();
		for(; fv != this->finite_vertices_end(); ++fv)
		{
			if(fv->info().get_flag())
			{ //flag is in Alpha
				result((unsigned int)fv->point().x(),(unsigned int)fv->point().y(),0,3)= 255;
			}
		}
		return result;
	}


	// (DEPRECATED)this function will draw segments between all border points
	std::vector<Segment_2> link_points_on_border(const std::vector<Vertex_handle>& vector_of_border_points)
	{


		std::vector<Segment_2> vector_of_segments;

		for(int i=0; i<vector_of_border_points.size()-1; ++i)
		{
			for(int j=i+1; j<vector_of_border_points.size(); ++j)
			{
				Segment_2 segment(vector_of_border_points[i]->point(),vector_of_border_points[j]->point());			
				vector_of_segments.push_back(segment);
			}
		}


		return vector_of_segments;
	}


	// routine to "navigate" in the segment
	Point_2 point_on_segment(Vector_2 v, Point_2 source, double u)
	{
		return source + u*v;
	}


	// add constraints in triangulation 
	// adaptation of the Needleman-Wunsch algorithm for global alignement
	// paper reference:
	// "A general method applicable to the search for similarities in the amino acid sequence of two proteins"
	std::vector<Segment_2> nw_add_constraints(const CImg <unsigned char>& image, int gap_score)
	{

		std::vector<Vertex_handle>   vector_of_border_points = set_on_border_2D_vertices(image);

		Gyroviz_border_points_dt2<Gt, Tds> border_dt2(vector_of_border_points);

		std::vector<Segment_2/*Vertex_handle*/> vector_of_vertex_segments = border_dt2.segments_out_of_dt2();

		Point_2 source_vertex;
		Point_2 end_vertex;

		Vertex_handle source_handle;
		Vertex_handle end_handle;
		Vector_2 v;
		std::vector<Segment_2> vector_of_constraints;

		double length_curr_segment = 0;
		int global_score = -1;// no constraint

		// the similarity matrix i will use is 
		// S(Border/Segment) = +1
		// S(Gap/Segment) = gap_score (default : -2)

		for(int i = 0; i</*vector_of_vertex_segments*/vector_of_vertex_segments.size()/*-1*/; ++i)
		{

			//source_handle = vector_of_vertex_segments[i];
			//end_handle = vector_of_vertex_segments[++i];
			//source_vertex = source_handle->point();
			//end_vertex    = end_handle->point();
			//Segment_2 s = (source_vertex, end_vertex);

			source_vertex =  vector_of_vertex_segments[i].source();
			end_vertex =  vector_of_vertex_segments[i].target();
			Segment_2 s = vector_of_vertex_segments[i];
			v = end_vertex - source_vertex;

			length_curr_segment = (int)ceil(sqrt(s.squared_length()));
			int current_gap_score = gap_score;

			for(int j = 0; j<length_curr_segment; ++j)
			{ //TO OPTIMIZE
				if(image((unsigned int)point_on_segment(v, source_vertex, j/length_curr_segment).x(),
					(unsigned int)point_on_segment(v, source_vertex, j/length_curr_segment).y(), 0, 0) == 255 &&
					image((unsigned int)point_on_segment(v, source_vertex, j/length_curr_segment).x(),
					(unsigned int)point_on_segment(v, source_vertex, j/length_curr_segment).y(), 0, 1) == 255 &&
					image((unsigned int)point_on_segment(v, source_vertex, j/length_curr_segment).x(),
					(unsigned int)point_on_segment(v, source_vertex, j/length_curr_segment).y(), 0, 2) == 255)
				{
					global_score = global_score + current_gap_score;
					current_gap_score = gap_score;
				}
				else
				{
					current_gap_score++;
					global_score = global_score - current_gap_score/*/gap_score*/; //TEST

				}
			}

			if(global_score >= 0)
			{
				this->insert_constraint(source_vertex,end_vertex/*source_handle,end_handle*/);
				vector_of_constraints.push_back(s);
				global_score = 0;
			}
			else
				global_score = 0;
		}

		return vector_of_constraints;
	}


	// add constraints in triangulation 
	// adaptation of the Needleman-Wunsch algorithm for local alignement
	// paper reference:
	// "Identification of common molecular subsequences"
	void sw_add_constraints()
	{
		Point_2 origin_vertex;
		Point_2 end_vertex;
	}



	static Gyroviz_constrained_triangle_soup list_cdt2_to_vector_twc(const std::list<Gyroviz_cdt2>& list_cdt2)
	{
		Gyroviz_constrained_triangle_soup result;
		
		int counter = 0; // cam index starts at one


		for (std::list<Gyroviz_cdt2>::const_iterator cdt2_it = list_cdt2.begin(); 
			cdt2_it != list_cdt2.end(); ++cdt2_it)
		{
			Finite_faces_iterator ff = cdt2_it->finite_faces_begin();
			for(; ff != cdt2_it->finite_faces_end(); ++ff)
			{
				counter++;
				
				Vertex_handle v1 = ff->vertex(0);
				Vertex_handle v2 = ff->vertex(1);
				Vertex_handle v3 = ff->vertex(2);
				Point_3 p1 = v1->info().get_point3();
				Point_3 p2 = v2->info().get_point3();
				Point_3 p3 = v3->info().get_point3();
				Triangle_3 t(p1,p2,p3);

				Gyroviz_triangle_with_cam triangle_wc(t,counter);
				
				result.push_back(triangle_wc);
			}
		}
		return result;
	}


	/// Get the region of interest, ignoring the outliers.
	/// This method is used to define the OpenGL arcball sphere.
	Sphere region_of_interest() const
	{
		// A good candidate is a sphere containing the dense region of the point cloud:
		// - center point is barycenter
		// - Radius is 2 * standard deviation
		float radius = 2.f * (float)m_standard_deviation;
		return Sphere(m_barycenter, radius*radius);
	}

	/// Update region of interest.
	/// Owner is responsible to call this function after modifying the triangulation.
	void update_bounding_box()
	{
		// Update bounding box and barycenter.
		// TODO: we should use the functions in PCA component instead.
		FT xmin,xmax,ymin,ymax,zmin,zmax;
		xmin = ymin = zmin =  1e38;
		xmax = ymax = zmax = -1e38;
		Vector_3 v = CGAL::NULL_VECTOR;
		FT norm = 0;
		assert(points_begin() != points_end());
		Finite_vertices_iterator fv = this->finite_vertices_begin();
		for(; fv != this->finite_vertices_end(); ++fv)
		{
			const Point_3& p = fv->info().get_point3();

			// update bbox
			xmin = (std::min)(p.x(),xmin);
			ymin = (std::min)(p.y(),ymin);
			zmin = (std::min)(p.z(),zmin);
			xmax = (std::max)(p.x(),xmax);
			ymax = (std::max)(p.y(),ymax);
			zmax = (std::max)(p.z(),zmax);

			// update barycenter
			v = v + (p - CGAL::ORIGIN);
			norm += 1;
		}
		//
		Point_3 p(xmin,ymin,zmin);
		Point_3 q(xmax,ymax,zmax);
		m_bounding_box = Iso_cuboid_3(p,q);
		//
		m_barycenter = CGAL::ORIGIN + v / norm;

		/// Compute standard deviation
		Geom_traits::Compute_squared_distance_3 sqd;
		FT sq_radius = 0;
		/*Finite_vertices_iterator*/ fv = this->finite_vertices_begin();
		for(; fv != this->finite_vertices_end(); ++fv)
		{
			const Point_3& p = fv->info().get_point3();
			sq_radius += sqd(p, m_barycenter);
		}
		sq_radius /= number_of_vertices();
		m_standard_deviation = CGAL::sqrt(sq_radius);
	}


	// draw 2D cdt vertices
	void gl_draw_2D_vertices(const unsigned char r, const unsigned char g,
		const unsigned char b, float size)
	{
		//DEBUG
		int number_of_vertices_2v = 0;

		::glPointSize(size);
		::glColor3ub(r,g,b);
		::glBegin(GL_POINTS);

		Finite_vertices_iterator fv = this->finite_vertices_begin();
		for(; fv != this->finite_vertices_end(); ++fv)
		{
			::glVertex2d(fv->point().x(),fv->point().y());


			number_of_vertices_2v++;
		}

		::glEnd();
	}


	// draw 2D only points near detected borders 
	void gl_draw_on_border_2D_vertices(const unsigned char r, const unsigned char g,
		const unsigned char b, float size/* TEST, const CImg <unsigned char>& image*/)
	{
		/*
		TEST 

		::glPointSize(size);
		::glColor3ub(r,g,b);
		::glBegin(GL_POINTS);

		std::vector<Point_2> vector_of_border_points = set_on_border_2D_vertices(image);
		for(int i=0; i<vector_of_border_points.size(); ++i)
		{
		::glVertex2d(vector_of_border_points[i].x(),vector_of_border_points[i].y());
		}
		::glEnd();*/

		//DEBUG
		int number_of_vertices_b = 0;

		::glPointSize(size);
		::glColor3ub(r,g,b);
		::glBegin(GL_POINTS);

		Finite_vertices_iterator fv  = this->finite_vertices_begin();
		for(; fv != this->finite_vertices_end(); ++fv)
		{
			if(fv->info().get_flag()){
				::glVertex2d(fv->point().x(),fv->point().y());

				number_of_vertices_b++;

			}		
		}

		::glEnd();

	}


	// draw 2D cdt constrained edges
	void gl_draw_2D_constrained_edges(const unsigned char r, const unsigned char g,
		const unsigned char b, float line_width/*, const CImg <unsigned char>& image, int gap_score*/)
	{
		::glColor3ub(r,g,b);
		::glLineWidth(line_width);
		::glBegin(GL_LINES);

		// DEBUG 
		int counter_constraints=0;

		//std::vector<Segment_2> vector_of_constraints = nw_add_constraints(image, gap_score);
		Finite_edges_iterator fe = this->finite_edges_begin();
		for(; fe != this->finite_edges_end(); ++fe)
		{
			if(fe->first->is_constrained(fe->second))
			{

				counter_constraints++;

				Point_2 p1 = fe->first->vertex(ccw(fe->second))->point();
				Point_2 p2 = fe->first->vertex(cw(fe->second))->point();
				::glVertex2d(p1.x(), p1.y());
				::glVertex2d(p2.x(), p2.y());
			}

		}
		//for(int i=0; i<vector_of_constraints.size(); ++i)
		//{
		//	::glVertex2d(vector_of_constraints[i].source().x(),vector_of_constraints[i].source().y());
		//	::glVertex2d(vector_of_constraints[i].target().x(),vector_of_constraints[i].target().y());
		//}
		::glEnd();

	}


	// draw 2D cdt constrained delaunay triangles
	void gl_draw_2D_constrained_delaunay_triangles(const unsigned char r, const unsigned char g,
		const unsigned char b, float line_width /*, CImg <unsigned char> image, int gap_score*/)
	{
		//::glColor3ub(r,g,b);
		//::glBegin(GL_TRIANGLES);
		//Finite_faces_iterator ff = this->finite_faces_begin();
		//for(; ff != this->finite_faces_end(); ++ff)
		//{
		//	Vertex_handle v1 = ff->vertex(0);
		//	Vertex_handle v2 = ff->vertex(1);
		//	Vertex_handle v3 = ff->vertex(2);
		//	::glVertex2d(v1->point().x(), v1->point().y());
		//	::glVertex2d(v2->point().x(), v2->point().y());
		//	::glVertex2d(v3->point().x(), v3->point().y());
		//} 
		//::glEnd();

		::glLineWidth(line_width);
		::glBegin(GL_LINES);

		Finite_edges_iterator fe = this->finite_edges_begin();
		for(; fe != this->finite_edges_end(); ++fe)
		{
			if(fe->first->is_constrained(fe->second))
				::glColor3ub(255,0,0);
			else
				::glColor3ub(r,g,b);

			Point_2 p1 = fe->first->vertex(ccw(fe->second))->point();
			Point_2 p2 = fe->first->vertex(cw(fe->second))->point();
			::glVertex2d(p1.x(), p1.y());
			::glVertex2d(p2.x(), p2.y());
		}
		::glEnd();
	} 


	// 3D projection of the tracked 2D points
	void gl_draw_soup_vertices(const unsigned char r, const unsigned char g,
		const unsigned char b, float size)
	{
		//DEBUG
		int number_of_vertices_sv = 0;

		::glPointSize(size);
		::glColor3ub(r,g,b);
		::glBegin(GL_POINTS);

		Finite_vertices_iterator fv = this->finite_vertices_begin();
		for(; fv != this->finite_vertices_end(); ++fv)
		{
			const Point_3& p = fv->info().get_point3();
			::glVertex3d(p.x(),p.y(),p.z());

			number_of_vertices_sv++;

		}
		::glEnd();
	}


	// draw 2D only points near detected borders 
	void gl_draw_on_border_3D_vertices(const unsigned char r, const unsigned char g,
		const unsigned char b, float size)
	{
		//DEBUG
		int number_of_vertices_bv = 0;

		::glPointSize(size);
		::glColor3ub(r,g,b);
		::glBegin(GL_POINTS);

		Finite_vertices_iterator fv  = this->finite_vertices_begin();
		for(; fv != this->finite_vertices_end(); ++fv)
		{
			if(fv->info().get_flag())
			{
				::glVertex3d(fv->info().get_point3().x(),fv->info().get_point3().y(),fv->info().get_point3().z());

				number_of_vertices_bv++;

			}		
		}
		::glEnd();
	}



	// 3D projection of the tracked 2D constrained edges
	void gl_draw_soup_constrained_edges(const unsigned char r, const unsigned char g,
		const unsigned char b, const float width/*, const CImg <unsigned char>& image*/)
	{
		::glLineWidth(width);
		::glColor3ub(r,g,b);
		::glBegin(GL_LINES);

		/*std::vector<Segment_2> vector_of_constraints = nw_add_constraints(image);*/
		Finite_edges_iterator fe = this->finite_edges_begin();
		for(; fe != this->finite_edges_end(); ++fe)
		{
			if(fe->first->is_constrained(fe->second))
			{
				Point_3 p1 = fe->first->vertex(ccw(fe->second))->info().get_point3();
				Point_3 p2 = fe->first->vertex(cw(fe->second))->info().get_point3();
				::glVertex3d(p1.x(), p1.y(),  p1.z());
				::glVertex3d(p2.x(), p2.y(),  p2.z());
			}

		}

		::glEnd();
	}


	// 3D projection of the tracked 2D constrained triangulation
	void gl_draw_soup_constrained_triangles(const unsigned char r, const unsigned char g,
		const unsigned char b){

			::glColor3ub(r,g,b);
			::glBegin(GL_TRIANGLES);

			Finite_faces_iterator ff = this->finite_faces_begin();
			for(; ff != this->finite_faces_end(); ++ff)
			{
				Vertex_handle v1 = ff->vertex(0);
				Vertex_handle v2 = ff->vertex(1);
				Vertex_handle v3 = ff->vertex(2);
				Point_3 p1 = v1->info().get_point3();
				Point_3 p2 = v2->info().get_point3();
				Point_3 p3 = v3->info().get_point3();
				::glVertex3d(p1.x(), p1.y(), p1.z());
				::glVertex3d(p2.x(), p2.y(), p2.z());
				::glVertex3d(p3.x(), p3.y(), p3.z());
			} 
			::glEnd();
	} 



};

#endif // _Gyroviz_cdt2_
