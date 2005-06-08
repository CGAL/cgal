// file: examples/Arrangement_2/example11.C


//#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_point_location/Arr_lm_random_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_grid_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_halton_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_middle_edges_generator.h>

typedef CGAL::Quotient<CGAL::MP_Float>					Number_type;
typedef CGAL::Cartesian<Number_type>					Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>				Traits_2;
typedef Traits_2::Point_2								Point_2;
typedef Traits_2::X_monotone_curve_2					Segment_2;
typedef CGAL::Arrangement_2<Traits_2>					Arrangement_2;
typedef Arrangement_2::Halfedge_handle					Halfedge_handle;
typedef CGAL::Arr_naive_point_location<Arrangement_2>	Naive_point_location;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> 
														Walk_point_location;
//typedef CGAL::Arr_random_landmarks_generator<Arrangement_2> Landmarks_generator;
//typedef CGAL::Arr_grid_landmarks_generator<Arrangement_2> Landmarks_generator;
//typedef CGAL::Arr_halton_landmarks_generator<Arrangement_2> Landmarks_generator;
typedef CGAL::Arr_middle_edges_landmarks_generator<Arrangement_2>
														Landmarks_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Landmarks_generator> 
														Lm_point_location;


typedef std::list<Point_2>		Points_list;
typedef Points_list::iterator	Point_iterator;


class My_observer : public CGAL::Arr_observer<Arrangement_2>
{
public:

	My_observer (Arrangement_2& arr) :
	  CGAL::Arr_observer<Arrangement_2> (arr)
	  {}

	  virtual void before_split_face (Face_handle,
		  Halfedge_handle e)
	  {
		  std::cout << "-> New face, caused by: " << e.curve() << std::endl;
	  }

	  virtual void after_modify_edge (Halfedge_handle e)
	  {
		  std::cout << "-> Existing edge " << e.curve() 
			  << " has just been modified." << std::endl;
	  }

	  virtual void after_assign ()
	  { 
		  std::cout << "after_assign"<< std::endl;
	  }

	  virtual void after_clear (Face_handle /* u */)
	  { 
		  std::cout << "after_clear"<< std::endl;
	  }

	  virtual void before_global_change ()
	  { 
		  std::cout << "before_global_change"<< std::endl;
	  }

	  virtual void after_global_change ()
	  {
		  std::cout << "after_global_change"<< std::endl;
	  }

	  virtual void after_create_vertex (Vertex_handle /* v */)
	  {
		  std::cout << "after_create_vertex"<< std::endl;
	  }

	  virtual void after_remove_vertex ()
	  {
		  std::cout << "after_remove_vertex"<< std::endl;
	  }
};

int main ()
{
	// Construct the arrangement using the specialized insertion functions.
	Arrangement_2  arr;
	My_observer    obs (arr);

	Segment_2       cv1 (Point_2(1.0, 0.0), Point_2(3.0, 2.0));
	Segment_2       cv2 (Point_2(4.0, -1.0), Point_2(3.0, -2.0));
	Segment_2       cv3 (Point_2(4.0, -1.0), Point_2(1.0, 0.0));
	Segment_2       cv4 (Point_2(1.0, 0.0), Point_2(4.0, 1.0));
	Segment_2       cv5 (Point_2(3.0, 2.0), Point_2(4.0, 1.0));
	Segment_2       cv6 (Point_2(6.0, 0.0), Point_2(4.0, -1.0));
	Segment_2       cv7 (Point_2(4.0, 1.0), Point_2(6.0, 0.0));

	Halfedge_handle h1 = arr.insert_in_face_interior (cv1, arr.unbounded_face());
	Halfedge_handle h2 = arr.insert_in_face_interior (cv2, arr.unbounded_face());
	Halfedge_handle h3 = arr.insert_at_vertices (cv3, h2, h1.twin());
	Halfedge_handle h4 = arr.insert_from_left_vertex (cv4, h1.twin());
	Halfedge_handle h5 = arr.insert_at_vertices (cv5, h1, h4);
	Halfedge_handle h6 = arr.insert_from_left_vertex (cv6, h3.twin());
	Halfedge_handle h7 = arr.insert_at_vertices(cv7, h5, h6);

	// Print the arrangement vertices.
	Arrangement_2::Vertex_const_iterator  vit;
	Arrangement_2::Vertex_const_handle    vh;
	int                                   i, j;

	std::cout << arr.number_of_vertices() << " vertices:" << std::endl;
	for (i = 1, vit = arr.vertices_begin(); 
		vit != arr.vertices_end(); vit++, i++)
	{
		vh = *vit;
		std::cout << '\t' << i << ": " << vh.point() << std::endl;
	}
	std::cout << std::endl;

	// Print the arrangement edges.
	Arrangement_2::Edge_const_iterator    eit;
	Arrangement_2::Halfedge_const_handle  hh;

	std::cout << arr.number_of_edges() << " edges:" << std::endl;
	for (i = 1, eit = arr.edges_begin(); eit != arr.edges_end(); eit++, i++)
	{
		hh = *eit;
		std::cout << '\t' << i << ": " << hh.curve() << std::endl;
	}
	std::cout << std::endl;

	// Print the arrangement faces.
	Arrangement_2::Face_const_iterator           fit;
	Arrangement_2::Face_const_handle             fh;
	Arrangement_2::Ccb_halfedge_const_circulator ccb;
	Arrangement_2::Holes_const_iterator          hoit;

	std::cout << arr.number_of_faces() << " faces." << std::endl;
	for (i = 1, fit = arr.faces_begin(); fit != arr.faces_end(); fit++, i++)
	{
		// Print the outer boundary of the face.
		fh = *fit;
		std::cout << '\t' << i << ": ";
		if (fh.is_unbounded())
		{
			std::cout << "Unbounded face." << std::endl;
		}
		else
		{
			ccb = fh.outer_ccb();
			std::cout << (*ccb).source().point();
			do
			{
				std::cout << " -> " << (*ccb).target().point();
				ccb++;
			} while (ccb != fh.outer_ccb());
			std::cout << std::endl;
		}

		// Print the holes.
		for (j = 1, hoit = fh.holes_begin(); hoit != fh.holes_end(); hoit++, j++)
		{
			std::cout << "\t\tHole " << i << ": ";
			ccb = *hoit;
			std::cout << (*ccb).source().point();
			do
			{
				std::cout << " -> " << (*ccb).target().point();
				ccb++;
			} while (ccb != *hoit);
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;

	// Perform point location.
	//Lm_point_location		lm_pl (arr);
	
	std::cout << " before creating Landmarks_generator" << std::endl;
	Landmarks_generator		lg(arr);
	std::cout << " after creating Landmarks_generator" << std::endl;
	Lm_point_location		lm_pl (arr, &lg);
	std::cout << " after creating Lm_point_location" << std::endl;

	Walk_point_location		walk_pl (arr);
	Naive_point_location	naive_pl (arr);
	CGAL::Object			lm_obj, walk_obj, naive_obj;
	Arrangement_2::Vertex_const_handle    lm_vh, walk_vh, naive_vh;
	Arrangement_2::Halfedge_const_handle  lm_hh, walk_hh, naive_hh;
	Arrangement_2::Face_const_handle      lm_fh, walk_fh, naive_fh;

	Point_2			q;
	Points_list		plist;
	Point_iterator	piter;

	plist.push_back(Point_2(0,0));
	plist.push_back(Point_2(1,0));
	plist.push_back(Point_2(3,2));
	plist.push_back(Point_2(5,0.5));
	plist.push_back(Point_2(4,1));
	plist.push_back(Point_2(3,1));
	plist.push_back(Point_2(3,-2));
	plist.push_back(Point_2(3,-3));
	plist.push_back(Point_2(2,-1));
	plist.push_back(Point_2(3,0));
	plist.push_back(Point_2(2,-2));
	plist.push_back(Point_2(-1,-1));
	plist.push_back(Point_2(-1,0));
	plist.push_back(Point_2(3,-1));

	for (piter = plist.begin(); piter != plist.end(); piter++ )
	{
		q = (*piter);

		std::cout << "Before locate. q = " << q << std::endl;
		lm_obj = lm_pl.locate (q);
		walk_obj = walk_pl.locate (q);
		naive_obj = naive_pl.locate (q);
		std::cout << "After locate." << std::endl;

		if (CGAL::assign (lm_fh, lm_obj))
		{    
			if (lm_fh.is_unbounded())
				std::cout << "Inside unbounded face." << std::endl;
			else
				std::cout << "Inside face." << std::endl;

			if (CGAL::assign (walk_fh, walk_obj))
			{
				if (CGAL::assign (naive_fh, naive_obj))
				{
					if ((lm_fh == naive_fh) &&	(lm_fh == walk_fh))
					{
						std::cout 
							<< "All point locations returned the same face"
							<< std::endl;
					}
					else if (naive_fh == walk_fh)
						std::cout << "Error in Lm pl." << std::endl;
					else if (naive_fh == lm_fh)
						std::cout << "Error in walk pl." << std::endl;
					else
						std::cout << "Error in naive pl." << std::endl;
				}
				else
				{
					std::cout << "Error: naive pl does not return the same value" 
						<< std::endl;
					CGAL_assertion(false);
				}
			}
			else
			{
				std::cout << "Error: walk pl does not return the same value" 
					<< std::endl;
				CGAL_assertion(false);
			}
		}
		else if (CGAL::assign (lm_hh, lm_obj))
		{
			std::cout << "On halfedge: " << lm_hh.curve() << std::endl;
			if (CGAL::assign (walk_hh, walk_obj))
			{
				if (CGAL::assign (naive_hh, naive_obj))
				{
					if ((lm_hh == naive_hh) &&	(lm_hh == walk_hh))
						std::cout 
							<< "All point locations returned the same edge"
							<< std::endl;
					else if (naive_hh == walk_hh)
						std::cout << "Error in Lm pl." << std::endl;
					else if (naive_hh == lm_hh)
						std::cout << "Error in walk pl." << std::endl;
					else
						std::cout << "Error in naive pl." << std::endl;				
				}
				else
				{
					std::cout << "Error: naive pl does not return the same value" 
						<< std::endl;
					CGAL_assertion(false);
				}
			}
			else
			{
				std::cout << "Error: walk pl does not return the same value" 
					<< std::endl;
				CGAL_assertion(false);
			}
		}
		else if (CGAL::assign (lm_vh, lm_obj))
		{
			std::cout << "On vertex: " << lm_vh.point() << std::endl;
			if (CGAL::assign (walk_vh, walk_obj))
			{
				if (CGAL::assign (naive_vh, naive_obj))
				{
					if ((lm_vh == naive_vh) &&	(lm_vh == walk_vh))
						std::cout 
							<< "All point locations returned the same edge"
							<< std::endl;
					else if (naive_vh == walk_vh)
						std::cout << "Error in Lm pl." << std::endl;
					else if (naive_vh == lm_vh)
						std::cout << "Error in walk pl." << std::endl;
					else
						std::cout << "Error in naive pl." << std::endl;				
				}
				else
				{
					std::cout << "Error: naive pl does not return the same value" 
						<< std::endl;
					CGAL_assertion(false);
				}
			}
			else
			{
				std::cout << "Error: walk pl does not return the same value" 
					<< std::endl;
				CGAL_assertion(false);
			}
		}
		else
		{
			std::cout << "Illegal point-location result." << std::endl;    
			CGAL_assertion(false);
		}
	}


	// Insert additional segments.
	Segment_2       s1 (Point_2(-1.0, 0.0), Point_2(0.0, 2.0));
	Segment_2       s2 (Point_2(-1.0, 0.0), Point_2(-1.0, -2.0));
	Segment_2       s3 (Point_2(0.0, 2.0), Point_2(1.0, 0.0));
	Segment_2       s4 (Point_2(-1.0, -2.0), Point_2(3.0, -2.0));

	arr_insert_non_intersecting (arr, lm_pl, s1);
	arr_insert_non_intersecting (arr, lm_pl, s2);
	arr_insert_non_intersecting (arr, lm_pl, s3);
	arr_insert_non_intersecting (arr, lm_pl, s4);

	//// Perform point location again.
	//obj = pl.locate (q);
	//if (CGAL::assign (fh, obj))
	//{
	//	if (fh.is_unbounded())
	//		std::cout << "Inside unbounded face." << std::endl;
	//	else
	//		std::cout << "Inside face." << std::endl;
	//}
	//else if (CGAL::assign (hh, obj))
	//{
	//	std::cout << "On halfedge: " << hh.curve() << std::endl;
	//}
	//else if (CGAL::assign (vh, obj))
	//{
	//	std::cout << "On vertex: " << vh.point() << std::endl;
	//}
	//else
	//{
	//	std::cout << "Illegal point-location result." << std::endl;    
	//}

	// Test the insertion function (iis2 and iis3 cause some overlaps).
	Segment_2       iis1 (Point_2(0.0, -3.0), Point_2(5.0, 2.0));

	arr_insert (arr, lm_pl, iis1);

	std::cout << "V = " << arr.number_of_vertices()
		<< ",  E = " << arr.number_of_edges() 
		<< ",  F = " << arr.number_of_faces() << std::endl;

	Segment_2       iis2 (Point_2(-0.0, -2.0), Point_2(5.0, -2.0));

	arr_insert (arr, lm_pl, iis2);

	std::cout << "V = " << arr.number_of_vertices()
		<< ",  E = " << arr.number_of_edges() 
		<< ",  F = " << arr.number_of_faces() << std::endl;

	Segment_2       iis3 (Point_2(2.0, 2.0), Point_2(5.0, 0.5));

	arr_insert (arr, lm_pl, iis3);

	std::cout << "V = " << arr.number_of_vertices()
		<< ",  E = " << arr.number_of_edges() 
		<< ",  F = " << arr.number_of_faces() << std::endl;

	for (piter = plist.begin(); piter != plist.end(); piter++ )
	{
		q = (*piter);

		std::cout << "Before locate. q = " << q << std::endl;
		lm_obj = lm_pl.locate (q);
		walk_obj = walk_pl.locate (q);
		naive_obj = naive_pl.locate (q);
		std::cout << "After locate." << std::endl;

		if (CGAL::assign (lm_fh, lm_obj))
		{    
			if (lm_fh.is_unbounded())
				std::cout << "Inside unbounded face." << std::endl;
			else
				std::cout << "Inside face." << std::endl;

			if (CGAL::assign (walk_fh, walk_obj))
			{
				if (CGAL::assign (naive_fh, naive_obj))
				{
					if ((lm_fh == naive_fh) &&	(lm_fh == walk_fh))
					{
						std::cout 
							<< "All point locations returned the same face"
							<< std::endl;
					}
					else if (naive_fh == walk_fh)
						std::cout << "Error in Lm pl." << std::endl;
					else if (naive_fh == lm_fh)
						std::cout << "Error in walk pl." << std::endl;
					else
						std::cout << "Error in naive pl." << std::endl;
				}
				else
				{
					std::cout << "Error: naive pl does not return the same value" 
						<< std::endl;
					CGAL_assertion(false);
				}
			}
			else
			{
				std::cout << "Error: walk pl does not return the same value" 
					<< std::endl;
				CGAL_assertion(false);
			}
		}
		else if (CGAL::assign (lm_hh, lm_obj))
		{
			std::cout << "On halfedge: " << lm_hh.curve() << std::endl;
			if (CGAL::assign (walk_hh, walk_obj))
			{
				if (CGAL::assign (naive_hh, naive_obj))
				{
					if ((lm_hh == naive_hh) &&	(lm_hh == walk_hh))
						std::cout 
							<< "All point locations returned the same edge"
							<< std::endl;
					else if (naive_hh == walk_hh)
						std::cout << "Error in Lm pl." << std::endl;
					else if (naive_hh == lm_hh)
						std::cout << "Error in walk pl." << std::endl;
					else
						std::cout << "Error in naive pl." << std::endl;				
				}
				else
				{
					std::cout << "Error: naive pl does not return the same value" 
						<< std::endl;
					CGAL_assertion(false);
				}
			}
			else
			{
				std::cout << "Error: walk pl does not return the same value" 
					<< std::endl;
				CGAL_assertion(false);
			}
		}
		else if (CGAL::assign (lm_vh, lm_obj))
		{
			std::cout << "On vertex: " << lm_vh.point() << std::endl;
			if (CGAL::assign (walk_vh, walk_obj))
			{
				if (CGAL::assign (naive_vh, naive_obj))
				{
					if ((lm_vh == naive_vh) &&	(lm_vh == walk_vh))
						std::cout 
							<< "All point locations returned the same edge"
							<< std::endl;
					else if (naive_vh == walk_vh)
						std::cout << "Error in Lm pl." << std::endl;
					else if (naive_vh == lm_vh)
						std::cout << "Error in walk pl." << std::endl;
					else
						std::cout << "Error in naive pl." << std::endl;				
				}
				else
				{
					std::cout << "Error: naive pl does not return the same value" 
						<< std::endl;
					CGAL_assertion(false);
				}
			}
			else
			{
				std::cout << "Error: walk pl does not return the same value" 
					<< std::endl;
				CGAL_assertion(false);
			}
		}
		else
		{
			std::cout << "Illegal point-location result." << std::endl;    
			CGAL_assertion(false);
		}
	}


	return (0);
}

