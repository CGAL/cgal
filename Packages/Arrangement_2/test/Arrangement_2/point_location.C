//#include "short_names.h"
//#include <CGAL/MP_Float.h>
//#include <CGAL/Quotient.h>
//#include <CGAL/Arr_observer.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_triangulation_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_point_location/Arr_lm_random_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_grid_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_halton_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_middle_edges_generator.h>

#include "Segment_reader.h"

typedef CGAL::Quotient<CGAL::MP_Float>                    Number_type;
typedef CGAL::Cartesian<Number_type>                      Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                Traits_2;
typedef Traits_2::Point_2                                 Point_2;
typedef Traits_2::Curve_2                                 Curve_2;
typedef Traits_2::X_monotone_curve_2                      Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                     Arrangement_2;
typedef Arrangement_2::Halfedge_handle                    Halfedge_handle;

typedef CGAL::Arr_naive_point_location<Arrangement_2>     
                                                    Naive_point_location;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> 
                                                    Walk_point_location;
typedef CGAL::Arr_triangulation_point_location<Arrangement_2> 
                                                    Trg_point_location;
typedef CGAL::Arr_trapezoid_ric_point_location<Arrangement_2> 
                                                    Ric_point_location;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2> 
                                                    Lm_point_location;
typedef CGAL::Arr_random_landmarks_generator<Arrangement_2>
                                                    Random_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Random_lm_generator> 
                                                    Lm_random_point_location;
typedef CGAL::Arr_grid_landmarks_generator<Arrangement_2>
                                                    Grid_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Grid_lm_generator> 
                                                    Lm_grid_point_location;
typedef CGAL::Arr_halton_landmarks_generator<Arrangement_2>
                                                    Halton_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Halton_lm_generator> 
                                                    Lm_halton_point_location;
typedef CGAL::Arr_middle_edges_landmarks_generator<Arrangement_2>
                                                    Mide_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Mide_lm_generator> 
                                                    Lm_mide_point_location;
// ===> Add new point location type here <===

typedef std::list<Point_2>                                Points_list;
typedef Points_list::iterator                             Point_iterator;
typedef std::list<Curve_2>                                Curve_list;
typedef std::vector<CGAL::Object>                         Objects_vector;
typedef Objects_vector::iterator                          Object_iterator;

// ===> Change the number of point-location startegies
//      when a new point location is added. <===
#define NUM_OF_POINT_LOCATION_STRATEGIES 9


int check_point_location (Arrangement_2 &arr, Points_list &plist)
{
  //init - all point locations
  Naive_point_location            naive_pl (arr);                 // 0
  Walk_point_location             walk_pl (arr);                  // 1
  Trg_point_location              trg_pl(arr);                    // 2
  Ric_point_location              ric_pl(arr);                    // 3
  Lm_point_location               lm_pl (arr);                    // 4
  Random_lm_generator             random_g(arr);    
  Lm_random_point_location        randon_lm_pl (arr, &random_g);  // 5
  Grid_lm_generator               grid_g(arr);      
  Lm_grid_point_location          grid_lm_pl (arr, &grid_g);      // 6
  Halton_lm_generator             halton_g(arr);
  Lm_halton_point_location        halton_lm_pl (arr, &halton_g);  // 7
  Mide_lm_generator               mide_g(arr);
  Lm_mide_point_location          mide_lm_pl (arr, &mide_g);      // 8
  // ===> Add new point location instance here. <===

  CGAL::Object                    obj;
  Objects_vector                  objs[NUM_OF_POINT_LOCATION_STRATEGIES];
  Object_iterator                 ob_iter[NUM_OF_POINT_LOCATION_STRATEGIES];
  Arrangement_2::Vertex_const_handle    vh[NUM_OF_POINT_LOCATION_STRATEGIES];
  Arrangement_2::Halfedge_const_handle  hh[NUM_OF_POINT_LOCATION_STRATEGIES];
  Arrangement_2::Face_const_handle      fh[NUM_OF_POINT_LOCATION_STRATEGIES];
  
  Point_2               q;
  Point_iterator        piter;

  //LOCATE the points in the list using all PL strategies

  std::cout << "Naive" << std::endl;
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = naive_pl.locate (q);
    objs[0].push_back(obj);
  }

  std::cout << "Walk" << std::endl;
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = walk_pl.locate (q);
    objs[1].push_back(obj);
  }

  std::cout << "Triangulation" << std::endl;
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = trg_pl.locate (q);
    objs[2].push_back(obj);
  }

  std::cout << "RIC" << std::endl;
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = ric_pl.locate (q);
    objs[3].push_back(obj);
  }

  std::cout << "Landmarks (vertices)" << std::endl;
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = lm_pl.locate (q);
    objs[4].push_back(obj);
  }

  std::cout << "Random LM" << std::endl;
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = randon_lm_pl.locate (q);
    objs[5].push_back(obj);
  }

  std::cout << "Grid LM" << std::endl;
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = grid_lm_pl.locate (q);
    objs[6].push_back(obj);
  }

  std::cout << "Halton LM" << std::endl;
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = halton_lm_pl.locate (q);
    objs[7].push_back(obj);
  }

  std::cout << "Middle edges LM" << std::endl;
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = mide_lm_pl.locate (q);
    objs[8].push_back(obj);
  }

  // ===> Add a call to operate the the new point location. <===

  //END LOCATION 
  int pls_num = NUM_OF_POINT_LOCATION_STRATEGIES;
  int pl_index;

  //init all obejct iterators
  for (pl_index=0; pl_index<pls_num; pl_index++)
  {
    ob_iter[pl_index] = objs[pl_index].begin();
  }
   
  //get size of objects
  unsigned int size = objs[0].size();
  std::cout <<"size is "<< size << std::endl;

  for (pl_index=0; pl_index<pls_num; pl_index++)
  {
    if (size != objs[pl_index].size())
    {
      std::cout << "Error: size of pl number "<<pl_index<<" is "
        <<objs[pl_index].size()<< std::endl;
      return (-1);
    }
  }

  //assign and check results
  unsigned int qi; //qi is the query point index
  for (qi=0; qi<size; qi++) 
  {
    //assign object to a face
    if (CGAL::assign (fh[0], ob_iter[0][qi]))
    {    
      for (int pl_index=1; pl_index<pls_num; pl_index++)
      {
        if (! CGAL::assign(fh[pl_index], ob_iter[pl_index][qi]))
        {
          std::cout << "Error in point location number " 
            << pl_index << "did no return a face"<< std::endl;
          return (-1);
        }
        else if (fh[pl_index] != fh[0])
        {
          std::cout << "Error: point location number " 
            << pl_index << "return a different face"<< std::endl;
          return (-1);         
        }
      }  
      //if (fh[0]->is_unbounded())
      //  std::cout << "Unbounded face." << std::endl;
      //else
      //  std::cout << "Face." << std::endl;
    }

    //assign object to a halfedge
    else if (CGAL::assign (hh[0], ob_iter[0][qi]))
    {
      for (int pl_index=1; pl_index<pls_num; pl_index++)
      {
        if (! CGAL::assign(hh[pl_index], ob_iter[pl_index][qi]))
        {
          std::cout << "Error in point location number " 
            << pl_index << "did no return halfedge"<< std::endl;
          return (-1);
        }
        else if (hh[pl_index] != hh[0])
        {
          std::cout << "Error: point location number " 
            << pl_index << "return a different halfedge"<< std::endl;
          return (-1);         
        }
      }
      //std::cout << "Halfedge: "<< hh[0]->curve() << std::endl;
    }

    //assign object to a vertex
    else if (CGAL::assign (vh[0], ob_iter[0][qi]))
    {
      for (int pl_index=1; pl_index<pls_num; pl_index++)
      {
        if (! CGAL::assign(vh[pl_index], ob_iter[pl_index][qi]))
        {
          std::cout << "Error in point location number " 
            << pl_index << "did no return a vertex"<< std::endl;
          return (-1);
        }
        else if (vh[pl_index] != vh[0])
        {
          std::cout << "Error: point location number " 
            << pl_index << "return a different vertex"<< std::endl;
          return (-1);         
        }
      }
      //std::cout << "Vertex: "<< vh[0]->point() << std::endl;
    }
   
    else
    {
      std::cout << "Illegal point-location result." << std::endl;    
      return (-1);
    }
  }
  
  return (0);
}

/*! */
int read_points(const char * points_filename, Points_list &plist)
{
  //read points from file into list
  std::ifstream inp_pnt_file(points_filename);
  if (!inp_pnt_file.is_open()) 
  {
    std::cerr << "Cannot open file " << points_filename << "!" << std::endl;
    return (-1);
  }

  int points_count;
  inp_pnt_file >> points_count;

  for (int i = 0; i < points_count; i++) {
    Number_type x, y;
    inp_pnt_file >> x >> y;
    Point_2 pnt(x, y);
    plist.push_back(pnt);
  }

  return 0;
}

int main (int argc, char * argv[])
{
 	//get arguments
	if (argc < 3) {
		std::cout << "Usage: " << argv[0] <<" curve_file pnt_file" << std::endl;
		std::cout << "curve_file  - the input curves file" << std::endl;
		std::cout << "pnt_file    - the input query points" << std::endl;
		exit(-1);
	}	

	const char * curves_filename = argv[1];
	const char * points_filename = argv[2];

  //read curves and insert them into the arrangement 
	Segment_reader<Traits_2> reader;
  CGAL::Bbox_2             bbox;
	Curve_list               curve_list;
	reader.read_data(curves_filename, std::back_inserter(curve_list), bbox);	

  //insert all curves into the arrangement
  Arrangement_2  arr;
  insert (arr, curve_list.begin(), curve_list.end());

  //read point and insert them into a list of points
  Points_list           plist;

  //read points from file into list
  if (read_points(points_filename, plist))
  {
    std::cerr << "ERROR in read_points."<<std::endl;
    exit(-1);
  }

  //check point location of points
  if (check_point_location(arr, plist))
  {
    std::cerr << "ERROR in check_point_location."<<std::endl;
    exit(-1);
  }
  
  std::cout << "Point-location test was finished O.K."<<std::endl;

  return (0);
}
