
//#define CONICS
#define SEGMENTS
//#define DRAW_QT
//#define SEGMENTS_IN_DOUBLE

#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Timer.h>

#ifdef SEGMENTS
#include <CGAL/Gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#elif defined (CONICS)
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#endif

#ifdef SEGMENTS_IN_DOUBLE
#include "Segment_reader_double.h"
#elif defined (SEGMENTS)
#include "Segment_reader.h"
#elif defined (CONICS)
#include "Conic_reader.h"
#endif

#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_point_location/Arr_lm_random_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_grid_generator.h>
#include <CGAL/Arr_point_location/Arr_lm_halton_generator.h>

#ifdef SEGMENTS
#include <CGAL/Arr_triangulation_point_location.h>
#include <CGAL/Arr_point_location/Arr_lm_middle_edges_generator.h>
#endif

//for Qt windows drawing
#ifdef DRAW_QT
#include <CGAL/IO/Qt_widget.h>
#include <qapplication.h>
typedef CGAL::Qt_widget                                 Window_stream;
#endif

#ifdef SEGMENTS
typedef CGAL::Gmpq                                      Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;

#elif defined (CONICS)
typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Rational                             Number_type;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef Rat_kernel::Point_2                             Rat_point_2;
typedef Rat_kernel::Segment_2                           Rat_segment_2;
typedef Rat_kernel::Circle_2                            Rat_circle_2;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, 
                                 Alg_kernel,
                                 Nt_traits>             Traits_2;
#endif

typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::Curve_2                               Curve_2;
typedef Traits_2::X_monotone_curve_2                    Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;
typedef Arrangement_2::Halfedge_handle                  Halfedge_handle;
typedef Arrangement_2::Edge_const_iterator              Edge_const_iterator;
typedef Arrangement_2::Vertex_const_iterator            Vertex_const_iterator;

typedef CGAL::Arr_naive_point_location<Arrangement_2>     
                                                    Naive_point_location;
typedef CGAL::Arr_walk_along_line_point_location<Arrangement_2> 
                                                    Walk_point_location;
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
#ifdef SEGMENTS
typedef CGAL::Arr_middle_edges_landmarks_generator<Arrangement_2>
                                                    Mide_lm_generator;
typedef CGAL::Arr_landmarks_point_location<Arrangement_2, Mide_lm_generator> 
                                                    Lm_mide_point_location;
typedef CGAL::Arr_triangulation_point_location<Arrangement_2> 
                                                    Trg_point_location;
#endif
// ===> Add new point location type here <===

typedef std::list<Point_2>                                Points_list;
typedef Points_list::iterator                             Point_iterator;
typedef std::list<Curve_2>                                Curve_list;
typedef std::vector<CGAL::Object>                         Objects_vector;
typedef Objects_vector::iterator                          Object_iterator;

// ===> Change the number of point-location startegies
//      when a new point location is added. <===
#ifdef SEGMENTS
#define NUM_OF_POINT_LOCATION_STRATEGIES 9
#else
#define NUM_OF_POINT_LOCATION_STRATEGIES 7
#endif

/*! */
int check_point_location (Arrangement_2 &arr, Points_list &plist)
{
  //init - all point locations
  CGAL::Timer            timer;

   timer.reset(); timer.start();
  Naive_point_location            naive_pl (arr);                 // 0
  timer.stop(); 
  std::cout << "Naive construction took " << timer.time() <<std::endl;

   timer.reset(); timer.start();
  Walk_point_location             walk_pl (arr);                  // 1
  timer.stop(); 
  std::cout << "Walk construction took " << timer.time() <<std::endl;

  timer.reset(); timer.start();
  Ric_point_location              ric_pl(arr);                    // 3
  timer.stop(); 
  std::cout << "Ric construction took " << timer.time() <<std::endl;

   timer.reset(); timer.start();
  Lm_point_location               lm_pl (arr);                    // 4
  timer.stop(); 
  std::cout << "Lm (vert) construction took " << timer.time() <<std::endl;

   timer.reset(); timer.start();
  Random_lm_generator             random_g(arr);    
  Lm_random_point_location        randon_lm_pl (arr, &random_g);  // 5
  timer.stop(); 
  std::cout << "Random lm construction took " << timer.time() <<std::endl;

   timer.reset(); timer.start();
  Grid_lm_generator               grid_g(arr);      
  Lm_grid_point_location          grid_lm_pl (arr, &grid_g);      // 6
  timer.stop(); 
  std::cout << "Grid lm construction took " << timer.time() <<std::endl;

   timer.reset(); timer.start();
  Halton_lm_generator             halton_g(arr);
  Lm_halton_point_location        halton_lm_pl (arr, &halton_g);  // 7
  timer.stop(); 
  std::cout << "Halton lm construction took " << timer.time() <<std::endl;

#ifdef SEGMENTS
  timer.reset(); timer.start();
  Mide_lm_generator               mide_g(arr);
  Lm_mide_point_location          mide_lm_pl (arr, &mide_g);      // 8
  timer.stop(); 
  std::cout << "Mide lm construction took " << timer.time() <<std::endl;

  timer.reset(); timer.start();
  Trg_point_location              trg_pl(arr);                    // 2
  timer.stop(); 
  std::cout << "Trg construction took " << timer.time() <<std::endl;
#endif

  // ===> Add new point location instance here. <===

  CGAL::Object                    obj;
  Objects_vector                  objs[NUM_OF_POINT_LOCATION_STRATEGIES];
  Object_iterator                 ob_iter[NUM_OF_POINT_LOCATION_STRATEGIES];
  Arrangement_2::Vertex_const_handle    vh_ref, vh_curr;
  Arrangement_2::Halfedge_const_handle  hh_ref, hh_curr;
  Arrangement_2::Face_const_handle      fh_ref, fh_curr;
  
  Point_2               q;
  Point_iterator        piter;

  //LOCATE the points in the list using all PL strategies

  //std::cout << "Time in seconds" <<std::endl; ;
  std::cout << std::endl << "Naive" ;
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = naive_pl.locate (q);
    objs[0].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << " location took " << timer.time() <<std::endl;

  std::cout << "Walk";
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = walk_pl.locate (q);
    objs[1].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << " location took " << timer.time() <<std::endl;

  std::cout << "RIC" ;
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = ric_pl.locate (q);
    objs[2].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << " location took " << timer.time() <<std::endl;

  std::cout << "Landmarks (vertices) " <<std::endl;;
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++)
  {
    q = (*piter);
    obj = lm_pl.locate (q);
    objs[3].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << " location took " << timer.time() <<std::endl;

  std::cout << "Random LM" ;
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = randon_lm_pl.locate (q);
    objs[4].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << " location took " << timer.time() <<std::endl;

  std::cout << "Grid LM" ;
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = grid_lm_pl.locate (q);
    objs[5].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << " location took " << timer.time() <<std::endl;

  std::cout << "Halton LM" ;
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = halton_lm_pl.locate (q);
    objs[6].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << " location took " << timer.time() <<std::endl;

#ifdef SEGMENTS
  std::cout << "Middle edges LM" ;
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = mide_lm_pl.locate (q);
    objs[7].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << " location took " << timer.time() <<std::endl;

  std::cout << "Triangulation" ;
  timer.reset(); 
  timer.start(); //START
  for (piter = plist.begin(); piter != plist.end(); piter++ )
  {
    q = (*piter);
    obj = trg_pl.locate (q);
    objs[8].push_back(obj);
  }
  timer.stop(); ///END
  std::cout << " location took " << timer.time() <<std::endl;
#endif

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
    if (CGAL::assign (fh_ref, ob_iter[0][qi]))
    {    
      for (int pl_index=1; pl_index<pls_num; pl_index++)
      {
        if (! CGAL::assign(fh_curr, ob_iter[pl_index][qi]))
        {
          std::cout << "Error in point location number " 
            << pl_index << " did no return a face"<< std::endl;
          return (-1);
        }
        else if (fh_curr != fh_ref)
        {
          std::cout << "Error: point location number " 
            << pl_index << " return a different face"<< std::endl;
          return (-1);         
        }
      }  
      //if (fh_ref->is_unbounded())
      //  std::cout << "Unbounded face." << std::endl;
      //else
      //  std::cout << "Face." << std::endl;
    }

    //assign object to a halfedge
    else if (CGAL::assign (hh_ref, ob_iter[0][qi]))
    {
      //std::cout << "Halfedge: "<< hh_ref->curve() << std::endl;
      for (int pl_index=1; pl_index<pls_num; pl_index++)
      {
        if (! CGAL::assign(hh_curr, ob_iter[pl_index][qi]))
        {
          std::cout << "Error in point location number " 
            << pl_index << " did no return halfedge"<< std::endl;
          return (-1);
        }
        else if ((hh_curr != hh_ref) && (hh_curr->twin() != hh_ref))
        {
          std::cout << "Error: point location number " 
            << pl_index << " return a different halfedge"<< std::endl;
          //std::cout << "Halfedge (curr): "<< hh_curr->curve() << std::endl;
          return (-1);         
        }
      }      
    }

    //assign object to a vertex
    else if (CGAL::assign (vh_ref, ob_iter[0][qi]))
    {
      for (int pl_index=1; pl_index<pls_num; pl_index++)
      {
        if (! CGAL::assign(vh_curr, ob_iter[pl_index][qi]))
        {
          std::cout << "Error in point location number " 
            << pl_index << " did no return a vertex"<< std::endl;
          return (-1);
        }
        else if (vh_curr != vh_ref)
        {
          std::cout << "Error: point location number " 
            << pl_index << " return a different vertex"<< std::endl;
          return (-1);         
        }
      }
      //std::cout << "Vertex: "<< vh_ref->point() << std::endl;
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

#ifdef DRAW_QT
/*! */
inline Window_stream & operator<<(Window_stream & os, Arrangement_2 & arr)
{
	Edge_const_iterator ei;
	os << CGAL::BLUE;
	for (ei = arr.edges_begin(); ei != arr.edges_end(); ++ei)
		os << (*ei).curve();
	Vertex_const_iterator vi;
	os << CGAL::RED;
	for (vi = arr.vertices_begin(); vi != arr.vertices_end(); ++vi)
		os << (*vi).point();
	return os;
}

/*! */
void  draw (Arrangement_2 & arr, Points_list &plist, 
            bool draw_points, CGAL::Bbox_2  bbox,
            int argc, char * argv[])
{
  //draw pm
  QApplication app(argc, argv);
  Window_stream m_window;
  app.setMainWidget(&m_window);
  double x_min=bbox.xmin(), x_max=bbox.xmax(), 
         y_min=bbox.ymin(), y_max=bbox.ymax();
  double x,y;

  //update the bouding box to fit the points
  if (draw_points)
  {
    Point_iterator  piter;
    for (piter = plist.begin(); piter != plist.end(); piter++ )
    {
      x = CGAL::to_double(piter->x());
      y = CGAL::to_double(piter->y());
      if (x < x_min) x_min = x;
      if (x > x_max) x_max = x;
      if (y < y_min) y_min = y;
      if (y > y_max) y_max = y;
    }
  }

  m_window.resize(400,400);
  // logical window size 
  m_window.set_window(x_min-10, x_max+10, y_min-10, y_max+10);   
  m_window.show();

  m_window.lock();
  m_window << CGAL::BackgroundColor(CGAL::WHITE) << CGAL::RED;
  m_window << arr;

  //Print the points in green 
  if (draw_points)
  {
    m_window << CGAL::GREEN;
    Point_iterator  piter;
    for (piter = plist.begin(); piter != plist.end(); piter++ )
    {
      m_window << (*piter);
    }
  }

  m_window.unlock();
  app.exec();    
}
#endif

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
  bool draw_arr = false;
  bool draw_points = false;

  if ( (argc >= 4) &&
    (strcmp (argv[3], "draw") == 0))
  {
    draw_arr = true;
    if ( (argc >= 5) &&
      (strcmp (argv[4], "draw") == 0))
    {
      draw_points = true;
    }
  }

  //read curves and insert them into the arrangement 
#ifdef SEGMENTS
  Segment_reader<Traits_2> reader;
#elif defined (CONICS)
 	Conic_reader<Traits_2>  reader;
#endif

  CGAL::Bbox_2             bbox;
  Curve_list               curve_list;
  reader.read_data(curves_filename, std::back_inserter(curve_list), bbox);  

  //insert all curves into the arrangement
  CGAL::Timer timer;
  timer.reset(); 
  timer.start(); //START
  Arrangement_2  arr;
  insert (arr, curve_list.begin(), curve_list.end());
   timer.stop(); ///END
  std::cout << "Arrangement aggregate construction took " 
    << timer.time() <<std::endl;

  //read point and insert them into a list of points
  Points_list           plist;

  //read points from file into list
  if (read_points(points_filename, plist))
  {
    std::cerr << "ERROR in read_points."<<std::endl<<std::endl;
    exit(-1);
  }

  //check point location of points
  if (check_point_location(arr, plist))
  {
    std::cerr << "ERROR in check_point_location."<<std::endl<<std::endl;
    exit(-1);
  }
  
  std::cout << "Point-location test was finished O.K."<<std::endl<<std::endl;

  if (draw_arr)
  {
#ifdef DRAW_QT
    draw (arr, plist, draw_points, bbox , argc, argv);
#endif
  }

  return (0);
}
