
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Timer.h>
#include <CGAL/IO/Dxf_reader.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>


typedef CGAL::Simple_cartesian<double>            K;
typedef std::list<std::pair<K::Point_2, double> > Dxf_polygon;
typedef std::list<Dxf_polygon>                    Dxf_polygons_list;

typedef std::pair<K::Point_2, K::FT>              Dxf_circle;                      
typedef std::list<Dxf_circle>                     Dxf_circles_list;
typedef CGAL::Dxf_reader<K>                       Dxf_reader;                  

typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > 
                                                      NT;
typedef CGAL::Simple_cartesian<NT>                    Kernel;
typedef Kernel::FT                                    FT;
typedef Kernel::Point_2                               Point_2;
typedef Kernel::Circle_2                              Circle_2;

typedef CGAL::Gps_circle_segment_traits_2<Kernel>     Traits_2;
typedef Traits_2::Point_2                             Arc_point_2;
typedef Traits_2::Curve_2                             Curve_2;
typedef Traits_2::X_monotone_curve_2                  X_monotone_curve_2;
typedef Traits_2::Polygon_2                           Circ_polygon;
typedef Traits_2::Polygon_with_holes_2                Circ_polygon_with_holes;
typedef std::vector<Circ_polygon>                     Circ_pgn_vec;
typedef std::vector<Circ_polygon_with_holes>          Circ_pgn_with_holes_vec;
typedef CGAL::General_polygon_set_2<Traits_2>         Gps;


int main(int argc, char **argv)
{
  if(argc < 2)
  {
    std::cerr<<"Missing DXF file"<<std::endl;
    exit(-1);
  }

  std::ifstream input_file (argv[1]);
  if(!input_file.is_open())
  {
    std::cerr<<"Failed to open the file"<<std::endl;
    exit(-1);
  }

  Traits_2 tr;
  Gps gps(tr);

  Dxf_polygons_list polygons;
  Dxf_circles_list  circles;
  Dxf_reader        reader;

  reader(input_file, polygons, circles);

  std::cout << "Read:  " << polygons.size() << " polygons, and " << circles.size() << " circles" << std::endl;


  Circ_pgn_vec circ_polygons;
  circ_polygons.reserve(circles.size());
  
  Circ_pgn_with_holes_vec circ_pgn_with_holes;
  circ_pgn_with_holes.reserve(polygons.size());


  for(Dxf_circles_list::iterator circ_iterator = circles.begin();
      circ_iterator != circles.end();
      ++circ_iterator)
  {
    const Dxf_circle& dxf_circle = *circ_iterator;
    
    FT center_x = dxf_circle.first.x();
    FT center_y = dxf_circle.first.y();
    FT rad = dxf_circle.second;

    Curve_2 circ(Point_2(center_x, center_y), rad);
    std::vector<CGAL::Object> obj_vec;
    
    obj_vec.reserve(2);
    tr.make_x_monotone_2_object()(circ, std::back_inserter(obj_vec));
    
    CGAL_assertion(obj_vec.size() == 2);
    X_monotone_curve_2 cv1, cv2;
    CGAL::assign(cv1, obj_vec[0]);
    CGAL::assign(cv2, obj_vec[1]);
    Circ_polygon pgn;
    pgn.push_back(cv1);
    pgn.push_back(cv2);
    circ_polygons.push_back(pgn);
  }

  circles.clear();

  for(Dxf_polygons_list::iterator it = polygons.begin(); it != polygons.end(); ++it)
  {
    Circ_polygon curr_pgn;
    for(Dxf_polygon::iterator pit = it->begin(); pit != it->end(); ++pit)
    {
      Dxf_polygon::iterator next = pit;
      ++next;

      Point_2 ps((pit->first.x()),
                 (pit->first.y()));
      Point_2 pt;
      if(next == it->end())
      {
        pt = Point_2((it->begin()->first.x()),
                     (it->begin()->first.y()));
      }
      else
      {
        pt = Point_2((next->first.x()),
                     (next->first.y()));
      }
      
      if(pit->second) 
      {
        const FT bulge = (pit->second);
        const FT common = (1 - CGAL::square(bulge)) / (4*bulge);
        const FT x_coord = ((ps.x() + pt.x())/2) + common*(ps.y() - pt.y());
        const FT y_coord = ((ps.y() + pt.y())/2) + common*(pt.x() - ps.x());
        const FT sqr_bulge = CGAL::square(bulge);
        const FT sqr_rad = CGAL::squared_distance(ps, pt) * (1/sqr_bulge + 2 + sqr_bulge) / 16; 

        CGAL_assertion(ps != pt);
       
        Circle_2 supp_circ;
        if(pit->second > 0)
          supp_circ = Circle_2(Point_2(x_coord, y_coord), sqr_rad);
        else
          supp_circ = Circle_2(Point_2(x_coord, y_coord), sqr_rad, CGAL::CLOCKWISE);

        Curve_2 circ_arc(supp_circ, 
                         Arc_point_2(ps.x(), ps.y()),
                         Arc_point_2(pt.x(), pt.y()));
        std::vector<CGAL::Object> obj_vec;
        tr.make_x_monotone_2_object()(circ_arc, std::back_inserter(obj_vec));
        for(unsigned int i=0; i<obj_vec.size(); ++i)
        {
          X_monotone_curve_2 cv;
          if(CGAL::assign(cv, obj_vec[i]))
            curr_pgn.push_back(cv);
        }
      }
      else
      {
        if( ps == pt)
          continue;

        curr_pgn.push_back(X_monotone_curve_2(ps, pt));
      }
    }

    if(curr_pgn.orientation() == CGAL::CLOCKWISE)
      curr_pgn.reverse_orientation();

    Circ_polygon_with_holes pgn_with_holes;
    gps.simplify(curr_pgn, pgn_with_holes);
    circ_pgn_with_holes.push_back(pgn_with_holes);    
  }
  polygons.clear();

  CGAL::Timer t;
  
  std::cout<<"Performing union\n";
  
  t.start();
  gps.join(circ_polygons.begin(), circ_polygons.end(),
           circ_pgn_with_holes.begin(), circ_pgn_with_holes.end());
  t.stop();
  
  std::cout<<"union time : "<< t.time()<<" seconds\n";

  std::cout<<"|V| = " << gps.arrangement().number_of_vertices()<<"\n";
  std::cout<<"|E| = " << gps.arrangement().number_of_edges()<<"\n";
  std::cout<<"|F| = " << gps.arrangement().number_of_faces()<<"\n";

  input_file.close();  
  return 0;
}
