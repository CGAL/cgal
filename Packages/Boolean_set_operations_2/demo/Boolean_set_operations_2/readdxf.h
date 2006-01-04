// Descriptions of the file format can be found at
// http://www.autodesk.com/techpubs/autocad/acad2000/dxf/
// http://www.tnt.uni-hannover.de/soft/compgraph/fileformats/docs/DXF.ascii

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Dxf_reader.h>
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>

typedef CGAL::Simple_cartesian<double> SC;
typedef SC::Point_2 double_Point_2;
typedef SC::Circle_2 double_Circle_2;

typedef std::list<std::pair<double_Point_2, double> > Dxf_Polygon;
typedef std::list<Dxf_Polygon>                            Polygons;
typedef std::list<double_Circle_2> double_Circles;
typedef std::vector<Polygon>                            Polygons_vec;
typedef CGAL::_One_root_point_2<Coord_type, true>       One_root_point_2;

extern Polygon_set                       red_set;
int readdxf(std::istream& input_file, Polygon_set* pgn_set,CGAL::Qt_widget* w )
{
  Polygons polygons;
  double_Circles circles;
  CGAL::Dxf_reader<SC> reader;

  reader(input_file, polygons, circles);

  std::cout << "Read:  " << polygons.size() << " polygons, and " << circles.size() << " circles" << std::endl;

  Traits tr;

  Polygons_vec circ_polygons;
  for(double_Circles::iterator circ_iterator = circles.begin();
      circ_iterator != circles.end();
      ++circ_iterator)
  {
    const double_Circle_2& dbl_circ = *circ_iterator;
    Coord_type sqr_rad = dbl_circ.squared_radius();
    Coord_type center_x = dbl_circ.center().x();
    Coord_type center_y = dbl_circ.center().y();

    
    Curve circ(Point(center_x, center_y), sqr_rad);
    std::vector<CGAL::Object> obj_vec;
    obj_vec.reserve(2);

    tr.make_x_monotone_2_object()(circ, std::back_inserter(obj_vec));
    CGAL_assertion(obj_vec.size() == 2);

    XCurve cv1, cv2;
    CGAL::assign(cv1, obj_vec[0]);
    CGAL::assign(cv2, obj_vec[1]);
    Polygon pgn;
    pgn.push_back(cv1);
    pgn.push_back(cv2);
    circ_polygons.push_back(pgn);
  }

  for(Polygons::iterator it = polygons.begin(); it != polygons.end(); it++)
  {
    Polygon curr_pgn;
    for(Dxf_Polygon::iterator pit = it->begin(); pit != it->end(); pit++)
    {
      
      Dxf_Polygon::iterator next = pit;
      ++next;

      Point ps(pit->first.x(), pit->first.y());;
      Point pt;
      if(next == it->end())
      {
        pt = Point(it->begin()->first.x(), it->begin()->first.y());
      }
      else
      {
        pt = Point(next->first.x(), next->first.y());
      }
      
      if(pit->second) 
      {
        Coord_type bulge(pit->second);
        Coord_type common((1 - CGAL::square(bulge))/(4*(bulge)));
        Coord_type x_coord((ps.x() + pt.x())/2 + common * (ps.x() - pt.x()));
        Coord_type y_coord((ps.y() + pt.y())/2 + common * (pt.y() - ps.y()));
        Coord_type sqr_rad(CGAL::squared_distance(ps, pt)*CGAL::square(1/(bulge)+ (bulge)) / 16);

        CGAL_assertion(ps != pt);
       
        Circle supp_circ(Point(x_coord, y_coord), sqr_rad);
        Curve circ_arc(supp_circ, 
                         One_root_point_2(ps.x(), ps.y()),
                         One_root_point_2(pt.x(), pt.y()));
        std::vector<CGAL::Object> obj_vec;
        tr.make_x_monotone_2_object()(circ_arc, std::back_inserter(obj_vec));
        for(unsigned int i=0; i<obj_vec.size(); ++i)
        {
          XCurve cv;
          if(CGAL::assign(cv, obj_vec[i]))
            curr_pgn.push_back(cv);
        }
      }
      else
      {
        if( ps == pt)
          continue;

        //curr_pgn.push_back(XCurve(ps, pt));
        XCurve cv;
        Curve curve_seg(ps, pt);
        std::vector<CGAL::Object> obj_vec;
        tr.make_x_monotone_2_object()(curve_seg, std::back_inserter(obj_vec));
        CGAL_assertion(obj_vec.size() == 1);
        if(CGAL::assign(cv, obj_vec[0]))
          curr_pgn.push_back(cv);
      }
    }
    if(!tr.is_valid_2_object()(curr_pgn))
    {
      std::cout<<"invalid polygon!!!\n";
      std::cout<<curr_pgn<<"\n";
   red_set.join(curr_pgn);
      return(0);
    }
    circ_polygons.push_back(curr_pgn);
  }

  pgn_set->clear();
  pgn_set->join(circ_polygons.begin(), circ_polygons.end());
  return 0;
}
