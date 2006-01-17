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
typedef std::vector<Polygon_with_holes>                 Polygons_vec;
typedef CGAL::_One_root_point_2<Coord_type, true>       One_root_point_2;

extern Polygon_set                       red_set;
int readdxf(std::istream& input_file,
            Polygon_set* pgn_set,
            CGAL::Qt_widget* w, 
            CGAL::Bbox_2& box,
            bool are_simple = true)
{
  bool first_curve = true;
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

    
    Circle circle = Circle(Point_2(center_x, center_y), sqr_rad);
    Curve circ(circle);
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
    circ_polygons.push_back(Polygon_with_holes(pgn));
    if(first_curve)
    {
      first_curve = false;
      box = dbl_circ.bbox();
    }
    else
      box = box + dbl_circ.bbox();

  }

  for(Polygons::iterator it = polygons.begin(); it != polygons.end(); it++)
  {
    Polygon curr_pgn;
    for(Dxf_Polygon::iterator pit = it->begin(); pit != it->end(); pit++)
    {
      
      Dxf_Polygon::iterator next = pit;
      ++next;

      Point_2 ps(pit->first.x(), pit->first.y());;
      Point_2 pt;
      if(next == it->end())
      {
        pt = Point_2(it->begin()->first.x(), it->begin()->first.y());
      }
      else
      {
        pt = Point_2(next->first.x(), next->first.y());
      }
      
      if(pit->second) 
      {
        const Coord_type bulge = (pit->second);
        const Coord_type common = (1 - CGAL::square(bulge)) / (4*bulge);
        const Coord_type x_coord = ((ps.x() + pt.x())/2) + common*(ps.y() - pt.y());
        const Coord_type y_coord = ((ps.y() + pt.y())/2) + common*(pt.x() - ps.x());
        const Coord_type sqr_bulge = CGAL::square(bulge);
        const Coord_type sqr_rad = CGAL::squared_distance(ps, pt) * (1/sqr_bulge + 2 + sqr_bulge) / 16; 

        CGAL_assertion(ps != pt);
       
        Circle supp_circ;
        if(pit->second > 0)
          supp_circ = Circle(Point_2(x_coord, y_coord), sqr_rad);
        else
          supp_circ = Circle(Point_2(x_coord, y_coord), sqr_rad, CGAL::CLOCKWISE);
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
          if(first_curve)
          {
            first_curve = false;
            box = cv.bbox();
          }
          else
            box = box + cv.bbox();
        }
      }
      else
      {
        if( ps == pt)
          continue;

        XCurve cv(ps, pt);
        curr_pgn.push_back(cv);
        if(first_curve)
        {
          first_curve = false;
          box = cv.bbox();
        }
        else
          box = box + cv.bbox();
      }
    }
    if(curr_pgn.orientation() == CGAL::CLOCKWISE)
      curr_pgn.reverse_orientation();
   /* if(!tr.is_valid_2_object()(curr_pgn))
    {
      std::cout<<"invalid polygon!!!\n";
      return(0);
    }*/
    if(!are_simple)
    {
      Polygon_with_holes pgn_with_holes;
      red_set.simplify(curr_pgn, pgn_with_holes);
      circ_polygons.push_back(pgn_with_holes);
    }
    else
    {
      circ_polygons.push_back(Polygon_with_holes(curr_pgn));
    } 
  }

  pgn_set->clear();
  pgn_set->join(circ_polygons.begin(), circ_polygons.end());
  return 0;
}
