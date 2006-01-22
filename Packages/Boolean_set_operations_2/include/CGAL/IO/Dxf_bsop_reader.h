
#ifndef DXF_BSOP_READER_H
#define DXF_BSOP_READER_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Dxf_reader.h>
#include <iostream>

#include <list>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gmpq.h>

CGAL_BEGIN_NAMESPACE

template <class Kernel_>
class Dxf_bsop_reader
{
  typedef CGAL::Simple_cartesian<double>            K;
  typedef std::list<std::pair<K::Point_2, double> > Dxf_polygon;
  typedef std::list<Dxf_polygon>                    Dxf_polygons_list;

  typedef std::pair<K::Point_2, K::FT>              Dxf_circle;                      
  typedef std::list<Dxf_circle>                     Dxf_circles_list;
  typedef CGAL::Dxf_reader<K>                       Dxf_reader;        

  typedef Kernel_                                       Kernel;
  typedef typename Kernel::FT    FT;

  typedef typename Kernel::Point_2                       Point_2;
  typedef typename Kernel::Circle_2                      Circle_2;


  typedef CGAL::Gps_circle_segment_traits_2<Kernel>     Traits_2;
  typedef typename Traits_2::Point_2                    Arc_point_2;
  typedef typename Traits_2::Curve_2                    Curve_2;
  typedef typename Traits_2::X_monotone_curve_2         X_monotone_curve_2;
  typedef typename Traits_2::Polygon_2                  Circ_polygon;
  typedef typename Traits_2::Polygon_with_holes_2       Circ_polygon_with_holes;

  typedef std::vector<Circ_polygon>                     Circ_pgn_vec;
  typedef std::vector<Circ_polygon_with_holes>          Circ_pgn_with_holes_vec;

  typedef CGAL::General_polygon_set_2<Traits_2>         Gps;

  FT convert_double(const double& x)
  {
    int   denom = 1000000;
    int   numer = static_cast<int> (x * denom + CGAL  ::sign(x)*0.5);

    while (numer != 0 && numer % 10 == 0)
    {
      numer /= 10;
      denom /= 10;
    }

    return (CGAL::Gmpq (numer, denom));
  }

public:
  template <class Out1, class Out2>
  void operator() (std::ifstream& input,
                   Out1 pgns,
                   Out2 pgns_with_holes,
                   bool simplify = true)
  {
  
    Gps gps;
    Traits_2 tr;

    Dxf_polygons_list polygons;
    Dxf_circles_list  circles;
    Dxf_reader        reader;

    reader(input, polygons, circles);

    for(Dxf_circles_list::iterator circ_iterator = circles.begin();
        circ_iterator != circles.end();
        ++circ_iterator)
    {
      const Dxf_circle& dxf_circle = *circ_iterator;
      
      FT center_x = convert_double(dxf_circle.first.x());
      FT center_y = convert_double(dxf_circle.first.y());
      FT rad = convert_double(dxf_circle.second);

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
      *pgns++ = pgn;
    }

    circles.clear();

    for(Dxf_polygons_list::iterator it = polygons.begin(); it != polygons.end(); ++it)
    {
      Circ_polygon curr_pgn;
      for(Dxf_polygon::iterator pit = it->begin(); pit != it->end(); ++pit)
      {
        Dxf_polygon::iterator next = pit;
        ++next;

        Point_2 ps((convert_double(pit->first.x())),
                  (convert_double(pit->first.y())));
        Point_2 pt;
        if(next == it->end())
        {
          pt = Point_2(convert_double(it->begin()->first.x()),
                      convert_double(it->begin()->first.y()));
        }
        else
        {
          pt = Point_2(convert_double(next->first.x()),
                      convert_double(next->first.y()));
        }
        
        if(pit->second) 
        {
          const FT bulge = convert_double(pit->second);
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

      if(!simplify)
      {
        Circ_polygon_with_holes pgn_with_holes(curr_pgn);
        *pgns_with_holes++ = pgn_with_holes;
      }
      else
      {
        Circ_polygon_with_holes pgn_with_holes;
        gps.simplify(curr_pgn, pgn_with_holes);
        *pgns_with_holes++ = pgn_with_holes;
      }
    }
    polygons.clear();
    }
};

CGAL_END_NAMESPACE
#endif
