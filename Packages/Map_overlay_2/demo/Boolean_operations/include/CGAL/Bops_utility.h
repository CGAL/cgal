#ifndef BOPS_UTILITY_H
#define BOPS_UTILITY_H

#include <CGAL/config.h> // needed for the LONGNAME flag
#include <CGAL/basic.h>
#include <vector>
//#include <CGAL/leda_real.h>
//#include <CGAL/leda_rational.h>
#include <CGAL/IO/Pm_Window_stream.h>

CGAL_BEGIN_NAMESPACE

template <class PM_, class NT_>
class Bops_utility {
public:
  typedef   NT_     NT;
  typedef   PM_     PM;

  typedef typename PM::Traits  Traits;
  typedef typename Traits::Point_2                             Point;
  typedef typename Traits::Curve_2                             Curve;
  typedef typename Traits::X_curve                             X_curve;

  template <class Bops>
  void  draw_and_locate_maps (Bops& bops , 
                              const PM& pm1, 
                              const PM& pm2, 
                              CGAL::Window_stream& W)
  {    
    for (;;) {
      double  x,y;
      
      int b = W.read_mouse(x,y);
      if (b==7) 
        break;
      
      W.set_node_width(3);
      W.set_line_width(1);
      W.clear();
      W<<CGAL::BLUE;
      W<<pm1;
      W<<CGAL::RED;
      W<<pm2;
      
      leda_color fg_col = leda_violet;
      W.set_color(fg_col);
      //W.set_fill_color(fg_col);
      
      if (b == 4 || b == 5 || b == 6){
        typename Bops::Faces_container      faces_result;
        typename Bops::Halfedges_container  halfedges_result;
        typename Bops::Vertices_container   vertices_result;
        
        if (b==4)
          bops.intersection(faces_result,halfedges_result,vertices_result);
        
        else if (b==5)       
          bops.Union(faces_result,halfedges_result,vertices_result);
        
        else if (b==6)
          bops.symmetric_difference(faces_result, halfedges_result, vertices_result);
        
        std::cout<<"number of faces is "<< faces_result.size()<<std::endl;
        std::cout<<"The number of resulting halfedges is "<<halfedges_result.size()<<std::endl;
        std::cout<<"The number of resulting vertices is "<<vertices_result.size()<<std::endl;
        
        W.set_line_width(3);
        for (list<typename PM::Face_const_handle>::iterator  f_iter = faces_result.begin(); 
             f_iter != faces_result.end(); ++f_iter){
          typename PM::Face_const_handle fh = *f_iter;
          if (fh->is_unbounded())
            continue;
          
          W.set_color(fg_col);
          //W.set_fill_color(fg_col);
          
          typename PM::Ccb_halfedge_const_circulator cc = fh->outer_ccb();
          /*leda_list<leda_point> points_list;
            do {
            leda_point p = cc->source()->point().to_point();
            
            points_list.push_back(p);
            } while (++cc != fh->outer_ccb());
            W.draw_filled_polygon(points_list);*/
          
          for (typename PM::Holes_const_iterator hit = fh->holes_begin(); hit != fh->holes_end(); ++hit) {
            leda_list<leda_point> holes_points_list;
            typename PM::Ccb_halfedge_const_circulator cc(*hit);
            W.set_fill_color(leda_black);
            W << CGAL::WHITE;
            do{
              W << cc->curve();
              
              //leda_point p = cc->source()->point().to_point();
              
              //points_list.push_back(p);
            } while (++cc != *hit);
            //W.draw_filled_polygon(holes_points_list);
          }
          
          do {
            //W <<CGAL::WHITE;
            W << (*cc).curve();
          } while (++cc != fh->outer_ccb());
        }
        
        for (list<typename PM::Halfedge_const_handle>::iterator h_iter = halfedges_result.begin(); 
             h_iter != halfedges_result.end(); ++h_iter, ++h_iter){
          typename PM::Halfedge_const_handle h = *h_iter;
          
          W.set_color(fg_col);
          W.set_node_width(4);
          W << h->curve();
        }
        
        for (list<typename PM::Vertex_const_handle>::iterator  v_iter = vertices_result.begin(); 
             v_iter != vertices_result.end(); ++v_iter){
          W.set_node_width(4);
          
          typename PM::Vertex_const_handle vh = *v_iter;
          //Vertex_const_iterator vh = v_iter;
          
          /*if (tmp_notf.get_first_halfedge_above(vh) != vh->incident_halfedges() && tmp_notf.get_second_halfedge_above(vh) != vh->incident_halfedges())
            W.set_color(leda_violet);
            else if (tmp_notf.get_first_halfedge_above(vh) != vh->incident_halfedges())
            W<<CGAL::BLUE;
            else if (tmp_notf.get_second_halfedge_above(vh) != vh->incident_halfedges())
          W << CGAL::RED;
          else
          W << CGAL::ORANGE;*/
          
          W.set_color(leda_violet);
          W << vh->point();
        }
      }
    }
  }
  
  void  scan_segment_planar_map(const char* filename, PM& pm)
  {
    std::vector<Curve> curves;
    
    std::ifstream file(filename);
    int num_curves;
    file >> num_curves;
    while (num_curves--) {
    //NT x,y;
      double x,y; //(actually int)
      file >> x >> y;
      NT nx(x),ny(y);
      Point s(nx,ny);
      file >> x >> y;
      NT mx(x),my(y);
      Point t(mx,my);
      
      curves.push_back(Curve(s,t));
    }
    
    file.close();
    
    Traits traits;
    CGAL::sweep_to_construct_planar_map_2(curves.begin(), 
                                          curves.end(),
                                          traits, pm);
  }
  
  void  scan_segment_pmwx(const char* filename, PM& pmwx)
  { 
    std::ifstream file(filename);
    int num_curves;
    file >> num_curves;
    while (num_curves--) {
      //NT x,y;
      double x,y; //(actually int)
      file >> x >> y;
      NT nx(x),ny(y);
      Point s(nx,ny);
      file >> x >> y;
      NT mx(x),my(y);
      Point t(mx,my);
      
      pmwx.insert(Curve(s,t));
    }
    
    file.close();
  }

  void scan_polyline_planar_map(const char * filename, PM & pm)
  {  
    std::vector<Curve> curves;
    std::ifstream file(filename);
    
    int  num_polylines;
    double  x,y; // only used to read file not used in computations
    X_curve  polyline;
    
    file >> num_polylines;
    while (num_polylines--) {
      int  num_x_curves;
      file >> num_x_curves;
      while (num_x_curves--) {
        file >> x >> y;
        Point s(x,y);
        polyline.push_back(s);
      }
      
      curves.push_back(polyline);
      
      polyline.clear();
    }
    
    file.close();

    Traits traits;
    CGAL::sweep_to_construct_planar_map_2(curves.begin(), 
                                          curves.end(),
                                          traits, pm);
  }
  
  void scan_polyline_pmwx(const char * filename, PM & pmwx)
  {  
    std::ifstream file(filename);
    
    int  num_polylines;
    double  x,y; // only used to read file not used in computations
    X_curve  polyline;
    
    file >> num_polylines;
    while (num_polylines--) {
      int  num_x_curves;
      file >> num_x_curves;
      while (num_x_curves--) {
        file >> x >> y;
        Point s(x,y);
        polyline.push_back(s);
      }
      
      pmwx.insert(polyline);
      
    polyline.clear();
    
    }
    file.close();
  }
  
  void  scan_seg_circ_planar_map(const char* filename, PM& pm)
  {
    std::vector<Curve> curves;
    
    std::ifstream file(filename);
    int num_curves;
    file >> num_curves;
    char   type;
    int    i_arc;
    
    for (i_arc = 0; i_arc < num_curves; ++i_arc)
      {
        // Read the arc type.
        file >> type;
        
        std::cout << "Scanning arc no. " << i_arc + 1;
        
        // A full circle (c) or a circular arc (a):
        if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
          {
            // Read the circle, using the format "x0 y0 r^2"
            NT  x0, y0, r2;
            
            file >> x0 >> y0 >> r2;
            
            typename Traits::Circle      circle (Point (x0, y0), r2, CGAL::CLOCKWISE);
            
            if (type == 'c' || type == 'C')
              {
                std::cout << " (full circle)." << std::endl;
                curves.push_back (Curve(circle));
              }
            else
              {
                std::cout << " (circular arc)." << std::endl;
                
                // Read the end points.
                NT  x1, y1, x2, y2;
                
                file >> x1 >> y1 >> x2 >> y2;
                
                Point      source (x1, y1);
                Point      target (x2, y2);
                
                curves.push_back (Curve (circle, source, target));
              }
          }
        // A segment (s):
        else if (type == 's' || type == 'S')
          {
            std::cout << " (segment)." << std::endl;
            // Read the end points.
            NT  x1, y1, x2, y2;
            
            file >> x1 >> y1 >> x2 >> y2;
            
            Point      source (x1, y1);
            Point      target (x2, y2);
            
            curves.push_back (Curve (typename Traits::Segment (source, target)));
          }
        else
          {
            std::cout << "Unknown arc type '" << type 
                      << "' - Stopping here." << std::endl;
            exit(1);
          }
      }  
    
    file.close();
    
    Traits traits;
    
    std::cout<<"Sweeping all input curves"<<std::endl;
    CGAL::sweep_to_construct_planar_map_2(curves.begin(), 
                                          curves.end(),
                                          traits, pm);
  }
  
  void  scan_seg_circ_pmwx(const char* filename, PM& pmwx)
  {
    std::ifstream file(filename);
    int num_curves;
    file >> num_curves;
    char   type;
    int    i_arc;
  
    for (i_arc = 0; i_arc < num_curves; ++i_arc)
      {
        // Read the arc type.
        file >> type;
        
        std::cout << "Inserting arc no. " << i_arc + 1;
        
        // A full circle (c) or a circular arc (a):
        if (type == 'c' || type == 'C' || type == 'a' || type == 'A')
          {
            // Read the circle, using the format "x0 y0 r^2"
            NT  x0, y0, r2;
            
            file >> x0 >> y0 >> r2;
            
            typename Traits::Circle    circle (Point (x0, y0), r2, CGAL::CLOCKWISE);

            if (type == 'c' || type == 'C')
              {
                std::cout << " (full circle)." << std::endl;
                pmwx.insert (Curve(circle));
              }
            else
              {
                std::cout << " (circular arc)." << std::endl;

                // Read the end points.
                NT  x1, y1, x2, y2;
                
                file >> x1 >> y1 >> x2 >> y2;
                
                Point      source (x1, y1);
                Point      target (x2, y2);
                
                pmwx.insert (Curve (circle, source, target));
              }
          }
        // A segment (s):
        else if (type == 's' || type == 'S')
          {
            std::cout << " (segment)." << std::endl;
            // Read the end points.
            NT  x1, y1, x2, y2;
            
            file >> x1 >> y1 >> x2 >> y2;

            Point      source (x1, y1);
            Point      target (x2, y2);
            
            pmwx.insert (Curve (typename Traits::Segment (source, target)));
          }
        else
          {
            std::cout << "Unknown arc type '" << type 
                      << "' - Stopping here." << std::endl;
            exit(1);
          }
      }
    
    file.close();
  }
  
#ifdef CGAL_LEDA_RATIONAL_H 
  void calc_window_size(const PM &pm, double &min_x, double &max_x, double &min_y)
  {
    typename PM::Vertex_const_iterator v_iter = pm.vertices_begin();
    
    max_x = min_x = CGAL::to_double(v_iter->point().xcoordD());
    min_y = CGAL::to_double(v_iter->point().ycoordD());
    
    for (v_iter = pm.vertices_begin(); 
         v_iter != pm.vertices_end(); ++v_iter){
      NT x = v_iter->point().xcoord(), y = v_iter->point().ycoord();
      
      double  dx, dy;
      dx = CGAL::to_double(x);
      dy = CGAL::to_double(y);
      if (dx > max_x) max_x = dx;
      if (dx < min_x) min_x = dx;
      if (dy < min_y) min_y = dy;
    }
  }

#else
  void calc_window_size(const PM &pm, double &min_x, double &max_x, double &min_y)
  {
    typename PM::Vertex_const_iterator v_iter = pm.vertices_begin();
    
    max_x = min_x = CGAL::to_double(v_iter->point().x());
    min_y = CGAL::to_double(v_iter->point().y());
    
    for (v_iter = pm.vertices_begin(); 
         v_iter != pm.vertices_end(); ++v_iter){
      NT x = v_iter->point().x(), y = v_iter->point().y();
      
      double  dx, dy;
      dx = CGAL::to_double(x);
      dy = CGAL::to_double(y);
      if (dx > max_x) max_x = dx;
      if (dx < min_x) min_x = dx;
      if (dy < min_y) min_y = dy;
    }
  }
#endif
};


CGAL_END_NAMESPACE

#endif


