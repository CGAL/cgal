#ifndef OVL_UTILITY_H
#define OVL_UTILITY_H

#include <CGAL/config.h> // needed for the LONGNAME flag
#include <CGAL/basic.h>
#include <vector>
//#include <CGAL/leda_real.h>
//#include <CGAL/leda_rational.h>
#include <CGAL/IO/Pm_Window_stream.h>

CGAL_BEGIN_NAMESPACE

template <class PM_, class NT_>
class Ovl_utility {
public:
  typedef   NT_     NT;
  typedef   PM_     PM;

  typedef typename PM::Traits  Traits;
  typedef typename Traits::Point_2                             Point;
  typedef typename Traits::Curve_2                             Curve;
  typedef typename Traits::X_curve                             X_curve;
  
  template <class MapOverlay>
  void  draw_and_locate_maps (const MapOverlay& ovl, 
                              const PM& pm1, 
                              const PM& pm2, 
                              CGAL::Window_stream& W)
  { 
    for (;;) {
      double  x,y;
      
      int b = W.read_mouse(x,y);
      if (b == 4) 
        break;
      
      W.clear();
      W.set_line_width(1);
      W << CGAL::BLUE;
      W << pm1;
      W << CGAL::RED;
      W << pm2;
      
      const PM& pm = ovl.subdivision();      
      
      NT z(x), w(y);
      typename Traits::Point p(z, w);
      
      typename PM::Locate_type lt;
      typename PM::Halfedge_const_handle e = pm.locate(p, lt);
      
      typename PM::Face_const_handle fh = e->face();
      if (lt == PM::UNBOUNDED_FACE) {
        std::cout << "UNBOUNDED" << endl;
      }
      else {
        W.set_node_width(4);
      
        typename PM::Face_const_handle f1 = ovl.change_notification()->get_first_face_above(fh);
        typename PM::Face_const_handle f2 = ovl.change_notification()->get_second_face_above(fh);
        
        //if (f1 == fh || f2 == fh)
        //  std::cout<<"NULL faces pointers\n";
        
        leda_color fg_col;
        if (f1 != fh && !(f1->is_unbounded()) && f2 != fh && !(f2->is_unbounded()) )
          fg_col = leda_violet;
        else if (f1 != fh && !(f1->is_unbounded()))
          fg_col = leda_blue;
        else if (f2 != fh && !(f2->is_unbounded()))
          fg_col = leda_red;
        else
          fg_col = leda_orange;
        
        W.set_color(fg_col);
        //W.set_fill_color(fg_col);

        typename PM::Ccb_halfedge_const_circulator cc = fh->outer_ccb();
        
        W.set_line_width(3);
        do {
          W << (*cc).curve();
          W << cc->source()->point();
        } while (++cc != fh->outer_ccb());
        
        for (typename PM::Holes_const_iterator hit = fh->holes_begin(); hit != fh->holes_end(); ++hit) {
          typename PM::Ccb_halfedge_const_circulator cc(*hit);
          do{
            W << cc->curve();
          } while (++cc != *hit);
        }
        
        for (typename PM::Vertex_const_iterator  v_iter = pm.vertices_begin(); 
             v_iter != pm.vertices_end(); ++v_iter){
          if (v_iter->get_first_vertex_above() != NULL && v_iter->get_second_vertex_above() != NULL)
            W.set_color(leda_violet);
          else if (v_iter->get_first_vertex_above() != NULL && v_iter->get_second_vertex_above() == NULL)
            W<<CGAL::BLUE;
          else if (v_iter->get_second_vertex_above() != NULL && v_iter->get_first_vertex_above() == NULL)
            W << CGAL::RED;
          else
            W << CGAL::ORANGE;
          
          W << v_iter->point();
        }
        
        if (lt == PM::EDGE){
          if (e.operator->()->get_first_halfedge_above() != NULL)
            std::cout<<"edge belongs to first subdivision\n";
          if (e.operator->()->get_second_halfedge_above() != NULL)
            std::cout<<"edge belongs to second subdivision\n";
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


