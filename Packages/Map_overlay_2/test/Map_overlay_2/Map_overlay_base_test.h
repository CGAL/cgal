#ifndef  MAP_OVERLAY_BASE_TEST
#define  MAP_OVERLAY_BASE_TEST

#include <CGAL/Cartesian.h>

#ifndef CGAL_PLANAR_MAP_2
#include <CGAL/Planar_map_2.h>
#endif

#ifndef CGAL_ARR_2_OVERLAY_DCEL_H
#include <CGAL/Map_overlay_default_dcel.h>
#endif

#ifndef CGAL_MAP_OVERLAY_DEFAULTNOTIFIER_H
#include <CGAL/Map_overlay_default_notifier.h>
#endif

#ifndef CGAL_MAP_OVERLAY_H
#include <CGAL/Map_overlay.h>
#endif
 
//#include <CGAL/Arrangement_2.h>

#include <CGAL/sweep_to_construct_planar_map_2.h>

// Quotient is included anyway, because it is used to read
// data files. Quotient can read both integers and fractions.
// leda rational will only read fractions.
//#include <CGAL/Quotient.h> 
#include <iostream>
#include <list>
#include <string>
#include <fstream>

// we use the namespace std for compatability with MSVC

template <class Subdivision, class Read_curve>
class Map_overlay_base_test
{
protected:
  typedef typename Subdivision::Traits                  Traits;
  
  typedef typename Traits::Point                                 Point;
  typedef typename Traits::X_curve                               X_curve;
  typedef typename Traits::Curve                                 Curve;
  
  //typedef CGAL::Map_overlay_default_dcel<Traits>                Dcel;
  //typedef CGAL::Planar_map_2<Dcel, Traits>                      PM;
  
  typedef typename Subdivision::Vertex                           Vertex;
  typedef typename Subdivision::Face                             Face;
  typedef typename Subdivision::Halfedge                         Halfedge;
  typedef typename Subdivision::Vertex_handle                    Vertex_handle;
  typedef typename Subdivision::Halfedge_handle                  Halfedge_handle;
  typedef typename Subdivision::Face_handle                      Face_handle;
  typedef typename Subdivision::Vertex_const_handle              Vertex_const_handle;
  typedef typename Subdivision::Halfedge_const_handle            Halfedge_const_handle;
  typedef typename Subdivision::Face_const_handle                Face_const_handle;
  typedef typename Subdivision::Vertex_iterator                  Vertex_iterator;
  typedef typename Subdivision::Vertex_const_iterator            Vertex_const_iterator;
  typedef typename Subdivision::Halfedge_iterator                Halfedge_iterator;
  typedef typename Subdivision::Halfedge_const_iterator          Halfedge_const_iterator;
  typedef typename Subdivision::Face_iterator                    Face_iterator;
  typedef typename Subdivision::Face_const_iterator              Face_const_iterator;
  typedef typename Subdivision::Ccb_halfedge_circulator          Ccb_halfedge_circulator;
  typedef typename Subdivision::Ccb_halfedge_const_circulator  
                                                          Ccb_halfedge_const_circulator;
  typedef typename Subdivision::Halfedge_around_vertex_const_circulator  
                                                    Halfedge_around_vertex_const_circulator;
 
  typedef typename Subdivision::Holes_iterator                    Holes_iterator;
  typedef typename Subdivision::Holes_const_iterator              Holes_const_iterator;
  typedef typename Subdivision::Locate_type                       Locate_type;
  typedef typename Subdivision::Traits_wrap                       Traits_wrap;
  
  typedef typename Subdivision::Planar_map                        Planar_map;
  
  //typedef CGAL::Arr_base_node<Curve>                     Base_node;
  //typedef CGAL::Map_overlay_default_dcel<Traits, 
  //  CGAL::Arr_2_vertex_base<Point>, 
  //  CGAL::Arr_2_halfedge_base<Base_node>, 
  //  CGAL::Arr_2_face_base>   Arr_Dcel;
  //typedef CGAL::Arrangement_2<Arr_Dcel,Traits,Base_node>        Arrangement;
  //typedef Arrangement::Change_notification                      PMWXChangeNotification; 
  
  //typedef Arrangement::Halfedge_iterator                        Arr_halfedge_iterator;
  
  typedef CGAL::Pm_naive_point_location<Planar_map>                       PmNaivePL;
  typedef CGAL::Pm_walk_along_line_point_location<Planar_map>             PmWalkPL;
  
  typedef CGAL::Map_overlay_default_notifier<Planar_map> 
                                                   MapOverlay_change_notification;
  typedef CGAL::Map_overlay_2<Subdivision,MapOverlay_change_notification>  MapOverlay; 
  
  //MapOverlay  map_overlay;  
  std::map<const void*, Vertex_const_handle>    vertices;
  std::map<const void*, Halfedge_const_handle>  halfedges;
  std::map<const void*, Face_const_handle>      faces;
  
public:
  Map_overlay_base_test() {}

protected:
  /****************************
   * Class Implementation
   ****************************/
 
  int get_next_int(std::ifstream& file)
  {
    int   num = 0;
    //NT            result(INT_MAX);
    std::string   s;
    char          c = 0;
    
    //file.set_ascii_mode();
    while (file)
      {
        // try to convert next token to integer
        file >> c;
        
        if (c=='#') // comment
          {
            std::getline(file, s);
          }
        else  // an integer number.
          {
            file.putback(c);
            
            file >> num;
            return num;
          }
      }
    
    return num;
  }
  
  const Point& leftmost(const Point &p1, const Point &p2) const
  { 
    CGAL::Comparison_result rx = CGAL::compare_lexicographically_xy(p1, p2);
    
    if (rx == CGAL::SMALLER)
      return p1;
    else
      return p2;
  }

  const Point& rightmost(const Point &p1, const Point &p2) const
  { 
    CGAL::Comparison_result rx = CGAL::compare_lexicographically_xy(p1, p2);
    
    if (rx == CGAL::LARGER)
      return p1;
    else
      return p2;
  }
  
  Point find_rightmost_intersection(const Curve& cv1, 
                                    const Curve& cv2, 
                                    const Point& point) const
  {
    Traits traits;
    
    Point xp1, xp2, p=point;
    while (traits.nearest_intersection_to_right(cv1, cv2, p, xp1, xp2)){
      p = xp2;  // if there is no overlap - put xp1, else put xp2.
    }

    return xp2;
  }
  
  template <class Container>
  void find_intersections(const Curve& cv1, 
                          const Curve& cv2, 
                          Container& points)
  {
    Traits  traits;
    
    Point p = leftmost(leftmost(traits.curve_source(cv1),
                                traits.curve_target(cv1)), 
                       leftmost(traits.curve_source(cv2),
                                traits.curve_target(cv2)));
    
    Point xp1, xp2;
    while (traits.do_intersect_to_right(cv1, cv2, p))
      {
        traits.nearest_intersection_to_right(cv1, cv2, p, xp1, xp2);
        
        points.push_back(xp1);
        if (xp1 != xp2)
          points.push_back(xp2);
        
        p = xp2;  // if there is no overlap - put xp1, else put xp2.
      }
  }

  void read_file_build_creator(std::ifstream& file, Planar_map& pm)
  {
    Curve curr_curve;
    std::list<Curve> curves;
    Read_curve read_curve;
    // 1. read polylines and build arrangement
    
  
    // read number of curves
    unsigned int num_curves = get_next_int(file);
    
    //cout<<"num_curves="<<num_curves<<endl;
    
    // read curves (test specific)
    while (num_curves--) {
      curr_curve = read_curve(file);
      
      curves.push_back(curr_curve);
	//arr.insert(curr_curve);
    }

    Traits traits;
    sweep_to_construct_planar_map_2(curves.begin(),curves.end(), traits, pm);
  }

  /**********************************************************
   * Class validity functions.
   **********************************************************/
  void  check_that_vertices_are_in_map_overlay(const Subdivision& pm1, 
                                               const Subdivision& pm2, 
                                               const MapOverlay& map_overlay)
  {
    Vertex_const_iterator v_iter;
    for (v_iter = pm1.vertices_begin(); v_iter != pm1.vertices_end(); ++v_iter){
      Locate_type lt;
      Halfedge_const_handle  h=map_overlay.subdivision().locate(v_iter->point(), lt);
      
      CGAL_assertion(lt == Subdivision::VERTEX);
      
      // marking the vertices handle.
      if (h->source()->point() == v_iter->point())
        vertices[h->source().operator->()] = h->source();
      else
        vertices[h->target().operator->()] = h->target();
    }

    for (v_iter = pm2.vertices_begin(); v_iter != pm2.vertices_end(); ++v_iter)
      {
        Locate_type lt;
        Halfedge_const_handle  h=map_overlay.subdivision().locate(v_iter->point(), lt);
        CGAL_assertion(lt == Subdivision::VERTEX);
        
        // marking the vertices handle.
        if (h->source()->point() == v_iter->point())
          vertices[h->source().operator->()] = h->source();
        else
          vertices[h->target().operator->()] = h->target();
      }
    
    std::cout<<"check_that_vertices_are_in_map_overlay -- passed"<<std::endl;
  }

  void check_that_intersections_are_in_map_overlay(const Subdivision& pm1, 
                                                   const Subdivision& pm2, 
                                                   const MapOverlay& map_overlay)
  {
    for (Halfedge_const_handle h_iter1 = pm1.halfedges_begin(); 
         h_iter1 != pm1.halfedges_end(); ++h_iter1, ++h_iter1)
      for (Halfedge_const_handle h_iter2 = pm2.halfedges_begin(); 
           h_iter2 != pm2.halfedges_end(); ++h_iter2, ++h_iter2)
        {
          std::list<Point>  points;
          find_intersections(h_iter1->curve(), h_iter2->curve(), points);
          
          for (std::list<Point>::iterator p_iter = points.begin(); 
               p_iter != points.end(); ++p_iter){
            Locate_type lt;
            Halfedge_const_handle  h=map_overlay.subdivision().locate(*p_iter, lt);
            CGAL_assertion(lt == Subdivision::VERTEX);
            
            // marking the vertices handle.
            if (h->source()->point() == *p_iter)
              vertices[h->source().operator->()] = h->source();
            else
              vertices[h->target().operator->()] = h->target();
          }
        }
    
    std::cout<<"Check_that_intersections_are_in_map_overlay -- passed"<<std::endl;
  }
  
  void  check_that_halfedges_are_in_map_overlay(const Subdivision& pm1, 
                                                const Subdivision& pm2, 
                                                const MapOverlay& map_overlay)
  {
    check_that_halfedges_are_in_map_overlay(pm1,map_overlay);
    check_that_halfedges_are_in_map_overlay(pm2,map_overlay);

    std::cout<<"Check_that_halfedges_are_in_map_overlay -- passed"<<std::endl;
  }
  
  void  check_that_halfedges_are_in_map_overlay(const Subdivision& pm, 
                                                const MapOverlay& map_overlay)
  {
    Traits traits;
    
    for (Halfedge_const_iterator h_iter = pm.halfedges_begin();
         h_iter != pm.halfedges_end(); ++h_iter, ++h_iter)
      {
        // At each iteration we move on consecutive intersection points 
        // on h_iter->curve() from left to right.
        
#ifdef OVL_DEBUG_TEST
        std::cout<<"walking along "<<h_iter->curve() << std::endl;
#endif
      
        Point s = leftmost(h_iter->source()->point(),h_iter->target()->point());
        Point t = rightmost(h_iter->source()->point(),h_iter->target()->point());
        
        bool finish_traversal = false;
        
        while (!finish_traversal)
          {
            Locate_type lt;
            Halfedge_const_handle  h = 
              map_overlay.subdivision().locate(s,lt);
        
            CGAL_assertion(lt == Subdivision::VERTEX);
            
            Vertex_const_handle vertex = (h->source()->point() == s) ?
              h->source() : h->target();
            
            Halfedge_around_vertex_const_circulator  halfedge=vertex->incident_halfedges();
            
            // finding the halfdge in the overlay corresponding to the
            // original halfedge h_iter.
            
#ifdef OVL_DEBUG_TEST
            std::cout<<"--- vertex is "<<vertex->point()<<std::endl;
#endif
            do {
#ifdef OVL_DEBUG_TEST
              std::cout<<"walking around vertex: "<< halfedge->curve() <<std::endl;
#endif
              // finding the first halfedge which is part of h_iter->curve() and 
              // continues from the right of vertex->point().
              if (CGAL::compare_lexicographically_xy(
                                  rightmost(halfedge->source()->point(),
                                            halfedge->target()->point()), 
                                  vertex->point()) == CGAL::LARGER && 
                  traits.curves_overlap(halfedge->curve(), h_iter->curve()))
                break;
            } while (++halfedge != vertex->incident_halfedges());
            
#ifdef OVL_DEBUG_TEST
            std::cout<<"After walking around vertex, halfedge is "<< 
              halfedge->curve() <<std::endl;
#endif  
            CGAL_assertion(CGAL::compare_lexicographically_xy(
                                 rightmost(halfedge->source()->point(),
                                           halfedge->target()->point()), 
                                 vertex->point()) == CGAL::LARGER &&
                           traits.curves_overlap(halfedge->curve(), h_iter->curve()));
            
            halfedges[halfedge.operator->()] = halfedge;          // marking the corresponding
            halfedges[halfedge->twin().operator->()] = halfedge->twin();  // halfedge from the overlay.
            
            Point p1,p2;  // the intersection points: p1=s and p2 is the next intersection 
            // point on h_iter->curve();
            CGAL_assertion(traits.nearest_intersection_to_right(halfedge->curve(),
                                                                h_iter->curve(),s,p1,p2));
            
#ifdef OVL_DEBUG_TEST
            std::cout<<"p1 and p2  "<<p1<<" "<<p2<<std::endl; 
            std::cout<<"rightmost(p1,p2) "<<rightmost(p1,p2)<<std::endl;
            std::cout<<"leftmost(p1,p2) "<<leftmost(p1,p2)<<std::endl;
            std::cout<<"s and t are "<< s <<" "<<t <<std::endl;
#endif
            
            //       CGAL_assertion(p1 == s);
            
            Point next_point = find_rightmost_intersection(halfedge->curve(),
                                                           h_iter->curve(), s);
            
            
            // Here we assume that h_iter->curve() is x-monotone due to the fact
            // it is a part of Planar map.
            CGAL_assertion(CGAL::compare_lexicographically_xy(s,next_point) == CGAL::SMALLER);
            
            if (next_point == t)  // reaching the right endpoint of h_iter->curve().
              finish_traversal=true;
            else
              s = next_point; //rightmost(p1,p2);
          }
      }
  }
  
  void  check_that_faces_are_in_map_overlay(const MapOverlay& map_overlay)
  {
    const MapOverlay_change_notification*  notifier = map_overlay.change_notification();
    //Pm_file_writer<Subdivision>  pm_writer1(std::cout, pm1);
    //Pm_file_writer<Subdivision>  pm_writer2(std::cout, pm2);
    //Pm_file_writer<Subdivision>  pm_writer(std::cout, map_overlay.subdivision());
    
    for (Face_const_iterator f_iter = map_overlay.subdivision().faces_begin(); 
         f_iter != map_overlay.subdivision().faces_end(); ++f_iter)
      {
        //cout<<"f_iter"<<endl;
        //write_face(f_iter);
        
        Face_const_handle  face_creator1 = notifier->get_first_face_above(f_iter);
        Face_const_handle  face_creator2 = notifier->get_second_face_above(f_iter);
        
        faces[face_creator1.operator->()] = face_creator1;
        faces[face_creator2.operator->()] = face_creator2;
        
        //cout<<"face_creator1"<<endl;
        //write_face(face_creator1);
        
        //cout<<"face_creator2"<<endl;
        //write_face(face_creator2);
        // asserting that all halfedges along the holes of f_iter points to the same faces.
        for (Holes_const_iterator hit = f_iter->holes_begin(); 
             hit != f_iter->holes_end(); ++hit) 
          {
            Ccb_halfedge_const_circulator hole_halfedge(*hit);
            do {
              // for debugging
              //cout<<"notifier->get_first_face_above(hole_halfedge)"<<endl;
              //write_face(notifier->get_first_face_above(hole_halfedge));
              //cout<<"notifier->get_second_face_above(hole_halfedge)"<<endl;
              //write_face(notifier->get_second_face_above(hole_halfedge));
              
              CGAL_assertion(notifier->get_first_face_above(hole_halfedge) == 
                             face_creator1);
              CGAL_assertion(notifier->get_second_face_above(hole_halfedge) == 
                             face_creator2);
              // a queue to hold all faces involed with hit.
            } while (++hole_halfedge != *hit);
          }
        
        if (!f_iter->is_unbounded()){
          Ccb_halfedge_const_circulator face_halfedge = f_iter->outer_ccb();
          
          do {
            // for debugging
            //cout<<"notifier->get_first_face_above(hole_halfedge)"<<endl;
            //write_face(notifier->get_first_face_above(face_halfedge));
            //cout<<"notifier->get_second_face_above(hole_halfedge)"<<endl;
            //write_face(notifier->get_second_face_above(face_halfedge));
          
            CGAL_assertion(notifier->get_first_face_above(face_halfedge)==face_creator1);
            CGAL_assertion(notifier->get_second_face_above(face_halfedge)==face_creator2);
            // a queue to hold all faces involed with hit.
          } while (++face_halfedge != f_iter->outer_ccb());
        }
    }
    
    std::cout<<"Check_that_faces_are_in_map_overlay -- passed"<<std::endl;
  }

  void  check_all_features_are_marked(const MapOverlay& map_overlay)
  {
    // Checking vertices.
    for (Vertex_const_iterator v_iter = map_overlay.subdivision().vertices_begin();
         v_iter != map_overlay.subdivision().vertices_end(); ++v_iter)
      CGAL_assertion(vertices.find(&*v_iter) != vertices.end());
    
    // Checking halfedges.
    for (Halfedge_const_iterator h_iter = map_overlay.subdivision().halfedges_begin();
         h_iter != map_overlay.subdivision().halfedges_end(); ++h_iter)
      {
        CGAL_assertion(halfedges.find(&*h_iter) != halfedges.end());
      }
    
    // Checking faces of creators.
    Face_const_iterator f_iter;
    for (f_iter = map_overlay.first_creator()->subdivision().faces_begin();
         f_iter != map_overlay.first_creator()->subdivision().faces_end(); ++f_iter)
      CGAL_assertion(faces.find(&*f_iter) != faces.end());

    for (f_iter = map_overlay.second_creator()->subdivision().faces_begin();
         f_iter != map_overlay.second_creator()->subdivision().faces_end(); ++f_iter)
      CGAL_assertion(faces.find(&*f_iter) != faces.end());

    std::cout<<"check_all_features_are_marked -- passed"<<std::endl;
  }
  
  /////////////// debugging //////////////
  //         void write_face(Face_const_handle f) {
  
  //      std::cout<<"writing face"<<std::endl;
  
  //      std::cout<<"pointer of face="<<f.operator->()<<std::endl;
  
  //      if (f->is_unbounded()){
  //        std::cout<<"UNBOUNDED"<<std::endl;
  //        std::cout<<"number halfedges on outer boundary"<<std::endl;
  //        std::cout<<"0"<<std::endl;
  //      }
  //      else {
  //        std::cout<<"outer ccb"<<std::endl;
  
  //        Ccb_halfedge_const_circulator first = f->outer_ccb(), iter = first;
  
  //        std::size_t n = 0;
  //        do {
  //          std::cout<<iter->curve()<<" ";
  //          n++;
  //          iter++;
  //        } while (iter != first);
  
  //        std::cout<<"number halfedges on outer boundary"<<std::endl;
  //        std::cout<< n <<std::endl;
  
  //        std::cout << std::endl;
  //        }
  //          }
  
  /****************************
   * Class Interface
   ****************************/
  
public:
  void start(char * filename,
             CGAL::Map_overlay_base<Subdivision,MapOverlay_change_notification>* ovl_alg)
  {
    // Read data from file. Build Arrangement.
    std::ifstream file(filename);
      
    PmWalkPL     pl_walk1, pl_walk2;
    Planar_map   pm1(&pl_walk1), pm2(&pl_walk2);
    
    read_file_build_creator(file, pm1);
    read_file_build_creator(file, pm2);
    
    Subdivision  subdivision1(pm1), subdivision2(pm2);
    MapOverlay   first_creator(subdivision1, ovl_alg);
    MapOverlay   second_creator(subdivision2, ovl_alg);

    PmWalkPL pl_walk_ovl;
    //MapOverlay map_overlay(&pl_walk_ovl);
    //Subdivision pmwx1(pm1), pmwx2(pm2);
    
    //if (ovl_alg)
    MapOverlay map_overlay(first_creator, second_creator, 
                           &pl_walk_ovl, ovl_alg);
    
    //else
      //map_overlay = MapOverlay(pmwx1, pmwx2);
    // DEBUG
    //print_vertices(arr);

    // debug
    //       Arr_2::Face_handle    f  = arr->unbounded_face();
    //       Arr_2::Holes_iterator it = f->holes_begin(),end=f->holes_end();
    //       Arr_2::Ccb            c  = *it;
    //const X_curve& cv = curr->curve();
    
    
    // Check validity of arrangement after insertion
    CGAL_assertion(map_overlay.subdivision().is_valid());
    
    // Check that input vertices are indeed in the map overlay
    check_that_vertices_are_in_map_overlay(subdivision1, subdivision2, map_overlay);
    
     // Check that intersections vertices are indeed in the map overlay
    check_that_intersections_are_in_map_overlay(subdivision1, subdivision2, map_overlay);

    // Check that halfedges are indeed in the map overlay.
    check_that_halfedges_are_in_map_overlay(subdivision1, subdivision2, map_overlay);
    
    // Check that faces are indeed in the map overlay.
    check_that_faces_are_in_map_overlay(map_overlay); 

    check_all_features_are_marked(map_overlay);
  }  

  void start(char * filename)
  {
    // Read data from file. Build Arrangement.
    std::ifstream file(filename);
      
    PmWalkPL     pl_walk1, pl_walk2;
    Planar_map   pm1(&pl_walk1), pm2(&pl_walk2);
    
    read_file_build_creator(file, pm1);
    read_file_build_creator(file, pm2);
    
    Subdivision  subdivision1(pm1), subdivision2(pm2);
    MapOverlay   first_creator(subdivision1);
    MapOverlay   second_creator(subdivision2);
    
    PmWalkPL pl_walk_ovl;
    MapOverlay map_overlay(first_creator, second_creator, &pl_walk_ovl);
    
    
    // Check validity of arrangement after insertion
    CGAL_assertion(map_overlay.subdivision().is_valid());
    
    // Check that input vertices are indeed in the map overlay
    check_that_vertices_are_in_map_overlay(subdivision1, subdivision2, map_overlay);
    
     // Check that intersections vertices are indeed in the map overlay
    check_that_intersections_are_in_map_overlay(subdivision1, subdivision2, map_overlay);

    // Check that halfedges are indeed in the map overlay.
    check_that_halfedges_are_in_map_overlay(subdivision1, subdivision2, map_overlay);
    
    // Check that faces are indeed in the map overlay.
    check_that_faces_are_in_map_overlay(map_overlay); 

    check_all_features_are_marked(map_overlay);
  }  
};


#endif









