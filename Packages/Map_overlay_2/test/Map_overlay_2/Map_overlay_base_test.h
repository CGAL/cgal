#ifndef  MAP_OVERLAY_BASE_TEST
#define  MAP_OVERLAY_BASE_TEST

#include <CGAL/Cartesian.h>
//#include <CGAL/Arr_2_bases.h>
//#include <CGAL/Arr_2_default_dcel.h>
#include <fstream>

//#ifndef CGAL_PLANAR_MAP_2
//#include <CGAL/Planar_map_2.h>
//#endif

#ifndef CGAL_ARR_2_OVERLAY_DCEL_H
#include <CGAL/Map_overlay_default_dcel.h>
#endif

#ifndef CGAL_MAP_OVERLAY_POST_PROC_NOTIFIER_H
#include <CGAL/Map_overlay_post_proc_notifier.h>
#endif

#ifndef CGAL_MAP_OVERLAY_H
#include <CGAL/Map_overlay.h>
#endif
 
//#include <CGAL/Arrangement_2.h>

#include <CGAL/sweep_to_construct_planar_map.h>

//for debugging
#include <CGAL/IO/Pm_file_writer.h>

// Quotient is included anyway, because it is used to read
// data files. Quotient can read both integers and fractions.
// leda rational will only read fractions.
#include <CGAL/Quotient.h> 

#include <list>
#include <string>


// we use the namespace std for compatability with MSVC

template <class Subdivision, class Read_curve>
class Map_overlay_base_test
{
protected:
  typedef typename Subdivision::Traits                  Traits;
  
  typedef typename Traits::Point                                 Point;
  typedef typename Traits::X_curve                               X_curve;
  typedef typename Traits::Curve                                 Curve;
  
  typedef CGAL::Map_overlay_default_dcel<Traits>                Dcel;
  
  typedef CGAL::Planar_map_2<Dcel, Traits>                      PM;
  
  typedef typename Subdivision::Vertex                                     Vertex;
  typedef typename Subdivision::Face                                       Face;
  typedef typename Subdivision::Halfedge                                   Halfedge;
  typedef typename Subdivision::Vertex_handle                              Vertex_handle;
  typedef typename Subdivision::Halfedge_handle                            Halfedge_handle;
  typedef typename Subdivision::Face_handle                                Face_handle;
  typedef typename Subdivision::Vertex_const_handle                        Vertex_const_handle;
  typedef typename Subdivision::Halfedge_const_handle                      Halfedge_const_handle;
  typedef typename Subdivision::Face_const_handle                          Face_const_handle;
  typedef typename Subdivision::Vertex_iterator                            Vertex_iterator;
  typedef typename Subdivision::Vertex_const_iterator                      Vertex_const_iterator;
  typedef typename Subdivision::Halfedge_iterator                          Halfedge_iterator;
  typedef typename Subdivision::Halfedge_const_iterator                    Halfedge_const_iterator;
  typedef typename Subdivision::Face_iterator                              Face_iterator;
  typedef typename Subdivision::Face_const_iterator                        Face_const_iterator;
  typedef typename Subdivision::Ccb_halfedge_circulator                    Ccb_halfedge_circulator;
  typedef typename Subdivision::Ccb_halfedge_const_circulator              Ccb_halfedge_const_circulator;
  typedef typename Subdivision::Holes_iterator                             Holes_iterator;
  typedef typename Subdivision::Holes_const_iterator                       Holes_const_iterator;
  typedef typename Subdivision::Locate_type                                Locate_type;
  typedef typename Subdivision::Traits_wrap                                Traits_wrap;
  
  //typedef CGAL::Arr_base_node<Curve>                     Base_node;
  //typedef CGAL::Map_overlay_default_dcel<Traits, 
  //  CGAL::Arr_2_vertex_base<Point>, 
  //  CGAL::Arr_2_halfedge_base<Base_node>, 
  //  CGAL::Arr_2_face_base>   Arr_Dcel;
  //typedef CGAL::Arrangement_2<Arr_Dcel,Traits,Base_node>        Arrangement;
  //typedef Arrangement::Change_notification                      PMWXChangeNotification; 
  
  //typedef Arrangement::Halfedge_iterator                        Arr_halfedge_iterator;
  
  typedef CGAL::Pm_naive_point_location<Subdivision>                       PmNaivePL;
  typedef CGAL::Pm_walk_along_line_point_location<Subdivision>             PmWalkPL;
  
  typedef CGAL::Map_overlay_post_proc_notifier<Subdivision>                MapOverlay_change_notification;
  typedef CGAL::Map_overlay<Subdivision,MapOverlay_change_notification>    MapOverlay; 
  
  //MapOverlay  map_overlay;  
  std::map<const void*, Vertex_const_handle>    vertices;
  std::map<const void*, Halfedge_const_handle>  halfedges;
  // std::map<void*, Face_handle>      faces;
  
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

  template <class Container>
  void find_intersections(const Curve& cv1, const Curve& cv2, Container& points)
  {
    Traits  traits;
    
    Point p = leftmost(leftmost(traits.curve_source(cv1),traits.curve_target(cv1)), 
                        leftmost(traits.curve_source(cv2),traits.curve_target(cv2)));
    
    //cout<<"p="<<p<<endl;
    
    Point xp1, xp2;
    while (traits.nearest_intersection_to_right(cv1, cv2, p, xp1, xp2)){
      //cout<<"xp1="<<xp1<<endl;
      points.push_back(xp1);
      if (xp1 != xp2)
        points.push_back(xp2);
      
      p = xp2;  // if there is no overlap - put xp1, else put xp2.
    }
  }

  void read_file_build_creator(std::ifstream& file, Subdivision& pm)
  {
    Curve curr_curve;
    std::list<Curve> curves;
    Read_curve read_curve;
    // 1. read polylines and build arrangement
    
  
    // read number of curves
    unsigned int num_curves = get_next_int(file);
    
    // read curves (test specific)
    while (num_curves--) {
      curr_curve = read_curve(file);
      
      curves.push_back(curr_curve);
	//arr.insert(curr_curve);
    }

    sweep_to_construct_planar_map(curves.begin(),curves.end(), pm);
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

    for (v_iter = pm2.vertices_begin(); v_iter != pm2.vertices_end(); ++v_iter){
      Locate_type lt;
      Halfedge_const_handle  h=map_overlay.subdivision().locate(v_iter->point(), lt);
      CGAL_assertion(lt == Subdivision::VERTEX);

      // marking the vertices handle.
      if (h->source()->point() == v_iter->point())
        vertices[h->source().operator->()] = h->source();
      else
        vertices[h->target().operator->()] = h->target();
    }
    
    cout<<"check_that_vertices_are_in_map_overlay -- passed"<<endl;
  }

  void check_that_intersections_are_in_map_overlay(const Subdivision& pm1, 
                                                   const Subdivision& pm2, 
                                                   const MapOverlay& map_overlay)
  {
    for (Halfedge_const_handle h_iter1 = pm1.halfedges_begin(); 
         h_iter1 != pm1.halfedges_end(); ++h_iter1, ++h_iter1)
      for (Halfedge_const_handle h_iter2 = pm2.halfedges_begin(); 
           h_iter2 != pm2.halfedges_end(); ++h_iter2, ++h_iter2){
        std::list<Point>  points;
        find_intersections(h_iter1->curve(), h_iter2->curve(), points);
        
        for (std::list<Point>::iterator p_iter = points.begin(); p_iter != points.end(); ++p_iter){
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

    cout<<"Check_that_intersections_are_in_map_overlay -- passed"<<endl;
  }

  void  check_that_halfedges_are_in_map_overlay(const Subdivision& pm1, 
                                                const Subdivision& pm2, 
                                                const MapOverlay& map_overlay)
  {
  }

  void  check_that_faces_are_in_map_overlay(const MapOverlay& map_overlay)
  {
    const MapOverlay_change_notification*  notifier = map_overlay.change_notification();
    //Pm_file_writer<Subdivision>  pm_writer1(std::cout, pm1);
    //Pm_file_writer<Subdivision>  pm_writer2(std::cout, pm2);
    //Pm_file_writer<Subdivision>  pm_writer(std::cout, map_overlay.subdivision());
    
    for (Face_const_iterator f_iter = map_overlay.subdivision().faces_begin(); 
         f_iter != map_overlay.subdivision().faces_end(); ++f_iter){
      //cout<<"f_iter"<<endl;
      //write_face(f_iter);
      
      Face_const_handle  face_creator1 = notifier->get_first_face_above(f_iter);
      Face_const_handle  face_creator2 = notifier->get_second_face_above(f_iter);
      
      //cout<<"face_creator1"<<endl;
      //write_face(face_creator1);

      //cout<<"face_creator2"<<endl;
      //write_face(face_creator2);
      // asserting that all halfedges along the holes of f_iter points to the same faces.
      for (Holes_const_iterator hit = f_iter->holes_begin(); 
           hit != f_iter->holes_end(); ++hit) {
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
          
          CGAL_assertion(notifier->get_first_face_above(face_halfedge) == face_creator1);
          CGAL_assertion(notifier->get_second_face_above(face_halfedge) == face_creator2);
          // a queue to hold all faces involed with hit.
        } while (++face_halfedge != f_iter->outer_ccb());
      }
    }

    cout<<"Check_that_faces_are_in_map_overlay -- passed"<<endl;
  }

  void  check_all_freatures_are_marked(const MapOverlay& map_overlay)
  {
    for (Vertex_const_iterator v_iter = map_overlay.subdivision().vertices_begin();
         v_iter != map_overlay.subdivision().vertices_end(); ++v_iter)
      CGAL_assertion(vertices.find(&*v_iter) != vertices.end());

    //for (Halfedge_const_iterator h_iter = map_overlay.subdivision().halfedges_begin();
    //     h_iter != map_overlay.subdivision().halfedges_end(); ++h_iter)
    //  CGAL_assertion(halfedges.find(&*h_iter) != halfedges.end());
    
    cout<<"check_all_freatures_are_marked -- passed"<<endl;
  }

  /**** debugging ***
        void write_face(Face_const_handle f) {
    
    std::cout<<"writing face"<<std::endl;
    
    std::cout<<"pointer of face="<<f.operator->()<<std::endl;
    
    if (f->is_unbounded()){
      std::cout<<"UNBOUNDED"<<std::endl;
      std::cout<<"number halfedges on outer boundary"<<std::endl;
      std::cout<<"0"<<std::endl;
    }
    else {
      std::cout<<"outer ccb"<<std::endl;
      
      Ccb_halfedge_const_circulator first = f->outer_ccb(), iter = first;

      std::size_t n = 0;
      do {
        std::cout<<iter->curve()<<" ";
        n++;
        iter++;
      } while (iter != first);

      std::cout<<"number halfedges on outer boundary"<<std::endl;
      std::cout<< n <<std::endl;
      
      std::cout << std::endl;
      }
      }*/
        
  /****************************
   * Class Interface
   ****************************/
public:
  void start(char * filename1, char * filename2, 
             CGAL::Map_overlay_base<Subdivision,MapOverlay_change_notification>* ovl_alg=0)
  {
    // Read data from file. Build Arrangement.
    std::ifstream file1(filename1);
    std::ifstream file2(filename2);
      
    PmWalkPL pl_walk1, pl_walk2;
    Subdivision   pm1(&pl_walk1), pm2(&pl_walk2);
    
    read_file_build_creator(file1, pm1);
    read_file_build_creator(file2, pm2);
    
    //MapOverlay first_creator(pm1);
    //MapOverlay second_creator(pm2);
    
    PmWalkPL pl_walk_ovl;
    MapOverlay map_overlay(&pl_walk_ovl);
    
    if (ovl_alg)
      map_overlay = MapOverlay(pm1, pm2, ovl_alg);
    else
      map_overlay = MapOverlay(pm1, pm2);
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
    check_that_vertices_are_in_map_overlay(pm1, pm2, map_overlay);
    
     // Check that intersections vertices are indeed in the map overlay
    check_that_intersections_are_in_map_overlay(pm1, pm2, map_overlay);

    // Check that halfedges are indeed in the map overlay.
    //    check_that_halfedges_are_in_map_overlay(pm1, pm2, map_overlay);
    
    // Check that faces are indeed in the map overlay.
    check_that_faces_are_in_map_overlay(map_overlay); 

    check_all_freatures_are_marked(map_overlay);
  }  
 
};


#endif















