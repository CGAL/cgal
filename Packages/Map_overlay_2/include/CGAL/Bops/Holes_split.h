#ifndef CGAL_HOLES_SPLIT_H
#define CGAL_HOLES_SPLIT_H

//#ifndef LEDA_RAT_RAY_H
//#include <LEDA/rat_ray.h>
//#endif


#ifndef CGAL_RAY_2_SEGMENT_2_INTERSECTION_H
#include <CGAL/Ray_2_Segment_2_intersection.h>
#endif

#ifndef CGAL_ASSERTIONS_H
#include <CGAL/assertions.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class Planar_map_, class Notifier_>
class Holes_split {
  //typedef enum { EPSILON=0.01 }; 
  static const double EPSILON=0.01;  //change it to work on VC++.
  
  typedef  Planar_map_    Planar_map;
  typedef  Notifier_      Notifier;
  
  typedef typename Planar_map_::Traits        Traits;
  typedef typename Traits::Point              Point;
  typedef typename Traits::X_curve            X_curve;
  typedef typename Traits::Curve              Curve;

  typedef typename Planar_map::Vertex                 Vertex;
  typedef typename Planar_map::Face                   Face;
  typedef typename Planar_map::Halfedge               Halfedge;
  typedef typename Planar_map::Vertex_handle          Vertex_handle;
  typedef typename Planar_map::Halfedge_handle        Halfedge_handle;
  typedef typename Planar_map::Face_handle            Face_handle;
  typedef typename Planar_map::Vertex_const_handle    Vertex_const_handle;
  typedef typename Planar_map::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Planar_map::Face_const_handle      Face_const_handle;
  typedef typename Planar_map::Vertex_iterator        Vertex_iterator;
  typedef typename Planar_map::Vertex_const_iterator  Vertex_const_iterator;
  typedef typename Planar_map::Halfedge_iterator      Halfedge_iterator;
  typedef typename Planar_map::Halfedge_const_iterator   Halfedge_const_iterator;
  typedef typename Planar_map::Face_iterator          Face_iterator;
  typedef typename Planar_map::Face_const_iterator    Face_const_iterator;
  typedef typename Planar_map::Ccb_halfedge_circulator    Ccb_halfedge_circulator;
  typedef typename Planar_map::Ccb_halfedge_const_circulator   Ccb_halfedge_const_circulator;
  typedef typename Planar_map::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
  typedef typename Planar_map::Holes_iterator          Holes_iterator;
  typedef typename Planar_map::Holes_const_iterator    Holes_const_iterator;
  typedef typename Planar_map::Locate_type             Locate_type;
  typedef typename Planar_map::Traits_wrap             Traits_wrap;
  
  typedef typename Point::R                    Kernel;
  typedef Ray_2<Kernel>                        Ray;
  //typedef leda_rat_ray                                      Ray;
  
  class less_xy_Vertex_handle {
  public:
    
    inline  bool operator()(Vertex_handle v1, Vertex_handle v2) const {
      return CGAL::compare_lexicographically_xy(v1->point(),v2->point()) == CGAL::SMALLER ;
    }
  };
  
  
  bool has_decomposing_edge(Vertex_handle v)
  {
    Halfedge_around_vertex_circulator circ = v->incident_halfedges();
    
    do {
      if (circ->decomposing())
        return true;
    } while(++circ!=v->incident_halfedges());

    return false;
  }
  
  template <class OutputIterator>
  bool get_decomposing_edge(Vertex_handle v, OutputIterator halfedges)
  {
    bool b=false;
    Halfedge_around_vertex_circulator circ = v->incident_halfedges();
    
    do {
      if (circ->decomposing()){
        b=true;
        *halfedges++ = circ;
      }
    } while(++circ!=v->incident_halfedges());
    
    return b;
  }
  
  void add_vertical_curve(Vertex_handle v, Planar_map& pm, Notifier& notf, bool up)
  {
    // debugging.
    //cout<< "In add_vertical_curve" <<endl;
    //cout<< v->point().xcoordD()<<","<<v->point().ycoordD() << endl;
    
    Locate_type lt;
    double eps = up? EPSILON: -EPSILON;
    Point pertrubed_p = Point(v->point().x(), 
                              v->point().y() + eps);
    //cout<<"pertrubed_p="<< pertrubed_p.xcoordD()<<","<<pertrubed_p.ycoordD() << endl;
    
    // we have to raise p since the function
    // vertical_ray_shoot will return the vertex of of in
    // pmwx (that's how it is defined when we shoot a vertical ray
    // from a feature that exists in the map.
    Halfedge_handle h = pm.vertical_ray_shoot(pertrubed_p, lt, up);
    
    //cout <<"h->curve()="<<h->curve()<<endl;
    
    if (h->face()->is_unbounded())
      return;  // what we should do here is using the bbox and adding 
               // a vertical curve upward (downward) will it hits the bbox.
    
    Ray ray(v->point(), pertrubed_p);
    Point hitting_point;
  
    Object obj=intersection(h->curve(), ray);
    if (!assign(hitting_point, obj))
      return;
    
    //if (!ray.intersection(h->curve(), hitting_point)) // only for LEDA
    //  return;
    //cout<<"hitting_point="<<hitting_point.xcoordD()<<","<< hitting_point.ycoordD() <<endl;
    
    notf.set_decomposing_edge(true);
    //pm.non_intersecting_insert_from_vertex(Curve(v->point(),hitting_point),
    //                                         v,up,&notf);
    
    Traits traits;
    X_curve c1,c2;
    
    traits.curve_split(h->curve(), c1, c2, hitting_point);
    Halfedge_handle split_h = pm.split_edge(h,c1,c2,&notf);
    CGAL_assertion(v->point() != hitting_point);
    pm.insert(Curve(v->point(),hitting_point),&notf);
    
    notf.set_decomposing_edge(false);
  }
  
  void  set_faces_in_holes(Planar_map& pm)
  {
    for (Face_iterator f_iter=pm.faces_begin(); 
         f_iter!=pm.faces_end(); ++f_iter){
      
      if (f_iter->is_unbounded())
        continue;
      
      for (Holes_iterator holes_it = f_iter->holes_begin(); 
           holes_it != f_iter->holes_end(); ++holes_it)  
        // set the face inside the hole to true.
        (*holes_it)->twin()->face()->set_in_hole(true); 
    }
  }
  
public:
  // This function splits the holes in pm.
  // Here, pm has a special structure of one face with possibly many holes.
  void  split_holes(Planar_map& pm, 
                    Notifier& notf)
  {    
    set_faces_in_holes(pm);
      
    for (Face_iterator f_iter=pm.faces_begin(); 
         f_iter!=pm.faces_end(); ++f_iter){
      Face_handle face= f_iter;
      
      if (face->is_unbounded())
        continue; // splitting only bounded faces.
      
      // keeping all Ccb in a seperate vector, since face is going to be changed.
      std::vector<Ccb_halfedge_circulator> holes(face->holes_begin(),face->holes_end());
      
      //for (Holes_iterator holes_it = face->holes_begin(); 
      //     holes_it != face->holes_end(); ++holes_it) {
      for (std::vector<Ccb_halfedge_circulator>::iterator holes_it = holes.begin(); 
           holes_it != holes.end(); ++holes_it) {
        Ccb_halfedge_circulator face_cc(*holes_it);
        
        // finding the leftmost and the rightmost points of the hole.
        std::vector<Vertex_handle> vertices_cc; 
        do {
          vertices_cc.push_back(face_cc->source());
        } while (++face_cc != *holes_it);
        
        std::sort(vertices_cc.begin(), vertices_cc.end(), less_xy_Vertex_handle() );
        
        Vertex_handle leftlow_most = vertices_cc[0];
        Vertex_handle lefthigh_most = vertices_cc[0];
        
        // updating lefthigh_most to have the highest value of y.
        for (unsigned int i=0; i < vertices_cc.size() - 1 &&
             vertices_cc[i]->point().x() == vertices_cc[i+1]->point().x(); ++i)
          lefthigh_most=vertices_cc[i+1];
        
        Vertex_handle righthigh_most = vertices_cc.back();
        Vertex_handle rightlow_most = vertices_cc.back();
        
        // updating rightlow_most to have the lowest value of y.
        for (unsigned int i=vertices_cc.size() - 1;  i > 1 &&
             vertices_cc[i]->point().x() == vertices_cc[i-1]->point().x(); --i)
          lefthigh_most=vertices_cc[i-1];

        /*Vertex_handle leftlow_most=face_cc->source();
          Vertex_handle lefthigh_most=face_cc->source();
          Vertex_handle rightlow_most=face_cc->source();
          Vertex_handle righthigh_most=face_cc->source();
          do {
          if (face_cc->source()->point().xcoord() < leftmost->point().xcoord())
          lefthigh_most=leftlow_most=face_cc->source();
          else if (face_cc->source()->point().xcoord() == leftmost->point().xcoord()){
          
          }
          if (face_cc->source()->point().xcoord() > rightmost->point().xcoord())
          rightmost=face_cc->source();
          } while (++face_cc != *holes_it);*/
        
        bool b_lefthigh=has_decomposing_edge(lefthigh_most);
        bool b_leftlow=has_decomposing_edge(leftlow_most);
        bool b_righthigh=has_decomposing_edge(righthigh_most);
        bool b_rightlow=has_decomposing_edge(rightlow_most);
          
        if (!b_lefthigh)
          add_vertical_curve(lefthigh_most, pm, notf, true);
        
        //bool b_leftlow=has_decomposing_edge(leftlow_most);
        if (!b_leftlow)
          add_vertical_curve(leftlow_most, pm, notf, false);
        
        //bool b_righthigh=has_decomposing_edge(righthigh_most);
        if (!b_righthigh)
          add_vertical_curve(righthigh_most, pm, notf, true);
        
        //bool b_rightlow=has_decomposing_edge(rightlow_most);
        if (!b_rightlow)
          add_vertical_curve(rightlow_most, pm, notf, false);
      }
    }
  } 
};

CGAL_END_NAMESPACE

#endif






