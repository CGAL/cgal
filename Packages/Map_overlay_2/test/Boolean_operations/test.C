#include <CGAL/config.h> // needed for the LONGNAME flag

#include <CGAL/Cartesian.h>
#include <CGAL/leda_rational.h>

#ifndef CGAL_PLANAR_MAP_2
#include <CGAL/Planar_map_2.h>
#endif

//#include <CGAL/Pm_with_intersections.h>

#include <CGAL/Arr_2_bases.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Arrangement_2.h>

#ifndef CGAL_BOP_DEFAULT_DCEL_H
#include <CGAL/Bop_default_dcel.h>
#endif

#ifndef CGAL_MAP_OVERLAY_DEFAULT_NOTIFIER_H
#include <CGAL/Map_overlay_default_notifier.h>
#endif

#ifndef CGAL_MAP_OVERLAY_H
#include <CGAL/Map_overlay.h>
#endif

#ifndef BOOLEAN_OPERATIONS_2_H
#include <CGAL/Boolean_operations_2.h>
#endif

#include <CGAL/Polygon_2.h>
#include <LEDA/polygon.h>
 
#include <CGAL/IO/Arr_iostream.h>

#include <CGAL/sweep_to_construct_planar_map_2.h>

// Quotient is included anyway, because it is used to read
// data files. Quotient can read both integers and fractions.
// leda rational will only read fractions.
#include <CGAL/Quotient.h> 

#include <list>
#include <string>
#include <fstream>

// Making sure test doesn't fail if LEDA is not installed
#if ! defined(CGAL_USE_LEDA) && \
      (CGAL_ARR_TEST_TRAITS == CGAL_POLYLINE_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_LEDA_TRAITS || \
       CGAL_ARR_TEST_TRAITS == CGAL_SEGMENT_CIRCLE_TRAITS)

int main(int argc, char* argv[])
{
  std::cout << "A try to run test with LEDA traits but LEDA is not installed.";
  std::cout << std::endl;
  std::cout << "Test is not performed.";
  std::cout << std::endl;

  return 0;
}
#else

// Choose traits
#include <CGAL/Arr_leda_segment_exact_traits.h>

// Quotient is included anyway, because it is used to read
// data files. Quotient can read both integers and fractions.
// leda rational will only read fractions.
//nclude <CGAL/Quotient.h> 


class Boolean_operations_test {
  struct Face_with_counter : public CGAL::Arr_2_face_base {
    Face_with_counter() : CGAL::Arr_2_face_base(), counter(-1) {}
    int counter;
  };
  
public:
  //typedef CGAL::Quotient<int>                NT;
  typedef leda_rational                      NT;
  typedef CGAL::Cartesian<NT>                R;
  typedef CGAL::Arr_segment_exact_traits<R>  Traits;
  //typedef CGAL::Arr_segment_exact_traits  Traits;
  typedef Traits::Point                      Point;
  typedef Traits::X_curve                    X_curve;
  typedef Traits::Curve                      Curve;

  typedef CGAL::Polygon_traits_2<R>          Polygon_traits;
  typedef std::list<Point>                   Polygon_Container;
  typedef CGAL::Polygon_2<Polygon_traits, Polygon_Container > Polygon;
  typedef std::list<Polygon>                 Polygon_list;

  typedef CGAL::Bop_default_dcel<Traits>                 Dcel;
  typedef CGAL::Planar_map_2<Dcel, Traits>               Planar_map;
  typedef CGAL::Map_overlay_default_notifier<Planar_map> 
                                             MapOverlay_change_notification;
  typedef CGAL::Map_overlay<Planar_map,MapOverlay_change_notification>  MapOverlay; 
  typedef CGAL::Boolean_operations_2<MapOverlay>                        Bops;
  
  typedef Planar_map::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
  typedef Planar_map::Ccb_halfedge_const_circulator   Ccb_halfedge_const_circulator;
  typedef Planar_map::Face_handle               Face_handle;
  typedef Planar_map::Vertex_iterator           Vertex_iterator;
  typedef Planar_map::Face_iterator             Face_iterator;
  typedef Planar_map::Holes_iterator            Holes_iterator;
  typedef Planar_map::Holes_const_iterator            Holes_const_iterator;
  typedef Planar_map::Face                      Face;
  typedef Planar_map::Halfedge                  Halfedge;
  typedef Planar_map::Vertex                    Vertex;
  
  typedef Bops::Faces_container                      Faces_container;
  typedef Bops::Halfedges_container                  Halfedges_container;
  typedef Bops::Vertices_container                   Vertices_container;
  
  //  typedef CGAL::Planar_map_with_intersections_2<Planar_map>   Pmwx;
  
  typedef CGAL::Arr_base_node<X_curve>          Base_node;
  typedef CGAL::Pm_dcel<CGAL::Arr_2_vertex_base<Point>,
                        CGAL::Arr_2_halfedge_base<Base_node >,
                        Face_with_counter >  Arr_Dcel;
  typedef CGAL::Arrangement_2<Arr_Dcel,Traits,Base_node >      Arrangement;

  typedef CGAL::Pm_walk_along_line_point_location<Planar_map>  PmWalkPL;
  typedef CGAL::Pm_naive_point_location<Planar_map>            PmNaivePL;
  
  
private:
  template <class Point>
  class less_xy {
  public:
    
    inline  bool operator()(const Point& p1, const Point& p2) const 
    { 
      return CGAL::compare_lexicographically_xy(p1,p2) == CGAL::SMALLER;
    }
  };
  
  //------------------------------------ Boolean operations functions.
  
  //generalized face_diff function, to acount for overlaps.
  int face_diff (Arrangement::Ccb_halfedge_circulator circ) {
    Traits t;
    int diff = 0;
    Arrangement::Overlap_circulator oc = circ->overlap_edges();
    do {
      if (circ->source()->point() == t.curve_source(oc->curve()) ) 
        diff--;     //we're inside, going outside
      else
        diff++;
    } while (++oc != circ->overlap_edges());

    return diff;
  } 

  // covering_DFS will compute for each face in how many polygons it is.
  // It is a recursive DFS function and will be called with the unbounded 
  // face after its counter has been initialized to 0.
  void covering_DFS(Arrangement::Face_handle f) {
    Arrangement::Ccb_halfedge_circulator start,circ;
    
    // Do a recursive step for all neighbours, if any exists.
    if (f->does_outer_ccb_exist()) {
      start = circ = f->outer_ccb();
      do {
        if (circ->twin()->face()->counter < 0) {
          int diff = face_diff(circ);
          circ->twin()->face()->counter = (f->counter + diff); 
          covering_DFS(circ->twin()->face());
        }
      } while (++circ != start);
    }
    
    // Do a recursive step for all holes, if any exists.
    Arrangement::Holes_iterator hit = f->holes_begin();
    for (; Arrangement::Holes_iterator(hit)!=
           Arrangement::Holes_iterator(f->holes_end()); ++hit) {
      start = circ = (*hit);
      do {
        if (circ->twin()->face()->counter < 0) {
          int diff = face_diff(circ);
          circ->twin()->face()->counter = (f->counter + diff); 
          covering_DFS(circ->twin()->face());        
        }
      } while (++circ != start);
    }
  } 
  
  // Construct the arrangement that will use for calculating the intersection.
  void insert_polygons(Arrangement &arr, Polygon_list &in_poly_list)
  {
    Polygon::Edge_const_iterator it;
    Arrangement::Curve_iterator         ci;
    Polygon_list::iterator       plit;
    
    // for each polygon in list
    for (plit = in_poly_list.begin(); plit != in_poly_list.end(); ++plit) {

      // Make sure polygons are oriented counterclockwise
      // to satisfy assumption in DFS function
      if ( ! plit->is_counterclockwise_oriented())
        plit->reverse_orientation();
      
      // insert polygon to arrangement
      for (it = plit->edges_begin(); it != plit->edges_end(); ++it)
      {
	ci = arr.insert(*it);
      }
    }
  }

  /*
  // Convert faces of the arrangement that are in the intersection
  // to polygons.
  void polygons_from_faces(Arrangement& arr,
                           std::list<Face_iterator>& face_it_list,
                           Polygon_list& poly_list)
  {
    std::list<Face_iterator>::iterator  lit;
    //Pmwx::Ccb_halfedge_circulator cc;
    Polygon                        poly;
    
    for (lit = face_it_list.begin(); lit != face_it_list.end(); ++lit) {
      
      poly.erase(poly.vertices_begin(), poly.vertices_end());
      Arrangement::Ccb_halfedge_circulator cc=(*lit)->outer_ccb();
      do {
        poly.push_back(cc->curve().source());
      } while (++cc != (*lit)->outer_ccb());
      poly_list.push_back(poly);
    }
    }*/
  
  // performs the extraction of data out of the processed arrangement
  // if covering = 0, will perform union
  // otherwise, if there are n polygons in the arrangement and covering == n
  // then will perform intersection
  template <class OutputIterator>
  void get_faces_with_covering(OutputIterator faces, 
                               Arrangement& arr,
                               int covering)
  {
    Arrangement::Face_handle uf = arr.unbounded_face();
    uf->counter = 0;
    covering_DFS(uf);
    
    //"collecting" the union boundary faces. 
    for(Arrangement::Face_iterator fit = arr.faces_begin(); 
        fit!=arr.faces_end(); ++fit) {
      
      // if the face is covered by 'covering' 
      if (fit->counter == covering) {
        //unions.push_back(fit);
        *faces = fit;
        ++faces;
      }
    }
  }   
  
   template <class InputIterator1, class InputIterator2, class Container>
  void get_points_below_vertex_with_covering(InputIterator1 points_begin,
                                             InputIterator1 points_end,
                                             InputIterator2 polygons_begin, 
                                             InputIterator2 polygons_end,
                                             Container& points, 
                                             int covering)
  {    
    for (InputIterator1 p_iter=points_begin; p_iter!= points_end; ++p_iter){
      int num_vertices_above=0;
      for (InputIterator2 poly_iter=polygons_begin; poly_iter!=polygons_end; ++poly_iter){
        bool found=false;
        for (Polygon::Vertex_const_iterator v_iter=poly_iter->vertices_begin();
             v_iter != poly_iter->vertices_end() && !found; ++v_iter){  // *e_iter is CGAL Segment.
          if (*v_iter == *p_iter)
            found=true;
        }
        if (found)
          ++num_vertices_above;
      }
      
      //cout << *p_iter <<" is below "<< num_vertices_above <<" vertices"<<endl;
      
      if (num_vertices_above==covering)
        points.insert(*p_iter);
      else {    // removing that point from points, if it exists
        typename Container::iterator iter = points.find(*p_iter);
        if (iter != points.end())
          points.erase(iter);
      }
    }
  }
  
  template <class InputIterator1, class InputIterator2, class Container>
  void get_points_below_edge_with_covering(InputIterator1 points_begin,
                                           InputIterator1 points_end,
                                           InputIterator2 polygons_begin, 
                                           InputIterator2 polygons_end,
                                           Container& points, 
                                           int covering)
  {
    Traits traits;
    
    for (InputIterator1 p_iter=points_begin; p_iter!= points_end; ++p_iter){
      int num_edges_above=0;
      for (InputIterator2 poly_iter=polygons_begin; poly_iter!=polygons_end; ++poly_iter){
        bool found=false;
        for (Polygon::Edge_const_iterator e_iter=poly_iter->edges_begin();
             e_iter != poly_iter->edges_end() && !found; ++e_iter){  // *e_iter is CGAL Segment.
          Curve cv(*e_iter);
          if (traits.curve_get_point_status(cv,*p_iter) == Traits::ON_CURVE)
            found=true;
        }
        if (found)
          ++num_edges_above;
      }
      
      //cout << *p_iter <<" is below "<< num_edges_above <<" edges"<<endl;
      
      if (num_edges_above==covering)
        points.insert(*p_iter);
      else {    // removing that point from points, if it exists
        typename Container::iterator iter = points.find(*p_iter);
        if (iter != points.end())
          points.erase(iter);
      }
    }
  }
  
  template <class InputIterator, class Container>
  void get_points_below_halfedge(InputIterator halfedges_begin, 
                                 InputIterator halfedges_end,
                                 Container& points)
  {
    InputIterator iter=halfdges_begin;
    for ( ;halfedge_iter!=halfedge_end; ++halfedge_iter,  ++halfedge_iter){
      points.insert(halfedge_iter->source()->point());
      points.insert(halfedge_iter->target()->point());
    }
    
  }
  
  // we can use this function, instead of computing exact covering, since 
  // there are not isolated vertices in our structures. Hence if a vertex is below 
  // a face - it has to be connected to at least one edge.
  template <class InputIterator1, class InputIterator2,class Container>
  void get_points_below_faces(InputIterator1 points_begin,
                              InputIterator1 points_end,
                              InputIterator2 faces_begin, 
                              InputIterator2 faces_end,
                              Container& points)
  {
    std::set<Point,less_xy<Point> >  endpoints(points_begin, points_end);
    
    for (InputIterator2 face_iter=faces_begin; face_iter!=faces_end; ++face_iter){
      if ((*face_iter)->does_outer_ccb_exist()){
        Arrangement::Ccb_halfedge_circulator circ=(*face_iter)->outer_ccb();
        do {
          // inserting the point only if it's an end point.
          if (endpoints.find(circ->source()->point()) != endpoints.end())
             points.insert(circ->source()->point());
        } while (++circ != (*face_iter)->outer_ccb());  
      }
      
      Arrangement::Holes_iterator hit = (*face_iter)->holes_begin();
      for ( ;hit!= (*face_iter)->holes_end(); ++hit) {
        Arrangement::Ccb_halfedge_circulator circ = (*hit);
        do {
          if (endpoints.find(circ->source()->point()) != endpoints.end())
             points.insert(circ->source()->point());
        } while (++circ != *hit);
      }
    }
  }
  
  // This function is written for the symmetric difference.
  // the two former functions: get_points_below_vertex_with_covering and 
  // get_points_below_edge_with_covering could have return a point inside a polygon 
  // since this check was not perfromed. Here we clean such points.
  template <class InputIterator, class Container>
  void  remove_points_below_faces(InputIterator polygons_begin,
                                  InputIterator polygons_end, 
                                  Container& points)
  {
    for (typename Container::iterator p_iter=points.begin(); 
         p_iter!=points.end(); ){
      InputIterator poly_iter;
      for (poly_iter=polygons_begin; poly_iter!=polygons_end; ++poly_iter){
        if (poly_iter->has_on_bounded_side(*p_iter)){
          points.erase(p_iter++);
          break;
        }
      }
      if (poly_iter==polygons_end)
        ++p_iter;
    }
  }
  
  template <class InputIterator1, class InputIterator2,class OutputIterator>
  void get_vertices_with_covering(InputIterator1 polygons_begin, 
                                  InputIterator1 polygons_end,
                                  InputIterator2 faces_begin,
                                  InputIterator2 faces_end,
                                  OutputIterator vertices, 
                                  const Arrangement&  arr,
                                  int covering)
  {
    std::set<Point,less_xy<Point> >  endpoints;
    std::set<Point,less_xy<Point> >  covering_endpoints; 
    //std::set<Point,less_xy<Point> >  common_endpoints; 
    // holds equal endpoints from different polygons.
    
    std::list<Point> intersections;
    
    // inserting to common points the endpoints which have covering number of vertices above.
    for (InputIterator1 poly_iter=polygons_begin; 
         poly_iter != polygons_end; ++poly_iter){
      for (Polygon::Vertex_iterator v_iter=poly_iter->vertices_begin(); 
           v_iter!=poly_iter->vertices_end(); ++v_iter)
        endpoints.insert(*v_iter);
    }
    
    get_points_below_vertex_with_covering(endpoints.begin(),endpoints.end(),
                                          polygons_begin,polygons_end,
                                          covering_endpoints,covering);
    
    get_points_below_edge_with_covering(endpoints.begin(),endpoints.end(),
                                        polygons_begin,polygons_end,
                                        covering_endpoints,covering);
    
    for (Arrangement::Vertex_const_iterator v_iter=arr.vertices_begin();
         v_iter!=arr.vertices_end(); ++v_iter){
      if (endpoints.find(v_iter->point()) == endpoints.end()) // intersection
        intersections.push_back(v_iter->point());
    }
    
    if (covering > 1){ // intersection or union.
      get_points_below_faces(endpoints.begin(),endpoints.end(),
                             faces_begin, faces_end,
                             covering_endpoints);

      std::copy(covering_endpoints.begin(),covering_endpoints.end(),vertices);
      std::copy(intersections.begin(),intersections.end(),vertices);
    }
    
    else if (covering == 1)  { // symmetric difference.
      remove_points_below_faces(polygons_begin,polygons_end,covering_endpoints);
      std::copy(covering_endpoints.begin(),covering_endpoints.end(),vertices);
    }
  }  
  
  // The interface for an intersection function
  template <class OutputIterator>
  void intersect_polygons(Arrangement arr,
                          OutputIterator faces)
  {
    //insert_polygons(arr, in_poly_list);
      
    // faces with a covering two are faces that are in the intersection
    // of the two polygons.
    //get_faces_with_covering(arr, faces, in_poly_list.size());
    
    //polygons_from_faces(arr, face_it_list, out_poly_list);
    
    //if (begin == faces) return 0; else return 1;
  }

  //----------------------------------- compating Planar_map features
  class Equal_Halfedge {
  public:
    Equal_Halfedge(const Arrangement::Halfedge& h) : h_(h) {}
    
    bool operator()(const Halfedge& halfedge) const
    {

      Traits  traits;
      
      /* // for debugging!
        cout<<"("<<CGAL::to_double(h_.curve().source().x())<<","<< CGAL::to_double(h_.curve().source().y())<<") ";
        cout<<"("<<CGAL::to_double(h_.curve().target().x())<<","<< CGAL::to_double(h_.curve().target().y())<<") "<<endl;
 
        cout<<"("<<CGAL::to_double(halfedge.curve().source().x())<<","<< CGAL::to_double(halfedge.curve().source().y())<<") ";
        cout<<"("<<CGAL::to_double(halfedge.curve().target().x())<<","<< CGAL::to_double(halfedge.curve().target().y())<<") "<<endl;
        cout<<endl;
        // end debugging. */

      return ((h_.curve() == halfedge.curve() || 
               traits.curve_flip(h_.curve()) ==  halfedge.curve()) &&
              h_.source()->point() == halfedge.source()->point() &&
              h_.target()->point() == halfedge.target()->point());
    }
    
  private:
    const Arrangement::Halfedge& h_;
  };
  
  class Equal_Ccb {
  public:
    Equal_Ccb(Arrangement::Ccb_halfedge_const_circulator outer_ccb1_) : 
      outer_ccb1(outer_ccb1_) {}
    
    bool operator()(Ccb_halfedge_const_circulator outer_ccb2) const
    {  
      Arrangement::Ccb_halfedge_const_circulator circ1=outer_ccb1;
      Ccb_halfedge_const_circulator circ2=outer_ccb2;
      
      unsigned int size1=0,size2=0;
      do {
        // for debugging!
        /*cout<<"("<<CGAL::to_double(circ1->curve().source().x())<<","<< 
          CGAL::to_double(circ1->curve().source().y())<<") ";
        cout<<"("<<CGAL::to_double(circ1->curve().target().x())<<","<< 
        CGAL::to_double(circ1->curve().target().y())<<") "<<endl;*/
        
        ++size1;
      } while (++circ1 != outer_ccb1);
      
      do {
        // for debugging!
        /*cout<<"("<<CGAL::to_double(circ2->curve().source().x())<<","<< 
          CGAL::to_double(circ2->curve().source().y())<<") ";
        cout<<"("<<CGAL::to_double(circ2->curve().target().x())<<","<< 
        CGAL::to_double(circ2->curve().target().y())<<") "<<endl;*/
      
        ++size2;
      } while (++circ2 != outer_ccb2);
      
      //CGAL_For_all(circ1, outer_ccb1) { ++size1; }   // does not compile.
      //CGAL_For_all(circ2, outer_ccb2) { ++size2; }
      
      if (size1 != size2)
        return false;
      
      //if (CGAL::circulator_size(outer_ccb1) !=   // Dees not compile!
      //    CGAL::circulator_size(outer_ccb2) )
      //  return false;
        
      Arrangement::Ccb_halfedge_const_circulator ccb1 = outer_ccb1;
      Ccb_halfedge_const_circulator ccb2 = outer_ccb2;
      
      do { 
        if ( Equal_Halfedge(*ccb1)(*ccb2) )
          break;
      } while (++ccb1 != outer_ccb1);
      
      if ( !Equal_Halfedge(*ccb1)(*ccb2) )
        return false;  // did not find the first halfedge.
      
      // now ccb1 plays as the new starting point of outer_ccb1.
      do {
        if (! Equal_Halfedge(*ccb1)(*ccb2) )
          return false;      
        ++ccb1;
      } while (++ccb2 != outer_ccb2); 
        
      return true;
    }
    
  private:
    Arrangement::Ccb_halfedge_const_circulator  outer_ccb1;
  };
  
  class Equal_Face {
  public:
    Equal_Face(const Arrangement::Face& f1_) 
      : f1(f1_) {}
    
    bool operator()(const Face& f2) const
    { 
      if (std::distance(f1.holes_begin(),f1.holes_end()) !=
          std::distance(f2.holes_begin(),f2.holes_end()) )
        return false;
      
      for (Arrangement::Holes_const_iterator holes_iter1 = 
             f1.holes_begin();
           holes_iter1 != f1.holes_end(); ++holes_iter1) {
        
        Equal_Ccb  pred(*holes_iter1);
        Holes_const_iterator iter = std::find_if(f2.holes_begin(),
                                                 f2.holes_end(), 
                                                 pred);
        
        if (iter == f2.holes_end())
          return false;
      }
      
      if (f1.is_unbounded()) {
        if (f2.is_unbounded())
          return true;
        else 
          return false;
      }
      
      if (f2.is_unbounded()) {
        if (f1.is_unbounded())
          return true;
        else 
          return false;
      }
      
      // f1 and f2 are bounded.
      
      Equal_Ccb  pred(f1.outer_ccb());
      
      return ( pred(f2.outer_ccb()) );
    }
  private:
    const Arrangement::Face& f1;
  };
  
  template <class InputIterator1, class InputIterator2>
  void compare_vertices(InputIterator1 begin1, InputIterator1 end1, 
                        InputIterator2 begin2, InputIterator2 end2)
  {
    // for debugging:
    cout<<"vertices sizes:"<<std::distance(begin1,end1)<<" "
        << std::distance(begin2,end2) <<endl;
    
    //for (InputIterator1 iter1=begin1; iter1!=end1; ++iter1)
    //cout << CGAL::to_double(iter1->x()) << " "<<
    //         CGAL::to_double(iter1->y()) << endl;
    
    //for (InputIterator2 iter2=begin2; iter2!=end2; ++iter2)
    // cout << CGAL::to_double(iter2->x()) << " "<<
    //         CGAL::to_double(iter2->y()) << endl;
    // end degugging.
   
    CGAL_assertion(std::distance(begin1,end1) == 
                   std::distance(begin2,end2));
    
    for (InputIterator1 iter1=begin1; iter1!= end1; ++iter1)
      CGAL_assertion(std::find(begin2,end2,*iter1) != end2);
      
  }
   
  template <class InputIterator1, class InputIterator2>
  void compare_faces(InputIterator1 begin1, InputIterator1 end1, 
                     InputIterator2 begin2, InputIterator2 end2)
  { 
    std::list<Arrangement::Face> faces1;
    std::list<Face> faces2;

    for (InputIterator1 iter1 = begin1; iter1 != end1; ++iter1){
      faces1.push_back(*(*iter1));
    }
    
    for (InputIterator2 iter2 = begin2; iter2 != end2; ++iter2)
      faces2.push_back(*(*iter2));
    
    // asserting both contianers have the same size of faces.
    CGAL_assertion(faces1.size() == faces2.size());
    
    for (std::list<Arrangement::Face>::iterator face_iter = faces1.begin();
         face_iter != faces1.end(); ++face_iter)
      {
        Equal_Face  pred(*face_iter);
        std::list<Face>::iterator f = std::find_if(faces2.begin(), 
                                                   faces2.end(),
                                                   pred);
        
        // asserting every face of faces1 is found in faces2.
        CGAL_assertion(f != faces2.end());
        
        faces2.erase(f);
      }
          
    // asserting we have covered all faces from faces2.
    CGAL_assertion(faces2.empty());
  }

  void  check_intersection(const Planar_map& pm1,
                           const Planar_map& pm2,
                           const Polygon& polygon1, 
                           const Polygon& polygon2)
  { 
    Polygon_list  polygons;
    polygons.push_back(polygon1);
    polygons.push_back(polygon2);
    
    Bops bop(pm1, pm2);
      
    Faces_container     bops_faces;
    Halfedges_container bops_halfedges;
    Vertices_container  bops_vertices;
    
    bop.intersection(bops_faces,bops_halfedges,bops_vertices);

    std::list<Arrangement::Face_iterator> faces;
    Arrangement arr;
    insert_polygons(arr,polygons);
    get_faces_with_covering(std::back_inserter(faces),arr,polygons.size());
    
    cout<<"Face sizes:"<<faces.size()<<" "<<bops_faces.size()<<endl;
    
    compare_faces(faces.begin(), faces.end(), 
                  bops_faces.begin(), bops_faces.end());

    std::cout << "Intersection: Faces --- O.K "<< std::endl;

    std::list<Point>   points;
    
    get_vertices_with_covering(polygons.begin(), 
                               polygons.end(),
                               faces.begin(), 
                               faces.end(),
                               std::back_inserter(points),
                               arr,
                               polygons.size());
 
    std::set<Point, less_xy<Point> > set_points(points.begin(), points.end());
    std::list<Point>  bops_points;
    for (Vertices_container::iterator iter=bops_vertices.begin();
         iter!=bops_vertices.end(); ++iter)
      bops_points.push_back((*iter)->point());
    
    compare_vertices(set_points.begin(), set_points.end(), 
                     bops_points.begin(), bops_points.end());
    
    std::cout << "Intersection: Vertices --- O.K "<< std::endl;
    //compare_halfedges(halfedges.begin(), halfedges.end(), 
    //              bops_halfedges.begin(), bops_halfedges.end());
  }

  void  check_union(const Planar_map& pm1,
                    const Planar_map& pm2,   //const Bops& bop, 
                    const Polygon& polygon1, 
                    const Polygon& polygon2)
  {
    Polygon_list  polygons;
    polygons.push_back(polygon1);
    polygons.push_back(polygon2);
    
    Bops bop(pm1, pm2);
     
    Faces_container     bops_faces;
    Halfedges_container bops_halfedges;
    Vertices_container  bops_vertices;
    
    bop.Union(bops_faces,bops_halfedges,bops_vertices);

    std::list<Arrangement::Face_iterator> union_faces[2];
    // faces[0] will contain the symmetric diff faces, and faces[1] will contain the intersection faces.
    Arrangement arr;
    insert_polygons(arr,polygons);

    // For union we have to define all covers which are not 0.
    unsigned int i=1;
    for ( ; i <= polygons.size(); ++i) 
      get_faces_with_covering(std::back_inserter(union_faces[i-1]),arr,i);
    
    std::list<Arrangement::Face_iterator> faces(union_faces[0].begin(), union_faces[0].end());
    std::copy(union_faces[1].begin(), union_faces[1].end(), std::back_inserter(faces));
    
    cout<<"Faces sizes:"<<faces.size()<<" "<<bops_faces.size()<<endl;
    
    compare_faces(faces.begin(), faces.end(), 
                  bops_faces.begin(), bops_faces.end());
    
    std::cout << "Union: Faces --- O.K "<< std::endl;

    std::list<Point>   points;
    //std::set<Point,less_xy<Point> >  union_points; 
    for (i=1; i <= polygons.size(); ++i) {
      get_vertices_with_covering(polygons.begin(), 
                                 polygons.end(),
                                 union_faces[i-1].begin(), 
                                 union_faces[i-1].end(),
                                 std::back_inserter(points),
                                 arr,
                                 i);

      //for ( std::list<Point>::iterator p_iter
    }
    
    std::set<Point, less_xy<Point> > set_points(points.begin(), points.end());
    
    std::list<Point>  bops_points;
    for (Vertices_container::iterator iter=bops_vertices.begin();
         iter!=bops_vertices.end(); ++iter)
      bops_points.push_back((*iter)->point());
    
    compare_vertices(set_points.begin(), set_points.end(), 
                     bops_points.begin(), bops_points.end());
    
    std::cout << "Union: Vertices --- O.K "<< std::endl;
    //compare_halfedges(halfedges.begin(), halfedges.end(), 
    //              bops_halfedges.begin(), bops_halfedges.end());

  }
 
  void  check_symmetric_difference(const Planar_map& pm1,
                                   const Planar_map& pm2,
                                   const Polygon& polygon1, 
                                   const Polygon& polygon2)
  {
    Polygon_list  polygons;
    polygons.push_back(polygon1);
    polygons.push_back(polygon2);

    Bops bop(pm1, pm2);
    
    Faces_container     bops_faces;
    Halfedges_container bops_halfedges;
    Vertices_container  bops_vertices;

    bop.symmetric_difference(bops_faces,bops_halfedges,bops_vertices);
    
    std::list<Arrangement::Face_iterator> faces;
    Arrangement arr;
    insert_polygons(arr,polygons);

    // For symmetric difference we have to define all covers to be 1.
    get_faces_with_covering(std::back_inserter(faces), arr, 1);

    cout<<"Face sizes:"<<faces.size()<<" "<<bops_faces.size()<<endl;
    
    compare_faces(faces.begin(), faces.end(), 
                  bops_faces.begin(), bops_faces.end());
    
    std::cout << "Symmetric Difference: Faces --- O.K "<< std::endl;

    std::list<Point>   points;
    
    get_vertices_with_covering(polygons.begin(), 
                               polygons.end(),
                               faces.begin(), 
                               faces.end(),
                               std::back_inserter(points),
                               arr,
                               1);
    
    std::set<Point, less_xy<Point> > set_points(points.begin(), points.end());
    std::list<Point>  bops_points;
    for (Vertices_container::iterator iter=bops_vertices.begin();
         iter!=bops_vertices.end(); ++iter)
      bops_points.push_back((*iter)->point());
    
    compare_vertices(set_points.begin(), set_points.end(), 
                     bops_points.begin(), bops_points.end());
    
    std::cout << "Symmetric Difference: Vertices --- O.K "<< std::endl;
    //compare_halfedges(halfedges.begin(), halfedges.end(), 
    //              bops_halfedges.begin(), bops_halfedges.end());

  }
  
//------------------------------- Scanning input functnions  
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

  void read_file_build_creator(std::ifstream& file, Planar_map& pm)
  { 
    // read number of curves
    unsigned int num_curves = get_next_int(file);
    
    cout<<"num_curves="<<num_curves<<endl;
   
    
    std::list<Curve> curves;
    //Curve curr_curve;
    // read curves (test specific)
    int  x1,y1,x2,y2;
    while (num_curves--) {
      //file >> curr_curve;
      file >>x1 >>y1 >> x2>> y2;
      curves.push_back(Curve(Point(x1,y1),Point(x2,y2)));
    }
    
    Traits traits;
    sweep_to_construct_planar_map_2(curves.begin(),curves.end(), traits, pm);
  }
  
public:

  void start(char * filename)
  {
    std::ifstream file(filename);
    
    PmWalkPL     pl_walk1, pl_walk2;
    Planar_map   pm1(&pl_walk1), pm2(&pl_walk2);
    
    read_file_build_creator(file, pm1);
    read_file_build_creator(file, pm2);

    pm1.unbounded_face()->set_ignore_bop(false); 
    pm2.unbounded_face()->set_ignore_bop(false);

    Polygon polygon1, polygon2;
    
    for (Vertex_iterator v_iter1=pm1.vertices_begin();
         v_iter1!=pm1.vertices_end(); ++v_iter1)
      polygon1.push_back(v_iter1->point());

    for (Vertex_iterator v_iter2=pm2.vertices_begin();
         v_iter2!=pm2.vertices_end(); ++v_iter2)
      polygon2.push_back(v_iter2->point());
    
    Faces_container     bops_faces;
    Halfedges_container bops_halfedges;
    Vertices_container  bops_vertices;
    
    //bop.intersection(bops_faces,bops_halfedges,bops_vertices); 
    //cout<<"intersetion: bops_faces.size() "<<bops_faces.size()<<endl;

    check_intersection(pm1,pm2,polygon1,polygon2);
    
    //bop.Union(bops_faces,bops_halfedges,bops_vertices);
    //cout<<"union: bops_faces.size() "<<bops_faces.size()<<endl;
    
    check_union(pm1,pm2,polygon1,polygon2);

    check_symmetric_difference(pm1,pm2,polygon1,polygon2);
  }
};

int main(int argc, char* argv[])
{
  
  Boolean_operations_test  test;
  
  if (argc < 1 || argc > 2) {
    std::cout << "usage: test data_file" << std::endl;
    exit(1);
  }

  test.start(argv[1]);
  return 0;
}

#endif // CGAL_ARR_TEST_LEDA_CONFLICT

