#ifndef CGAL_USE_LEDA

// this demo uses LEDA
 
#else

#include "Polygon_intersection.h"

/////////////////////////////////////////////////////////////////////
// Intersecting Polygons
////////////////////////////////////////////////////////////////////

// we assume the polygons were inserted in CCW order!!

// The face_diff function defines if to increment or decrement
// the number of polygons when moving over an edge from one face to the other.
// Method: if the direction of circ is identical to the curve on it
// then we are inside a polygon (the polygons were inserted ccw) and therefore 
// we need to decrement, otherwise we are outside and need to increment

// The present example can be altered or generalized for other boolean
// operations.

// Commented out the regular face_diff that does not deal with overlaps
// as an educative example.
/*
int face_diff( Ccb_halfedge_circulator circ) {
  Traits t;
  if (circ->source()->point() == t.curve_source(circ->edge_node()->curve()) ) 
    return -1;     //we're inside, going outside
  else
    return 1;
}
*/

//generalized face_diff function, to acount for overlaps.
int face_diff (Ccb_halfedge_circulator circ) {
  Traits t;
  int diff = 0;
  Arr_2::Overlap_circulator oc = circ->overlap_edges();
  do {
    if (circ->source()->point() == t.curve_source(oc->x_curve()) ) 
    diff--;     //we're inside, going outside
  else
    diff++;
  } while (++oc != circ->overlap_edges());

  return diff;
} 

// covering_DFS will compute for each face in how many polygons it is.
// It is a recursive DFS function and will be called with the unbounded 
// face after its counter has been initialized to 0.
void covering_DFS(Face_handle f) {
  Ccb_halfedge_circulator start,circ;

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
  Arr_2::Holes_iterator hit = f->holes_begin();
  for (; Arr_2::Holes_iterator(hit)!=Arr_2::Holes_iterator(f->holes_end());
       ++hit)
  {
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
void insert_polygons(Arr_2 &arr, Polygon_list &in_poly_list)
{
  Polygon::Edge_const_iterator it;
  Arr_2::Curve_iterator        ci;
  Polygon_list::iterator       plit;

  // for each polygon in list
  for (plit = in_poly_list.begin(); plit != in_poly_list.end(); plit++) {

    // Make sure polygons are oriented counterclockwise
    // to satisfy assumption in DFS function
    if ( ! plit->is_counterclockwise_oriented())
      plit->reverse_orientation();

    // insert polygon to arrangement
    for (it = plit->edges_begin(); it != plit->edges_end(); it++)
      {
	ci = arr.insert(*it);
      }
  }
}

// Convert faces of the arrangement that are in the intersection
// to polygons.
void polygons_from_faces(Arr_2& arr,
			 std::list<Face_iterator>& face_it_list,
			 Polygon_list& poly_list)
{
  std::list<Face_iterator>::iterator  lit;
  //Arr_2::Ccb_halfedge_circulator cc;
  Polygon                        poly;
  
  for (lit = face_it_list.begin(); lit != face_it_list.end(); lit++) {

    poly.erase(poly.vertices_begin(), poly.vertices_end());
    Arr_2::Ccb_halfedge_circulator cc=(*lit)->outer_ccb();
    do {
      poly.push_back(cc->curve().source());
    } while (++cc != (*lit)->outer_ccb());
    poly_list.push_back(poly);
  }
}

// performs the extraction of data out of the processed arrangement
// if covering = 0, will perform union
// otherwise, if there are n polygons in the arrangement and covering == n
// then will perform intersection
void get_faces_with_covering(Arr_2& arr, 
			     std::list<Face_iterator>& unions, 
			     int covering)
{
  Face_handle uf = arr.unbounded_face();
  uf->counter = 0;
  covering_DFS(uf);
  
  //"collecting" the union boundary faces. 
  for(Face_iterator fit = arr.faces_begin(); fit!=arr.faces_end(); ++fit) {
    
    // if the face is covered by 'covering' 
    if (fit->counter == covering) {
      unions.push_back(fit);
    }
  }
}                                                                          

// The interface for an intersection function
bool intersect_polygons(Polygon_list &in_poly_list,
                        Polygon_list &out_poly_list)
{
  Arr_2               arr(new Arr_naive_pl());
  std::list<Face_iterator> face_it_list;

  insert_polygons(arr, in_poly_list);
  
  // faces with a covering two are faces that are in the intersection
  // of the two polygons.
  get_faces_with_covering(arr, face_it_list, in_poly_list.size());

  polygons_from_faces(arr, face_it_list, out_poly_list);

  if (out_poly_list.empty()) return 0; else return 1;
}

#endif
