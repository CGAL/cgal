// Copyright (c) 2016  INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id:  $
//
//
// Author(s)     : Iordan Iordanov

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H

#include <vector>
#include <map>

#include <CGAL/Square_root_2_field.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_ds_vertex_base_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_ds_face_base_2.h>
#include <CGAL/Hyperbolic_octagon_word_4.h>
//#include <CGAL/Hyperbolic_octagon_word_8.h>


namespace CGAL {

  template <  class GT,
  class TDS = Triangulation_data_structure_2<
  Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<GT>,
  Periodic_4_hyperbolic_triangulation_ds_face_base_2<GT> 
  >
  >
  class Periodic_4_hyperbolic_Delaunay_triangulation_2 :
  public Periodic_4_hyperbolic_triangulation_2<GT, TDS> {

    typedef Periodic_4_hyperbolic_Delaunay_triangulation_2<GT, TDS>   Self;
    typedef Periodic_4_hyperbolic_triangulation_2<GT, TDS>            Base;

  public:

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2  
    using Base::cw;
    using Base::ccw;
    using Base::geom_traits;
#endif

    typedef typename Base::Locate_type                         Locate_type;
    typedef typename Base::Geometric_traits                    Geometric_traits;
    typedef typename Base::Triangulation_data_structure        Triangulation_data_structure;
    //typedef typename Base::Int                                 Int;
    typedef typename Base::Offset                              Offset;
    typedef typename Base::Circle_2                            Circle_2;
    typedef Circle_2                                           Circle;
    typedef typename Base::Point_2                             Point_2;
    typedef Point_2                                            Point;
    typedef typename Base::Segment_2                           Segment_2;
    typedef Segment_2                                          Segment;
    typedef typename Base::Triangle_2                          Triangle_2;
    typedef Triangle_2                                         Triangle;

    typedef typename Base::Periodic_point                      Periodic_point;
    typedef typename Base::Periodic_segment                    Periodic_segment;
    typedef typename Base::Periodic_triangle                   Periodic_triangle;  

    typedef typename Base::Vertex                              Vertex;
    typedef typename Base::Edge                                Edge;
    typedef typename Base::Face                                Face;

    typedef typename Base::Vertex_handle                       Vertex_handle;
    typedef typename Base::Face_handle                         Face_handle;

    typedef typename Base::size_type                           size_type;
    typedef typename Base::difference_type                     difference_type;

    typedef typename Base::Face_iterator                       Face_iterator;
    typedef typename Base::Edge_iterator                       Edge_iterator;
    typedef typename Base::Vertex_iterator                     Vertex_iterator;
    typedef typename Base::Face_circulator                     Face_circulator;
    typedef typename Base::Edge_circulator                     Edge_circulator;
    typedef typename Base::Vertex_circulator                   Vertex_circulator;
    typedef typename Base::Line_face_circulator                Line_face_circulator;

    typedef Face_iterator                                      All_faces_iterator;
    typedef Edge_iterator                                      All_edges_iterator;
    typedef Vertex_iterator                                    All_vertices_iterator;

    typedef Face_iterator                                      Finite_faces_iterator;
    typedef Edge_iterator                                      Finite_edges_iterator;
    typedef Vertex_iterator                                    Finite_vertices_iterator;

    typedef Periodic_4_hyperbolic_triangulation_unique_vertex_iterator_2<Self>
    Unique_vertex_iterator;

  private:
    typedef typename GT::FT                                  FT;
    typedef std::pair< Vertex_handle, Offset >               Virtual_vertex;
    typedef std::map<Vertex_handle, Virtual_vertex>          Virtual_vertex_map;
    typedef typename Virtual_vertex_map::const_iterator      Virtual_vertex_map_it;
    typedef std::map<Vertex_handle, std::vector<Vertex_handle > > 
    Virtual_vertex_reverse_map;
    typedef typename Virtual_vertex_reverse_map::const_iterator
    Virtual_vertex_reverse_map_it;
    typedef Triple< Vertex_handle, Vertex_handle, Vertex_handle >
    Vertex_triple;

  public:
    typedef Periodic_4_hyperbolic_triangulation_triangle_iterator_2<Self>
    Periodic_triangle_iterator;
    typedef Periodic_4_hyperbolic_triangulation_segment_iterator_2<Self>
    Periodic_segment_iterator;
    typedef Periodic_4_hyperbolic_triangulation_point_iterator_2<Self>
    Periodic_point_iterator;

    typedef Point                                            value_type;
    typedef const value_type&                                const_reference;
    typedef Tag_false                                        Weighted_tag;

  //protected:
    //Geometric_traits              _gt;
    //Triangulation_data_structure  _tds; 
    //Circle                        _domain;



  protected:
    mutable std::vector<Vertex_handle> v_offsets;

  public:

    Periodic_4_hyperbolic_Delaunay_triangulation_2(Geometric_traits gt) : 
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>(gt) { }  

    Periodic_4_hyperbolic_Delaunay_triangulation_2(
      const Circle_2 domain = Circle_2(Point_2(FT(0),FT(0)), FT(1*1)), 
      const Geometric_traits &gt = Geometric_traits() ) :
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>(domain, gt) { }

    Periodic_4_hyperbolic_Delaunay_triangulation_2(const Periodic_4_hyperbolic_Delaunay_triangulation_2& tr) :
    Periodic_4_hyperbolic_triangulation_2<GT, TDS>(tr) { }

    Vertex_handle insert(const Point  &p,
     Face_handle start = Face_handle() );

    Vertex_handle insert(const Point& p,
     Locate_type lt,
     Face_handle loc, int li );

  template < class InputIterator >
    std::ptrdiff_t insert(InputIterator first, 
      InputIterator last, 
      bool is_large_point_set = true);

    bool locally_Delaunay(const Face_handle&, int);
    void propagating_flip(Face_handle&, int);
    void restore_Delaunay(Vertex_handle);
    void flip_single_edge(Face_handle, int);
    bool flippable(Face_handle, int);

  protected:
    class Conflict_tester;




    bool _side_of_octagon( const Face_handle& fh, const Offset& offset) const {
      int cnt = 0;
      typename GT::Side_of_fundamental_octagon side;
      for (int j = 0; j < 3; j++) {
        Offset o = offset.inverse().append(fh->vertex(j)->get_offset());
        Point  p = o.apply( fh->vertex(j)->point() );
        if ( side(p) == CGAL::ON_UNBOUNDED_SIDE ) {
          if ( p.y() + (CGAL_PI / FT(8))*p.x() > 0 ) {
            cnt++;
          } else {
          }
        }
      }
      return (cnt == 0);
    }

};  // class Periodic_4_hyperbolic_Delaunay_triangulation_2




template<class Gt, class Tds>
bool Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::flippable(Face_handle f, int i)
{
    /*
  Face_handle nb = f->neighbor(i);
  int j = nb->index(f);

  const Point *p[4];

  p[0] = &f->vertex(i)->point();      // i
  p[1] = &nb->vertex(j)->point();     // opposite
  p[2] = &f->vertex(ccw(i))->point(); // ccw
  p[3] = &f->vertex(cw(i))->point();  // cw

  if (f->has_zero_offsets() && nb->has_zero_offsets()) {
    if (orientation(*p[0], *p[1], *p[2]) == LEFT_TURN)
      return false;
    if (orientation(*p[0], *p[1], *p[3]) == RIGHT_TURN)
      return false;
  } else {
    Offset off[4];
    off[0] = f->offset(i);
    off[1] = f->neighbor_offset(j).append(nb->offset(j)); 
    off[2] = f->offset(ccw(i));
    off[3] = f->offset(cw(i));

    if (orientation(*p[0], *p[1], *p[2], off[0], off[1], off[2]) == LEFT_TURN)
      return false;
    if (orientation(*p[0], *p[1], *p[3], off[0], off[1], off[3]) == RIGHT_TURN)
      return false;
  }
    */
  
  return true;
}


template<class Gt, class Tds>
void Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
flip_single_edge(Face_handle f, int i)
{
  CGAL_triangulation_precondition(f != Face_handle());
  CGAL_triangulation_precondition(i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition(this->dimension() == 2);

  CGAL_triangulation_precondition(flippable(f, i));

  Face_handle nb = f->neighbor(i);
  int j = nb->opposite_index(f);

/*
  cout << endl;
  cout << "Flipping edge between faces " << f->get_number() << " and " << nb->get_number() << endl;
  cout << "Offsets of f: " << f->offset(0) << ", " << f->offset(1) << ", " << f->offset(2);
  cout << ", offsets of nb: " << nb->offset(0) << ", " << nb->offset(1) << ", " << nb->offset(2) << endl;
  cout << "Neighbors before: " << f->neighbor(cw(i))->get_number() << ", " << nb->neighbor(ccw(j))->get_number();
  cout << ", " << nb->neighbor(cw(j))->get_number() << ", " << f->neighbor(ccw(i))->get_number() << endl;
  */
  f->store_offsets();
  nb->store_offsets(f->neighbor_offset(i));

  this->tds().flip(f, i);

  i = ccw(i);
  nb = f->neighbor(i);
  j = nb->opposite_index(f);

  f->restore_offsets();
  nb->restore_offsets();
/*
  cout << "Now f is " << f->get_number() << " and nb is " << nb->get_number() << endl;
  cout << "Offsets of f: " << f->offset(0) << ", " << f->offset(1) << ", " << f->offset(2);
  cout << ", offsets of nb: " << nb->offset(0) << ", " << nb->offset(1) << ", " << nb->offset(2) << endl;
  cout << "Neighbors after: " << f->neighbor(cw(i))->get_number() << ", " << nb->neighbor(ccw(j))->get_number();
  cout << ", " << nb->neighbor(cw(j))->get_number() << ", " << f->neighbor(ccw(i))->get_number() << endl;
  cout << endl;
*/
  
  for (int j = 0; j < 3; j++) {
  	f->vertex(j)->remove_offset();
  	nb->vertex(j)->remove_offset();
  }


/*
  if (f->has_zero_offsets() && nb->has_zero_offsets()) {
    _tds.flip(f, i);
    return;
  } else {
      int j = nb->opposite_index(f);
      
      Offset nofi = f->neighbor_offset(i);
      Offset nofj = nb->neighbor_offset(j);

      cout << "nofi = " << nofi << ", nofj = " << nofj << endl;

      Vertex_handle vh[] = { f->vertex(i),
                             f->vertex(ccw(i)),
                             nb->vertex(j),
                             f->vertex(cw(i))    };
      
      Offset o[]         = { f->offset(i),
                             f->offset(ccw(i)),
                             f->neighbor_offset(i).append(nb->offset(j)), 
                             f->offset(cw(i))    };

      _tds.flip(f, i);
      
      i = ccw(i);
      nb = f->neighbor(i);
      j = nb->opposite_index(f);

      f->set_offsets();
      nb->set_offsets();

      f->set_offset(i, o[1]);      
      f->set_offset(ccw(i), o[2]);      
      f->set_offset(cw(i), o[0]); 

      cout << "------- face " << f->get_number() << " -------" << endl;
      for (int c = 0; c < 3; c++) {
        cout << "  neighbor " << c << " is face " << f->neighbor(c)->get_number() << " with n-offset " << f->neighbor_offset(c) << endl; 
      }
      cout << endl;
	
    

      nb->set_offset(j, o[3]);      
      nb->set_offset(cw(j), o[2]);
      nb->set_offset(ccw(j), o[0]);
      cout << "------- face " << nb->get_number() << " -------" << endl;
      for (int c = 0; c < 3; c++) {
        cout << "  neighbor " << c << " is face " << nb->neighbor(c)->get_number() << " with n-offset " << nb->neighbor_offset(c) << endl; 
      }
      cout << endl;
  } */


    }



template < class Gt, class Tds >
    void
    Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
    restore_Delaunay(Vertex_handle v)
    {
      Face_handle f = v->face();
      Face_handle next;
      int i;
      Face_handle start(f);
      do {
        i = f->index(v);
      next = f->neighbor(ccw(i));  // turn ccw around v
      propagating_flip(f, i);
      f = next;
    } while(next != start);
  }


template < class Gt, class Tds >
  void
  Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
  propagating_flip(Face_handle& f, int i)
  {
    Face_handle nb = f->neighbor(i);

    if (locally_Delaunay(f, i)) {
      return;
    }

    this->flip_single_edge(f, i);
    propagating_flip(f, i);
    i = nb->index(f->vertex(i));
    propagating_flip(nb, i);
  }



template < class Gt, class Tds >
  bool
  Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
  locally_Delaunay(const Face_handle &f, int nbi)
  {
    CGAL_BRANCH_PROFILER("locally_Delaunay(), simplicity check failures", tmp);

    Face_handle nb = f->neighbor(nbi);

    bool simplicity_criterion = f->has_zero_offsets() && nb->has_zero_offsets();

    const Point *p[4];
    Offset off[4];
    for (int index = 0; index < 3; index++)
    {
      p[index]   = &f->vertex(index)->point();
      off[index] = f->offset(index);
    }

    p[3]   = &nb->opposite_vertex(f)->point();
    off[3] = f->neighbor_offset(nbi).append(nb->opposite_offset(f)); //bu_append(f->neighbor_offset(nbi), nb->opposite_offset(f));

    Oriented_side os;
    os = this->side_of_oriented_circle( *p[0],  *p[1],  *p[2],  *p[3],
                                        off[0], off[1], off[2], off[3]);

    return (ON_POSITIVE_SIDE != os);
  }



#define INSERT_WITH_FLIPS 0


template < class Gt, class Tds >
  inline
  typename Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
  Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
  insert(const Point  &p,  Face_handle start) {

    Vertex_handle v;

    typedef typename Gt::Side_of_fundamental_octagon Side_of_fundamental_octagon;

    Side_of_fundamental_octagon check = Side_of_fundamental_octagon();
    CGAL::Bounded_side side = check(p);

    if (side != CGAL::ON_UNBOUNDED_SIDE) {
      Offset loff;
      if ( start == Face_handle() ) {
        Locate_type lt;
        int li;
        start = this->euclidean_visibility_locate(p, lt, li, loff);
        //cout << "Query point is located in face " << start->get_number() << " with loff = " << loff << " and LOCATE_TYPE = " << lt << endl;
        if (lt == Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::VERTEX) {
          return Vertex_handle();
        }
      }

      if (INSERT_WITH_FLIPS == 1) {
        cout << endl << endl << "======= LAST INSERTION STARTS HERE =======" << endl << endl;

        cout << "Point inserted in face " << start->get_number() << ", vertices: " << endl;
        for (int i = 0; i < 3; i++) { 
          cout << start->vertex(i)->idx() << " with offset " << start->offset(i) << ", " << endl;
        }
        if (!loff.is_identity()) {
         start->apply(loff);
         cout << " --- face inverted --- " << endl;
       }


       int next_number = this->number_of_faces();
       Vertex_handle ov[] = {  start->vertex(0),
        start->vertex(1),
        start->vertex(2) };

    // Save the offsets in the vertices
        start->store_offsets();

        Vertex_handle new_vertex = this->insert_in_face(p, start);
        new_vertex->set_point( p );    	
        new_vertex->set_idx( this->number_of_vertices() );

        cout << "New vertex has id " << new_vertex->idx() << endl;

    // Assign offsets
        Face_circulator iface = this->tds().incident_faces(new_vertex), end(iface);
        if (iface != 0) {
          do {
        // Restore offsets from the vertices
            iface->restore_offsets();

            if (iface->get_number() == -1) {
              iface->set_number(next_number++);
            }
          } while (++iface != end);

          cout << "New faces created: " << endl << "------------------------" << endl;
          do {
            cout << "Face " << iface->get_number() << ", vertices: ";
            for (int i = 0; i < 3; i++) {
              cout << iface->vertex(i)->idx() << " with offset " << iface->offset(i) << ", ";
            } 
            cout << endl;
            cout << "Neighbor offsets: " << iface->neighbor_offset(0) << ", " << iface->neighbor_offset(1) << ", " << iface->neighbor_offset(2) << endl;
          } while (++iface != end);
          cout << "--------- end ----------" << endl << endl;
        }

    // Delete offsets from the vertices -- they are kept in the faces now
        for (int i = 0; i < 3; i++) {
         ov[i]->remove_offset();
       }

       restore_Delaunay(new_vertex);

       CGAL_triangulation_assertion(this->is_valid());

       return new_vertex;
     

     } else { // means do not insert_with_flips, i.e., insert with star
      
      // Ok, we really insert the point now.
      // First, find the conflict region.
      std::vector<Face_handle> faces;
      
      this->find_conflicts(start, p, loff, std::back_inserter(faces));

         //CGAL_assertion(faces.size() > 0);
      Side_of_fundamental_octagon side;

      // Saving the offsets into the vertices (temporarily)
      for (int i = 0; i < faces.size(); i++) {
        faces[i]->store_offsets();
      }

      v = this->insert_in_hole(p, faces.begin(), faces.end());

      //cout << "Done inserting in the hole!" << endl;
      Face_circulator ifc = this->tds().incident_faces(v), done(ifc);
      do {
        ifc->restore_offsets(loff);
        ifc->tds_data().clear();
        ifc->make_canonical();
        //cout << "  Face " << ifc->get_number() << ": vertices " << ifc->vertex(0)->idx() << ", " << ifc->vertex(1)->idx() << ", " << ifc->vertex(2)->idx() << " with offsets " << ifc->offset(0) << ", " << ifc->offset(1) << ", " << ifc->offset(2) << ", " << endl;
      } while (++ifc != done);

      Vertex_circulator ivc = this->tds().incident_vertices(v), done_v(ivc);
      do {
      	ivc->remove_offset();
      } while (++ivc != done_v);

        CGAL_triangulation_assertion(this->is_valid());

        //cout << endl << "=============================" << endl << endl;

        cout << "Vertices in triangulation: " << this->number_of_vertices() << endl;

        return v;
      }
    }

    CGAL_triangulation_assertion(this->is_valid());

    return Vertex_handle();
  }



template < class GT, class Tds >
  class Periodic_4_hyperbolic_Delaunay_triangulation_2<GT,Tds>::Conflict_tester
  {
  // stores a pointer to the triangulation,
  // a point, and an offset
    const Self *t;
    Point p;

  public:
  /// Constructor
    Conflict_tester(const Self *_t) : t(_t), p(Point()) {}
    Conflict_tester(const Point &pt, const Self *_t) : t(_t), p(pt) { }

  /** The functor
    *
    * gives true if the circumcircle of c contains p
    */
    bool operator()(const Face_handle f, const Offset &off) const {
      return (t->_side_of_circle(f, p, off) == ON_BOUNDED_SIDE);
    }

    bool operator()(const Face_handle f, const Point& pt, const Offset &off) const {
      return (t->_side_of_circle(f, pt, off) == ON_BOUNDED_SIDE);
    }

    int compare_weight(Point, Point) const
    {
      return 0;
    }

    bool test_initial_face(Face_handle f, const Offset &off) const
    {
      if (!(operator()(f, off)))
        CGAL_triangulation_assertion(false);
      return true;
    }

    void set_point(const Point &_p) {
      p = _p;
    }

    const Point &point() const {
      return p;
    }

  };


template < class Gt, class Tds >
  inline
  typename Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
  Periodic_4_hyperbolic_Delaunay_triangulation_2<Gt, Tds>::
  insert(const Point  &p, Locate_type lt, Face_handle loc, int li)
  {
    std::cout << "Inseritng now! Right!" << std::endl;
    return Vertex_handle();
  }


} // namespace CGAL

#endif // CGAL_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
