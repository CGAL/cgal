//Here the input file will be a triangle. 
//Hence, all the insertion functions of Planar_map 
// will be invoked (in the order: insert_in_face_interior, 
//insert_from_vertex,insert_at_vertices). 
// Each insertion of one segment will invoke add_edge 
// of the notifier. The first will also 
// invoke add_hole, and the last call will invoke split_face.

#include <CGAL/Cartesian.h>

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Pm_segment_traits_2.h>

#include <CGAL/Planar_map_2.h>
#include <CGAL/Planar_map_2/Pm_change_notification.h>

#include <CGAL/assertions.h>

#include <CGAL/Polygon_2.h>

#include <iostream>
#include <set> 

CGAL_BEGIN_NAMESPACE         

template <class X_curve, class NT>
class Pm_halfedge_plus_length : public Pm_halfedge_base<X_curve>
{
public:
  typedef Pm_halfedge_base<X_curve>       halfedge_base;
  typedef Pm_halfedge_plus_length         halfedge_plus_length;
  typedef const halfedge_plus_length      const_halfedge_plus_length;
  typedef const_halfedge_plus_length*     const_pointer;
  typedef const_halfedge_plus_length&     const_ref;

  Pm_halfedge_plus_length() : halfedge_base() {}

  Pm_halfedge_plus_length(const X_curve& cv) : halfedge_base(cv) {}

  virtual ~Pm_halfedge_plus_length() {}
  
  const NT& length() const { return length_; }
  NT& length() { return length_; }
  
private:
  NT length_;
};

template <class NT>
class Pm_face_plus_area : public Pm_face_base 
{
public:
  typedef Pm_face_base              face_base;
  typedef Pm_face_plus_area         face_plus_area;
  typedef const face_plus_area      const_face_plus_area;
  typedef const_face_plus_area*     const_pointer;
  typedef const_face_plus_area&     const_ref;

private:
  struct less_face {
    bool operator() (const_pointer f1, const_pointer f2) const {
      return f1 < f2;
    }
  };

public:
  typedef std::set<const_pointer,less_face>::iterator   Inner_faces_iterator;
  typedef std::set<const_pointer,less_face>::const_iterator
                                                    Inner_faces_const_iterator;
  
  Pm_face_plus_area() : face_base() {}

  virtual ~Pm_face_plus_area() {}
  
  void insert_inner_face(const_pointer f)
    { inner_faces.insert(f); }
  
  Inner_faces_iterator inner_faces_begin() { return inner_faces.begin(); }
  Inner_faces_iterator inner_faces_end() { return inner_faces.end(); }

  Inner_faces_const_iterator inner_faces_begin() const
    { return inner_faces.begin(); }
  Inner_faces_const_iterator inner_faces_end() const
    { return inner_faces.end(); }

  const NT& full_area() const { return full_area_; }
  NT& full_area() { return full_area_; }
  
  const NT& area() const { return area_; }
  NT& area() { return area_; }
  
private:
  std::set<const_pointer, less_face>  inner_faces;
  NT full_area_;
  NT area_;
};

template <class Traits, class NT>
class Pm_measures_dcel : public Pm_dcel<Pm_vertex_base<typename Traits::Point>,
  Pm_halfedge_plus_length<typename Traits::X_curve, NT>,
  Pm_face_plus_area<NT> > 
{
public:  // CREATION
  
  Pm_measures_dcel() {}
  
};

template <class Planar_map_, class R>
class Notification : public Pm_change_notification<Planar_map_> {
public:
  typedef Planar_map_ Planar_map;
  typedef typename Planar_map::Traits Traits;
  
private:
  typedef CGAL::Polygon_traits_2<R>                     Polygon_traits;
  typedef std::list<typename Traits::Point>             Polygon_Container;
  typedef CGAL::Polygon_2<Polygon_traits, Polygon_Container> Polygon;
  typedef typename R::FT                                NT;
  
  /* struct less_face{
    bool operator()(typename Planar_map::Face_handle f1, 
                    typename Planar_map::Face_handle f2) const {
      return &*f1 < &*f2;
    }
    };*/

  NT  face_area(typename Planar_map::Face_handle face)
  {
    // std::set<typename Planar_map::Face_handle, less_face>  inner_faces;
    // inner_face is the set of all faces inside face (faces that adjacent to
    // its holes).
    for (typename Planar_map::Holes_iterator hit = face->holes_begin(); 
         hit != face->holes_end(); ++hit) {
      typename Planar_map::Ccb_halfedge_circulator cc(*hit);
      do{
        if (cc->twin()->face() != face){
          face->insert_inner_face(&*(cc->twin()->face()));
          
          // if (cc->twin()->face()->full_area() == 0)
          //   cc->twin()->face() = calc_full_area(cc->twin()->face());
        }   
      } while (++cc != *hit);
    }
    
    NT  inner_faces_area=0;
    for (typename Pm_face_plus_area<NT>::Inner_faces_iterator iter =
           face->inner_faces_begin();
         iter != face->inner_faces_end(); ++iter)
      inner_faces_area += (*iter)->full_area();

    if (face->is_unbounded())
      face->full_area() = 0;
    else {
      Polygon p;
      
      typename Planar_map::Ccb_halfedge_circulator outer_cc =
        face->outer_ccb();
      do {
        p.push_back(outer_cc->target()->point());
      } while(++outer_cc != face->outer_ccb());
      
      face->full_area() = p.area();
    }
    
    //cout<<"p.area()="<<to_double(face->full_area())<<
    //  "  inner_faces_area="<<to_double(inner_faces_area)<<endl;
    
    return face->full_area() - inner_faces_area;
  }
  
  bool  has_common_edge(typename Planar_map::Face_handle face1, 
                        typename Planar_map::Face_handle face2)
  {
    if (face1->is_unbounded() || face2->is_unbounded())
      return true;

    typename Planar_map::Ccb_halfedge_circulator cc1 = face1->outer_ccb();
    typename Planar_map::Ccb_halfedge_circulator cc2 = face2->outer_ccb();

    do{
      do {
        if (cc1->face() == cc2->twin()->face())
          return true;
      } while (++cc2 != face2->outer_ccb());
    } while (++cc1 != face1->outer_ccb());

    return false;
  }

  bool  edge_in_face(typename Planar_map::Halfedge_handle h, 
                     typename Planar_map::Face_handle face)
  {
    if (face->is_unbounded())
      return true;

    Planar_map face_pm;

    typename Planar_map::Ccb_halfedge_circulator cc = face->outer_ccb();

    do{
      face_pm.insert(cc->curve());
    } while (++cc != face->outer_ccb());
    
    typename Planar_map::Locate_type  lt1, lt2;
    face_pm.locate(h->source()->point(),lt1);
    face_pm.locate(h->target()->point(),lt2);
    
    return  (lt1 == Planar_map::FACE && lt2 == Planar_map::FACE);
  }

public:

  virtual void add_edge(const typename Traits::X_curve &, 
                        typename Planar_map::Halfedge_handle e, 
                        bool /* original_direction */, bool overlap = false)
  {
    (void) overlap;
    e->length() = squared_distance(e->source()->point(), e->target()->point());
    e->twin()->length() = e->length();
  }

  virtual void split_edge(typename Planar_map::Halfedge_handle orig_edge, 
                          typename Planar_map::Halfedge_handle new_edge,
                          const typename Traits::X_curve& c1,
                          const typename Traits::X_curve& c2)
  {
    CGAL_assertion_msg(orig_edge->curve() == c1,
                       "first part holds the first curve");
    CGAL_assertion_msg(new_edge->curve() == c2,
                       "second part holds the second curve");

    new_edge->length() = squared_distance(new_edge->source()->point(), 
                                          new_edge->target()->point());
    new_edge->twin()->length() = new_edge->length();
       
    // notice that orig edge is updated except for its length.
    CGAL_assertion_msg(new_edge->length() + 
                       squared_distance(orig_edge->source()->point(), 
                                        orig_edge->target()->point()) ==
                       orig_edge->length(), 
      "the sum of length of both parts equal the length of the original part");

    orig_edge->length() -= new_edge->length();
    orig_edge->twin()->length() = orig_edge->length();
  }

  virtual void split_face(typename Planar_map::Face_handle orig_face, 
                          typename Planar_map::Face_handle new_face)
  {
    //cout<<"---- in split_face"<< endl;

    CGAL_assertion_msg(has_common_edge(orig_face, new_face),
                       "two splited faces have a common edge");
    
    new_face->area() = face_area(new_face);
    
    //cout<<"orig_face area:"<<to_double(orig_face->area())<<endl;
    //cout<<"new_face original area:"<<to_double(new_face->area())<<endl;
    //cout<<"orig_face area:"<<to_double(face_area(orig_face))<<endl;
    

    // asserting that both parts areas equals the original part area.
    CGAL_assertion(orig_face->area() ==
                   new_face->area() + face_area(orig_face));

    orig_face->area() -= new_face->area();
    
    //cout<<"orig_face area:"<<to_double(orig_face->area())<<endl;
    //cout<<"new_face area:"<<to_double(new_face->area())<<endl;
  }

  virtual void add_hole(typename Planar_map::Face_handle in_face, 
                        typename Planar_map::Halfedge_handle new_hole)
  {
    CGAL_assertion_msg(edge_in_face(new_hole, in_face), 
                       "edge of hole is inside face");
  }
  
};

CGAL_END_NAMESPACE         

typedef CGAL::Quotient<CGAL::MP_Float>          NT;
typedef CGAL::Cartesian<NT>                     R;
typedef CGAL::Pm_segment_traits_2<R>            Traits;
typedef CGAL::Pm_measures_dcel<Traits, NT>      Dcel;
typedef CGAL::Planar_map_2< Dcel, Traits >      Planar_map;
typedef Traits::Point                           Point;
typedef Traits::X_curve                         Curve;
typedef Planar_map::Halfedge_handle             Halfedge_handle;
typedef Planar_map::Halfedge_const_iterator     Halfedge_const_iterator;
typedef Planar_map::Locate_type                 Locate_type;


void  check_lengths(const Planar_map& pm)
{
  for (Halfedge_const_iterator h_iter = pm.halfedges_begin(); 
       h_iter != pm.halfedges_end(); ++h_iter)
    CGAL_assertion(h_iter->length() ==
                   CGAL::squared_distance(h_iter->source()->point(), 
                                          h_iter->target()->point()));
}

/*void  check_areas(const Planar_map& pm)
{
  for (Face_const_iterator f_iter = pm.faces_begin(); 
       f_iter != pm.faces_end(); ++f_iter)
    CGAL_assertion(f_iter->area() == face_area(f_iter));
}*/

int main()
{
  Planar_map Pm;
  CGAL::Notification<Planar_map, R>  notf;
  
  int n; 
  std::cin >> n;
  while (n--) {
    double x1, y1, x2, y2;
    std::cin >> x1 >> y1 >> x2 >> y2;
    
    //std::cout << "Inserting ("<< x1 <<","<< y1 <<"--"<< x2 <<","<< y2 <<")"
    //<<std::endl;
    Halfedge_handle hh = Pm.insert(Curve(Point(x1,y1),Point(x2,y2)), &notf);
    std::cout << "Inserted ("<< hh->curve() <<")"<<std::endl;
  }
  CGAL_assertion(Pm.is_valid());

  // General test
  // This used to be tst21, (Shai, Aug. 03, 2000)

  std::cout << "Faces"  << Pm.number_of_faces() << std::endl;
  std::cout << "Halfedges " << Pm.number_of_halfedges() << std::endl;
  std::cout << "Vertices " << Pm.number_of_vertices() << std::endl;

  check_lengths(Pm);
  
  return 0;
}
