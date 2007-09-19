#ifndef CGAL_SPHERE_3_TABLE_H
#define CGAL_SPHERE_3_TABLE_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_key.h>
#include <CGAL/Tools/Log.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

template <class Traits_t>
class Sphere_3_table: public CGAL::Kinetic::Ref_counted<Sphere_3_table<Traits_t> > {
public:
  typedef Traits_t Traits;
  
  typedef typename Traits::Geom_traits Geom_traits;

  typedef Sphere_key Key;
  typedef typename Geom_traits::Sphere_3 Sphere_3;
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Plane_3 Plane_3;
  typedef typename Geom_traits::Vector_3 Vector_3;
  typedef typename Geom_traits::Segment_3 Segment_3;
  typedef typename Geom_traits::Line_3 Line_3;
  typedef typename Geom_traits::Line_2 Line_2;
  typedef typename Geom_traits::Point_2 Point_2;
  typedef typename Geom_traits::Vector_2 Vector_2;
  typedef typename Geom_traits::Circle_2 Circle_2;
  typedef typename Geom_traits::Segment_2 Segment_2;
  typedef typename std::vector<Sphere_3> Spheres;

  struct Equal_centers_exception {
  };

  CGAL_GETOBJECT(Geom_traits, geom_traits, return t_);

  template <class It> 
  Sphere_3_table(It bs, It es): has_temp_(false){
    spheres_.reserve(std::distance(bs, es)+3);
    initialize_1();
    spheres_.insert(spheres_.end(), bs, es);
    initialize_2();
  }


  std::ostream &write(std::ostream &out) const ;

  void set_temp_sphere(const Sphere_3 &s) {
    has_temp_=true;
    spheres_[0]=s;
  }

  // really just debugging
  Key new_sphere(const Sphere_3 &s) {
    spheres_.push_back(s);
    return Key(spheres_.size()-4);
  }

  CGAL_GETNR(unsigned int, size, return spheres_.size()-3);

  CGAL_CONST_ITERATOR(Sphere_3, sphere_3, typename Spheres::const_iterator,
		      return spheres_.begin()+3,
		      return spheres_.end());
  
  Sphere_3 operator[](Key k) const {
    //CGAL_precondition(static_cast<unsigned int> (k.index()+4) < spheres_.size());
    //CGAL_precondition(k.index()+4 >= 0);
    CGAL_check_bounds(static_cast<unsigned int>(k.internal_index()), static_cast<unsigned int>(0U), spheres_.size());
    return spheres_[k.internal_index()];
  }


  CGAL_SIZE(sphere_3s, return spheres_.size());

  struct Sphere_key_iterator_t{
    typedef Sphere_key_iterator_t This;
    typedef Key value_type;
    typedef const Key &reference;
    typedef const Key* pointer;
    typedef size_t difference_type;
    typedef std::random_access_iterator_tag iterator_category;
    Sphere_key_iterator_t(){}
    Sphere_key_iterator_t(int i): k_(i){}
    value_type operator*() const {return k_;}
    pointer operator->() const {return &k_;}
    Sphere_key_iterator_t operator++() {
      k_= Key(k_.input_index()+1);
      return *this;
    }
    Sphere_key_iterator_t operator++(int) {
      Sphere_key_iterator_t ret=*this;
      operator++();
      return ret;
    }
    difference_type operator-(Sphere_key_iterator_t o) const {
      return k_.internal_index() - o.k_.internal_index();
    }
    void operator+=(int i) {
      k_= Key(k_.input_index()+i);
    }
    CGAL_COMPARISONS1(k_);
  private:
    Key k_;
  };


  CGAL_CONST_ITERATOR(Sphere_key, sphere_key, Sphere_key_iterator_t, 
		return Sphere_key_iterator_t(0),
		return Sphere_key_iterator_t(spheres_.size()-3));

  
  FT discriminant(Key i) const;

  CGAL_GET(Bbox_3, bbox_3, return bbox_);


  /* 
     linear constructions----------------------------------------------------
  */

  CGAL_GET(FT, inf, return center(Key(Key::UB))[1]);

  // the point described by the vertex (a,b) should be on the positive side
  Plane_3 separating_plane(Key a, Key b) const ;

  // point from the second to the first
  Plane_3 equipower_plane(Key a, Key b) const ;
  
  Point_3 equipower_point(Key a, Key b) const ;

  Line_3 equipower_line(Key a, Key b) const ;

  Point_3 center(Key ind) const;
  FT squared_radius(Key ind) const;

  Sphere_3 sphere(Key ind) const;

  void initialize_1();
  void initialize_2();


  /*
    event data, not really appropriate, but it is the best place to put it
   */
  struct Event_pair_data {
    typedef typename Traits::Event_point_3 EP;
    Event_pair_data():index_(-1){}
    int index_;
    EP events_[2];
  };
  struct Triple_data {
    Event_pair_data srr_events_;
  };
  struct Pair_data {
    typedef typename Traits::Event_point_3 EP;
    struct KC_pair{
      typedef KC_pair This;
      KC_pair(Key k, Coordinate_index c): k_(k), c_(c){};
      Key k_;
      Coordinate_index c_;
      CGAL_COMPARISONS2(k_, c_);
    };
    
    std::map<KC_pair, Event_pair_data> cxr_events_;
  };

  

  /*UPair_data &upair_data(const Sphere_key_upair &t) {
    return upair_data_[t];
    }*/

  Pair_data &pair_data(const Sphere_key_pair &t) {
    return pair_data_[t];
  }

  Triple_data &triple_data(const Sphere_key_triple &t) {
    return triple_data_[t];
  }




  std::map< Sphere_key_triple, Triple_data> triple_data_;
  std::map< Sphere_key_pair, Pair_data> pair_data_;
  //std::map< Sphere_key_upair, Pair_data> pair_data_;

  Spheres spheres_;
  bool has_temp_;

  Geom_traits t_;
  typename Geom_traits::Intersect_3 di_;
  Bbox_3 bbox_;
};

CGAL_OUTPUT1(Sphere_3_table);


CGAL_AOS3_END_INTERNAL_NAMESPACE

#include <CGAL/Arrangement_of_spheres_3/Sphere_3_table_impl.h>


#endif
