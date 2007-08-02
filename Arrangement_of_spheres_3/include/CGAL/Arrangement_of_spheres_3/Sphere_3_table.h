#ifndef CGAL_SPHERE_3_TABLE_H
#define CGAL_SPHERE_3_TABLE_H
#include <CGAL/Arrangement_of_spheres_3_basic.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Arrangement_of_spheres_3/Sphere_key.h>
#include <CGAL/Tools/Log.h>

CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

CGAL_AOS3_TEMPLATE
class Sphere_3_table: public CGAL::Kinetic::Ref_counted<Sphere_3_table CGAL_AOS3_TARG > {
public:
#ifdef CGAL_AOS3_USE_TEMPLATES
  CGAL_AOS3_TRAITS;
public:
  typedef typename Traits::Geom_traits Geom_traits;
#else
  typedef Arrangement_of_spheres_3_geom_traits Geom_traits;
#endif
  typedef Sphere_key Key;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Sphere_3 Sphere_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::FT FT;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Point_3 Point_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Plane_3 Plane_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Vector_3 Vector_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Segment_3 Segment_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Line_3 Line_3;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Line_2 Line_2;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Point_2 Point_2;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Vector_2 Vector_2;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Circle_2 Circle_2;
  typedef CGAL_AOS3_TYPENAME Geom_traits::Segment_2 Segment_2;
  typedef CGAL_AOS3_TYPENAME std::vector<Sphere_3> Spheres;

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

  CGAL_CONST_ITERATOR(Sphere_3, sphere_3, CGAL_AOS3_TYPENAME Spheres::const_iterator,
		return spheres_.begin()+3,
		return spheres_.end());

  Sphere_3 operator[](Key k) const {
    //CGAL_precondition(static_cast<unsigned int> (k.index()+4) < spheres_.size());
    //CGAL_precondition(k.index()+4 >= 0);
    CGAL_check_bounds(k.internal_index(), 0, spheres_.size());
    return spheres_[k.internal_index()];
  }


  CGAL_SIZE(sphere_3s, return spheres_.size());

  struct Sphere_key_iterator_t{
    typedef Sphere_key_iterator_t This;
    typedef Key value_type;
    typedef const Key &reference_type;
    typedef const Key* pointer_type;
    typedef size_t difference_type;
    Sphere_key_iterator_t(){}
    Sphere_key_iterator_t(int i): k_(i){}
    value_type operator*() const {return k_;}
    pointer_type operator->() const {return &k_;}
    Sphere_key_iterator_t operator++() {
      k_= Key(k_.input_index()+1);
      return *this;
    }
    Sphere_key_iterator_t operator++(int) {
      Sphere_key_iterator_t ret=*this;
      operator++();
      return ret;
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

  CGAL_GET(FT, inf, return center(Key(Key::TR))[1]);

  // the point described by the vertex (a,b) should be on the positive side
  Plane_3 separating_plane(Key a, Key b) const ;

  // point from the second to the first
  Plane_3 equipower_plane(Key a, Key b) const ;
  
  // point from the second to the first
  Point_3 equipower_point(Key a, Key b) const ;

  Point_3 center(Key ind) const;
  FT squared_radius(Key ind) const;

  Sphere_3 sphere(Key ind) const;

  void initialize_1();
  void initialize_2();

  Spheres spheres_;
  bool has_temp_;

  Geom_traits t_;
  CGAL_AOS3_TYPENAME Geom_traits::Intersect_3 di_;
  Bbox_3 bbox_;
};
#ifdef CGAL_AOS3_USE_TEMPLATES
CGAL_OUTPUT1(Sphere_3_table);
#else
CGAL_OUTPUT(Sphere_3_table);
#endif

CGAL_AOS3_END_INTERNAL_NAMESPACE

#ifdef CGAL_AOS3_USE_TEMPLATES
#include <CGAL/Arrangement_of_spheres_3/Sphere_3_table_impl.h>

#endif

#endif
