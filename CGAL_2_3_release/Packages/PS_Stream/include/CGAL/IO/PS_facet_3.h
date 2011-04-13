// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/IO/PS_facet_3.h
// package       : PS_Stream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_PS_FACET_3_H
#define CGAL_PS_FACET_3_H

#include <CGAL/basic.h>

#include <algorithm>
#include <map>
#include <utility>
#include <vector>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/IO/PS_edge_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Aff_transformation_3.h>

//Include for the polygon
#include <CGAL/Polygon_2.h>

//Include for the Arrangement
#include <CGAL/Arithmetic_filter.h>
#include <CGAL/leda_real.h>
#include <CGAL/double.h>
#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Arrangement_2.h>

CGAL_BEGIN_NAMESPACE

enum FILLING{NO_FILL,WIRED_CULLBACK_FACING,UNIFORM_FILL,NORMAL_FILL};
enum PROJECTION{XY,YZ,ZX};

typedef double coord_type;
//typedef leda_real coord_type;
typedef CGAL::Cartesian< coord_type > CT;
typedef CGAL::Cartesian < leda_real > R;
typedef Point_2< CT > Point2;
typedef Point_3< CT > Point3;
typedef Plane_3< CT > Plane3;
typedef Line_3 < CT > Line3;
typedef Vector_3< CT > Vector3;
typedef Vector_2< CT > Vector2;
typedef Aff_transformation_3< CT > Transformation;

//The polygon for the xy_plan_overlap got Point in leda_real
typedef Point_2< R > Point2_Leda;
typedef Polygon_traits_2< R > Traits_poly;
typedef Traits_poly::Segment_2 Segment2;
typedef Traits_poly::Vector_2 Vector2_Leda;
typedef vector<Point2_Leda> Container_poly;
typedef Polygon_2<Traits_poly,Container_poly> Poly;

//Arrangement is in LEDA_REAL
typedef CGAL::Cartesian< leda_real > Arr_CT;
typedef Arr_CT   Arr_Rep;
typedef CGAL::Arr_segment_exact_traits<Arr_Rep>           Traits_arr;
typedef Traits_arr::Point                                 Arr_Point2;
typedef Traits_arr::X_curve                               X_curve;
typedef Traits_arr::Curve                                 Curve;
typedef CGAL::Arr_base_node<Curve>                     Base_node;
typedef CGAL::Arr_2_default_dcel<Traits_arr>               Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits_arr,Base_node >    Arr_2;
typedef Arr_2::Curve_node   Curve_node;

class PS_facet_3 {

public:

  //Exception due to a bad index 
  class BadIndex : public std::exception {};
  //Constructors

  //Constructor with a list of point3
  //Precondition : the list must contain the points in the
  //                trigonometric positive way
  PS_facet_3(vector<Point3> &V,Color edge_color = RED ,Color
	     filling_color = BLACK, FILLING
	     f=NO_FILL,int number=0,bool m=false);
  //Constructor with a list of PS_edge_3
  //Precondition : the list must contain the edges in the
  //               trigonometric positive  
  PS_facet_3(vector<PS_edge_3> &V,Color color= BLACK,FILLING
	     f=NO_FILL,int number=0,bool m=false);
  
  // Copy constructor
  PS_facet_3(const PS_facet_3 &face);

  //Inline functions

  //Accessors
  //Return the number of edge of the facet.(Visible+Invisible)
  inline int number_of_edge() const {return _number_of_edge;}
  //Return the facet color
  inline Color color() const {return _face_color;}
  //Return the gray level of the facet
  inline double gray_level() const {return _grey_level;}
  //Return the PS_edge_3 vector
  inline vector<PS_edge_3> get_Vect_Arete() const { return Face_Arete;}
  //Return the type of filling of the facet
  inline FILLING filling() const {return _filling;}
  //Return the type of projection of the facet
  inline PROJECTION projection() const {return _projection;}
  //Return the mark of the facet
  inline bool mark() const {return _mark;}
  //Return the map of the facet
  inline map<Curve_node *,PS_edge_3 *> map() const {return mapp;}
  //Return the number of the facet
  inline int number_of_the_facet() const {return _number_of_the_facet;}
  
  //Settings
  //Modify the type of filling of the facet
  inline void set_filling(FILLING f) {_filling=f;}
  //Modify the color of the facet
  inline void set_color(Color c) {_face_color=c;}
  //Modify the gray level.
  inline void set_gray_level(double c) {_grey_level=c;}
  //Modify the type of projection of the facet
  inline void set_projection(PROJECTION p) {_projection=p;}
  //Modify the mark of the facet
  inline void set_mark(bool m) {_mark=m;}
  //Modify the number of the facet
  inline void set_number_of_the_facet(int i) {_number_of_the_facet=i;}
  
  // Operators
  
  // Access operator for the i edge of the facet
  const PS_edge_3 operator[](int i) const throw(BadIndex);
  // Assignment/Access operator for the i edge of the facet
  PS_edge_3 operator[](int i) throw (BadIndex);
  // Assignment operator
  PS_facet_3 operator= (const PS_facet_3 face);
  //Assignment/Access operator according a PS_edge_3(Useful to define
  //a PS_edge_3 like a facet(and also a point3 to be a facet)
  PS_facet_3 operator= (const PS_edge_3 ar);
  
  friend ostream& operator<<(ostream& os,const PS_facet_3 &face);

  // Others functions

  //Allows to add a new edge in the face
  void add(PS_edge_3 &arete){Face_Arete.push_back(arete);_number_of_edge++;}
 
  //Return the minimum or the maximun of the edges of the facet
  //according the x,y,z 
  coord_type xmin();  coord_type ymin();  coord_type zmin(); 
  coord_type xmax();  coord_type ymax();  coord_type zmax(); 

  //Return the maximum(minimum) of the x,y,z maximum(minimum)
  coord_type max(); coord_type min();
  //Check if the two facet P(this) and Q overlap in x
  bool is_x_overlap(PS_facet_3& Q);
  //Check if the two facet P(this) and Q overlap in y
  bool is_y_overlap(PS_facet_3& Q);
  //Check if the two facet P(this) and Q overlap in z
  bool is_z_overlap(PS_facet_3& Q);
  //Check if the facet P(this) is behind the facet Q
  bool is_back_side(PS_facet_3& Q);
  //Check if the facet P(this) is in front of the facet Q
  bool is_front_side(PS_facet_3& Q);
 
  //Compute the second two test of the Depth Sort Algorithm
  bool two_test(PS_facet_3& Q);
  //Compute the first five test of the Depth Sort Algorithm
  bool five_test(PS_facet_3& Q);
  //test to know if there is a cycle
  bool has_a_cycle(PS_facet_3& Q);
  
  //Check the type of projection for the facet
  void check_projection();
  //Make the projection of a Point3 according the PROJECTION variable
  Arr_Point2 projection(Point3 q);
  //Give the point3 which correspond to the point2 projected
  Point3 lift(Arr_Point2 p);
  //Make the arrangement of the face with the Line3 l which cut the face
  int make_arrangement(PS_edge_3 &l);
  //To recover the facet of the arrangement
  vector<PS_facet_3> create_facets(PS_edge_3 &l);
  //Check if the projection of the two PS_facet_3 overlap in the plan XY
  bool is_xy_plan_overlap(PS_facet_3 &Q);
  //Make the transformation to z-axis view
  void transformation(Transformation &t); 

  //Return the cross product of the two vectors
  leda_real cp(Vector2_Leda p,Vector2_Leda q) {
    return (p.x() * q.y() - q.x() * p.y());
  }
  
  //Return the determinant of the 3 point p,q and r
  leda_real det(Point2_Leda p,Point2_Leda q,Point2_Leda r) {
    return cp(p-q,p-r);
  } 

  //Return the sign of x  
  int sign(double x) {
    if (x>0) {return 1;}
    else if (x<0) { return -1;}
    else {return 0;}
}
  

private:
 
  //Vector of edge of the facet
  vector<PS_edge_3> Face_Arete;
  //Facet color
  Color _face_color;
  //Number of edge of the facet;
  int _number_of_edge;
  //Type of filling of the facet
  //NO_FILL,UNIFORM_FILL,NORMAL_FILL
  FILLING _filling;
  //Gray level of the facet
  double _grey_level;
  //Type of projection for the facet
  PROJECTION _projection;
  //Arrangement of the face
  Arr_2 arr;
  //Map 
  std::map<Curve_node *,PS_edge_3 *> mapp;
  //Number of the facet : use to avoid cutting two facets that comes
  //from the same facet 
  int _number_of_the_facet;  
  //Mark to avoid cycle
  bool _mark;
};

CGAL_END_NAMESPACE

#endif // CGAL_PS_FACET_3_H
