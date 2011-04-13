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
// file          : src/PS_facet_3.C
// package       : PS_Stream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#include <CGAL/IO/PS_facet_3.h>

CGAL_BEGIN_NAMESPACE

PS_facet_3::PS_facet_3(vector<Point3> &V,Color
		       edge_color,Color face_color,FILLING f,int number,bool m)
  : _face_color(face_color), _number_of_edge(V.size()), _filling(f),
    _grey_level(0), _number_of_the_facet(number), _mark(m)
{
  unsigned int i=0;
  for(i=0;i<V.size()-1;i++) { 
    PS_edge_3 arete(V[i],V[i+1],edge_color); 
    Face_Arete.push_back(arete);
  }
  PS_edge_3 last_arete(V[V.size()-1],V[0],edge_color);
  Face_Arete.push_back(last_arete);
  arr=Arr_2();
  check_projection();
}

PS_facet_3::PS_facet_3(vector<PS_edge_3> &V,Color color,FILLING f,int number=0,bool m) : 
  Face_Arete(V),_face_color(color),_number_of_edge(V.size()),_filling(f),_grey_level(0),_number_of_the_facet(number),_mark(m)
{check_projection();}

// Copy constructor
PS_facet_3::PS_facet_3(const PS_facet_3 &face){
Face_Arete=face.get_Vect_Arete();
_face_color=face.color();
_number_of_edge=face.number_of_edge();
_filling=face.filling();
_grey_level=face.gray_level();
_projection=face.projection();
_mark=face.mark();
mapp=face.map(); 
_number_of_the_facet=face.number_of_the_facet();
}

// Acces operator for the i edge of the facet
const PS_edge_3::PS_edge_3 PS_facet_3::operator[](int i) const throw(BadIndex) {
  if ((i>=0) && (i<_number_of_edge)) return Face_Arete[i];
  else throw BadIndex();
}

// Assignment/Access operator for the i edge of the facet
PS_edge_3::PS_edge_3 PS_facet_3::operator[](int i) throw(BadIndex) {
   if ((i>=0) && (i<_number_of_edge)) return Face_Arete[i];
  else throw BadIndex();
}

// Assignment operator
PS_facet_3 PS_facet_3::operator= (const PS_facet_3 face) {

Face_Arete=face.get_Vect_Arete();
_face_color=face.color();
_number_of_edge=face.number_of_edge();
_filling=face.filling();
_grey_level=face.gray_level();
_projection=face.projection();
mapp=face.map();
_mark=face.mark(); 
_number_of_the_facet=face.number_of_the_facet();
return *this;
}

 //Assignment/Access operator according a PS_edge_3(Useful to define
//a PS_edge_3 like a facet(and also a point3 to be a facet)
PS_facet_3 PS_facet_3::operator= (const PS_edge_3 ar) {
  
  //Il faut vider le vecteur Face_arete car ici on ne doit avoir
  //qu'une arete dans le vecteur // A FAIRE
  Face_Arete.clear();
  Face_Arete.push_back(ar);
  _face_color=ar.color();
  _number_of_edge=1;
  return *this;
}

ostream& operator<<(ostream& os, const PS_facet_3& face) {
  os << "Facet number : " << face.number_of_the_facet() << endl;
  os << "Facet color : " << face.color() << " | Number of edge : " << face.number_of_edge();
  os << "\n"; 
  for (int i=0;i<face.number_of_edge();i++){
    os << "Edge " << i << " :" << face[i] << endl;
  }
return os;
} 

coord_type PS_facet_3::xmin() {
  coord_type xmin= Face_Arete[0].xmin();
  int i;

  for(i=1;i<_number_of_edge;i++) {
     if (xmin > Face_Arete[i].xmin())
       xmin= Face_Arete[i].xmin();
  }
  return xmin;
}

coord_type PS_facet_3::ymin() {
  coord_type ymin= Face_Arete[0].ymin();
  int i;

  for(i=1;i<_number_of_edge;i++) {
     if (ymin > Face_Arete[i].ymin())
       ymin= Face_Arete[i].ymin();
  }
  return ymin;
}

coord_type PS_facet_3::zmin() {
  coord_type zmin= Face_Arete[0].zmin();
  int i;

  for(i=1;i<_number_of_edge;i++) {
     if (zmin > Face_Arete[i].zmin())
       zmin= Face_Arete[i].zmin();
  }
  return zmin;
}



coord_type PS_facet_3::xmax() {
  coord_type xmax= Face_Arete[0].xmax();
  int i;

  for(i=1;i<_number_of_edge;i++) {
     if (xmax < Face_Arete[i].xmax())
       xmax= Face_Arete[i].xmax();
  }
  return xmax;
}

coord_type PS_facet_3::ymax() {
  coord_type ymax= Face_Arete[0].ymax();
  int i;

  for(i=1;i<_number_of_edge;i++) {
     if (ymax < Face_Arete[i].ymax())
       ymax= Face_Arete[i].ymax();
  }
  return ymax;
}

coord_type PS_facet_3::zmax() {
  coord_type zmax= Face_Arete[0].zmax();
  int i;

  for(i=1;i<_number_of_edge;i++) {
     if (zmax < Face_Arete[i].zmax())
       zmax = Face_Arete[i].zmax();
  }
  return zmax;
}

coord_type PS_facet_3::max() {
  
  coord_type zmax=(*this).zmax();
  coord_type ymax=(*this).ymax();
  coord_type xmax=(*this).xmax();
 
  vector<coord_type> v;
  v.push_back(xmax);
  v.push_back(ymax);
  v.push_back(zmax);

  sort(v.begin(),v.end());

  return v[2];
}

coord_type PS_facet_3::min() {
  coord_type zmin=(*this).zmin();
  coord_type ymin=(*this).ymin();
  coord_type xmin=(*this).xmin();
 
  vector<coord_type> v;
  v.push_back(xmin);
  v.push_back(ymin);
  v.push_back(zmin);

  sort(v.begin(),v.end());
  return (v[0]);
}

//Check if P(this) is behind Q
//P is behind Q if all the P vertex are behind the Q plane 
//Precondition : only with visible facet

bool PS_facet_3::is_back_side(PS_facet_3& Q) {
Plane3 planQ(Q[0].first_point(),Q[0].second_point(),Q[Q.number_of_edge()-1].first_point());

  for (int i=0;i<_number_of_edge;i++) {
    if ((planQ.has_on_positive_side(Face_Arete[i].first_point())) ||
	(planQ.has_on_positive_side(Face_Arete[i].second_point()))) {
      return false;
    } 
  }
  return true;
}

 // front_side check if P(this) is in front of Q
/* Precondition : only with visible facet */
bool PS_facet_3::is_front_side(PS_facet_3& Q) {
  Plane3 planQ(Q[0].first_point(),Q[0].second_point(),Q[Q.number_of_edge()-1].first_point());
  
  for (int i=0;i<_number_of_edge;i++) {
    if ((planQ.has_on_negative_side(Face_Arete[i].first_point())) ||
	(planQ.has_on_negative_side(Face_Arete[i].second_point()))) {
      return false;
    }
  }
  return true;
}

bool PS_facet_3::is_z_overlap(PS_facet_3& Q) {
  
  if ((zmax() < Q.zmin()) || (zmin() > Q.zmax())) {
 return false;
  }
  return true;
}

bool PS_facet_3::is_x_overlap(PS_facet_3& Q) {
  
  if ((xmax() < Q.xmin()) || (xmin() > Q.xmax())) {
    return false;
  }
  return true;
}

bool PS_facet_3::is_y_overlap(PS_facet_3& Q) {
  
  if ((ymax() < Q.ymin()) || (ymin() > Q.ymax())) {
    return false;
  }
  return true;
}



bool PS_facet_3::has_a_cycle(PS_facet_3& Q){
  return ((_mark==true) && (Q.mark()==true));
}


bool PS_facet_3::two_test(PS_facet_3& Q){
   
  if (is_back_side(Q)) {return true;}
  else if (Q.is_front_side(*this)) {return true;}
  else return false;
}

//Check if P(this) can be draw before Q
//Return true if one of the five test is true

bool PS_facet_3::five_test(PS_facet_3& Q){
   
  if(!is_z_overlap(Q)) {return true;}
  else if (!is_x_overlap(Q)) {return true;}
  else if (!is_y_overlap(Q)) {return true;}
  else if (is_back_side(Q)) {return true;}
  else if (Q.is_front_side(*this)) {return true;}
  else if (!is_xy_plan_overlap(Q)) {return true;}
  else return false;
}

void PS_facet_3::check_projection() {

Plane3 P(Face_Arete[0].first_point(),Face_Arete[0].second_point(),Face_Arete[_number_of_edge-1].first_point());

coord_type a=P.a();
coord_type c=P.c();

Vector3 v=P.orthogonal_vector();
 
 if((a==0) && (c==0)) {_projection=ZX;}
 else if(v.z()==0) {_projection=YZ;}
 else {_projection=XY;}
}

Arr_Point2 PS_facet_3::projection(Point3 q){
  
  switch(_projection) {
  case ZX : 
    {
      Arr_Point2 p(q.x(),q.z());
      return p;
    }
  case YZ : 
    {
      Arr_Point2 p(q.y(),q.z());
      return p; 
    }
  default :
    {
      
      Arr_Point2 p(q.x(),q.y());
      return p; 
    }
  }
}

Point3 PS_facet_3::lift(Arr_Point2 p) {

Plane3 P(Face_Arete[0].first_point(),Face_Arete[0].second_point(),Face_Arete[_number_of_edge-1].first_point());

coord_type a=P.a();
coord_type b=P.b();
coord_type c=P.c();
coord_type d=P.d();

  switch(_projection) {
  case ZX : 
    {
       Point3 q(CGAL::to_double(p.x()),
	       (-d-a* CGAL::to_double(p.x())-c*CGAL::to_double(p.y()) )/b,
		CGAL::to_double(p.y()));   

      return q;
    }
  case YZ : 
    {
    Point3 q( ((-d-b*CGAL::to_double(p.x())-c*CGAL::to_double(p.y()))/a),
	      CGAL::to_double(p.x()),
	      CGAL::to_double( p.y()));  

      return q; 
    }
  default :
    {
      Point3 q(CGAL::to_double( p.x()),
	      CGAL::to_double(p.y()),
	      ((-d-a*CGAL::to_double(p.x())-b*CGAL::to_double(p.y()))/c));   

     return q; 
    }
  }
}

 //ATTENTION A PRENDRE LA PARTIE ENTIERE CAR IL FAUT UN INT pour la ligne
// Make the arrangement of the face
int PS_facet_3::make_arrangement(PS_edge_3 &l) {
  arr=Arr_2();
  check_projection();
 
  Arr_Point2 p1 = projection(l.first_point());
  Arr_Point2 q1 = projection(l.second_point());
  Curve *c1 = new Curve(p1,q1);
  Arr_2::Curve_iterator cit = arr.insert(*c1);
  pair<Curve_node *,PS_edge_3 *> pair1( &(*cit),&l);
  mapp.insert(pair1);
  Arr_Point2 initial = projection(Face_Arete[0].first_point());
  Arr_Point2 first = (*&initial);
  int i;
  for(i=0;i<_number_of_edge-1;i++) {
    Arr_Point2 second = projection(Face_Arete[i].second_point());
    Curve *c2 = new Curve(first,second);
    cit = (*this).arr.insert(*c2);
    pair<Curve_node *,PS_edge_3 *> pair2(&(*cit),&(Face_Arete[i]));
    mapp.insert(pair2);
    first = second;
   }
   Curve *cf = new Curve(first,initial);
   cit = (*this).arr.insert(*cf);
   pair<Curve_node *,PS_edge_3 *> pairf(&(*cit),&(Face_Arete[_number_of_edge-1]));
   mapp.insert(pairf);
  if (arr.number_of_faces() >2) {return 1;}
  else return 0;
}

//To recover the facet of the arrangement
vector<PS_facet_3> PS_facet_3::create_facets(PS_edge_3 &l) {
  
  int nb_face=1;
  vector<PS_facet_3> facets;
  Arr_2::Face_iterator fit;
  if (arr.number_of_faces() >2) {
    std::map<Curve_node *,PS_edge_3 *>::iterator mi;
    for(fit=arr.faces_begin();fit!=arr.faces_end();++fit) {
      //the face is bounded and there is no hole inside the face
      if ( (!(*fit).is_unbounded()) && ( (*fit).holes_begin()== (*fit).holes_end()) ) {
	Arr_2::Ccb_halfedge_circulator initial=fit->outer_ccb(); 
	Arr_2::Ccb_halfedge_circulator e=fit->outer_ccb(); 
	//We have to check if the edge is not an edge create by the
	//cutting of the facet, because it is not possible to determine
	//if the facet is inner or outer.So if it is the case ,we take
	//another edge
	mi = mapp.find(&( * (e->edge_node()->curve_node())));
	while (mi->second == &l) {
	  e++;initial++;	
	  mi = mapp.find(&( * (e->edge_node()->curve_node())));
	}
	
	//Here I am sure that the edge that I check is not provide by an
	//outer face or a hole,and moreover is not on the curve of the
	//cutting line.
	mi = mapp.find(&( *(e->edge_node()->curve_node())));
	
	Arr_Point2 source = (*(e->source())).point();
	Arr_Point2 target = (*(e->target())).point(); 
	Point2 source_double(CGAL::to_double(source.x()),CGAL::to_double(source.y()));
	Point2 target_double(CGAL::to_double(target.x()),CGAL::to_double(target.y()));
	
	Arr_Point2 curve_source = (*mi->first).curve().source();
	Arr_Point2 curve_target = (*mi->first).curve().target(); 
	Point2 curve_source_double(CGAL::to_double(curve_source.x()),CGAL::to_double(curve_source.y()));
	Point2 curve_target_double(CGAL::to_double(curve_target.x()),CGAL::to_double(curve_target.y()));
	e++;
	if( mi->second == &l) {
#ifdef DEBUG
       cout << "I have found that I compare mi->second : " <<
	 mi->second << " | and the cutting line : -> CUTTING LINE " << &l << endl;
       ///////A REVOIR -> il faut tester avec le prochain mi->second
       //-> e++;mi = mapp.find(&( * (e->edge_node()->curve_node())));
#endif
     
     }
     Vector2 v1(target_double-source_double);
     Vector2 v2(curve_target_double-curve_source_double);

     //It is necessary to test if the halfedge and the PS_edge_3
     //corresponding got the same orientation. If not,the face is an
     //outer face and must be reject.To check this we make the scalar
     //product between the two vectors (the vector of the edge that we 
     //check and the curve_node of the arrangement which the edge
     //comes from).
     //If >0 the two edges got the same orientation so we keep the facet
     //else they don't have the same orientation so we reject the facet
     if (v1*v2>0) {
       //It is necessary to reinitialize the edge iterator.if not It
       //miss the first edge use to know if the facet got a good orientation
       e=initial;
       std::vector<PS_edge_3> v_edge;	 
       
       do {
	 mi = mapp.find(&( * (e->edge_node()->curve_node())));
	 Arr_Point2 p1 = (*(e->source())).point();
	 Arr_Point2 p2 = (*(e->target())).point();
	 Point3 p1_3D = lift(p1);
	 Point3 p2_3D = lift(p2);
	 PS_edge_3 edge(p1_3D,p2_3D, (*(mi->second)).color(),(*(mi->second)).visibility());
	 v_edge.push_back(edge);
	 e++;
       }
       while (e != initial);
       //When the new facet is create we give it the number of the
       //mother facet.
       PS_facet_3 f(v_edge,color(),filling(),number_of_the_facet());
       facets.push_back(f);
     }
      }
      nb_face++;
    }
    return facets;
  }
  else {
    return facets;
  }
} 

bool PS_facet_3::is_xy_plan_overlap(PS_facet_3 &Q) {
  
  PS_facet_3 f1(*this);
  PS_facet_3 f2(Q);
  
  f1.set_projection(XY);
  f2.set_projection(XY);
  
  Poly poly1;
  Poly poly2;
  
  for (int i=0;i<f1.number_of_edge();i++) {
    Point2_Leda p(f1[i].first_point().x(),f1[i].first_point().y());
    poly1.push_back(p);
  }
 
  for (int i=0;i<f2.number_of_edge();i++) {
    Point2_Leda p(f2[i].first_point().x(),f2[i].first_point().y());
    poly2.push_back(p);	   
  }
  
  Poly::Edge_const_iterator eit1;
  Poly::Edge_const_iterator eit2;
  
  for(eit1=poly1.edges_begin();eit1!=poly1.edges_end();eit1++) {
    for (eit2=poly2.edges_begin();eit2!=poly2.edges_end();eit2++){
      Segment2 poly1_segment=*eit1;
      Segment2 poly2_segment=*eit2;
      if (do_intersect(poly1_segment,poly2_segment)) {
	
	//Il faut tester les determinants
	// P et Q      
	Point2_Leda P(poly1_segment.source());
	Point2_Leda Q(poly1_segment.target());
	// R et S      
	Point2_Leda R(poly2_segment.source());
	Point2_Leda S(poly2_segment.target());
	
	int signPQR = sign(CGAL::to_double(det(P,Q,R)));
	int signPQS = sign(CGAL::to_double(det(P,Q,S)));
	int signRSQ = sign(CGAL::to_double(det(R,S,Q)));
	int signRSP = sign(CGAL::to_double(det(R,S,P)));
	
	if ((signPQR!=signPQS) && (signRSQ!=signRSP) && 
	    (signPQR!=0) && (signPQS!=0) && (signRSQ!=0) &&
	    (signRSP!=0)) {
	  return true;
	}
	// else I do nothing because I'am in a degenerated case
      }
    }
  }
  
  //Test to know if one of the two face is include in the other
  // Poly::Vertex_const_iterator vit1=poly1.vertices_begin();
  //Poly::Vertex_const_iterator vit2=poly2.vertices_begin();
  //cout << "JE SUIS AVT TEST DU UNBOUNDED " << endl;
  
  // PAS BON A REVOIR !!!! CAS A TESTER OBLIGATOIREMENT
  //  if ( (!(poly1.has_on_unbounded_side(*vit2))) ||
  //       (!(poly2.has_on_unbounded_side(*vit1)))) {
  //    cout << "J'ai renvoye vrai dans unbounded " <<endl;
  //    return true;
  //  }
  return false;
}

void PS_facet_3::transformation(Transformation &t) {
  for (int i=0;i<_number_of_edge;i++) {
    Face_Arete[i].transformation(t);
  }
}

CGAL_END_NAMESPACE
