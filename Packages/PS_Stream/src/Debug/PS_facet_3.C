#ifndef PS_FACET_3_C
#define PS_FACET_3_C

#include <CGAL/IO/PS_facet_3.h>
using namespace std;
CGAL_BEGIN_NAMESPACE

// PS_facet_3::~PS_facet_3() {
//   delete(Face_Arete);
//   delete(arr);
//   delete(mapp);
// }
  

PS_facet_3::PS_facet_3(vector<Point3> &V,Color
		       edge_color,Color face_color,FILLING f,int number,bool m) : _face_color(face_color), _number_of_edge(V.size()),_filling(f),_grey_level(0),_number_of_the_facet(number),_mark(m) {
  
#ifdef DEBUG
  cout << "########### FACET " << _number_of_the_facet << " ###############" << endl;
  cout << "Number of vertices : " << V.size() << endl;
  cout << "Display of the edges : " << endl; 
#endif
  
  unsigned int i=0;
  for(i=0;i<V.size()-1;i++) { 
    PS_edge_3 arete(V[i],V[i+1],edge_color); 
    
#ifdef DEBUG
    cout << "V["<< i <<"],V["<< (i+1) <<"] : " << V[i] << " -> " << &(V[i])
	 << " | " << V[i+1] << " -> " << &(V[i+1]) << endl;  
#endif
    Face_Arete.push_back(arete);
  }
  PS_edge_3 last_arete(V[V.size()-1],V[0],edge_color);
  
#ifdef DEBUG
  cout << "V["<< (V.size()-1) << "],V[0] : " << V[V.size()-1]
       <<" -> " << &(V[V.size()-1]) << " | " << V[0] << " -> "  << &(V[0]) << endl;
#endif

  Face_Arete.push_back(last_arete);
  arr=Arr_2();
  check_projection();
#ifdef DEBUG
  cout << endl;
#endif
  
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
#ifdef DEBUG
      cout << "P n'est pas entierement derriere Q " << endl;
#endif
      return false;
    } 
  }
#ifdef DEBUG
 cout << "P est entierement derriere Q " << endl;
#endif

  return true;
}

 // front_side check if P(this) is in front of Q
/* Precondition : only with visible facet */
bool PS_facet_3::is_front_side(PS_facet_3& Q) {
  Plane3 planQ(Q[0].first_point(),Q[0].second_point(),Q[Q.number_of_edge()-1].first_point());
  
  for (int i=0;i<_number_of_edge;i++) {
    if ((planQ.has_on_negative_side(Face_Arete[i].first_point())) ||
	(planQ.has_on_negative_side(Face_Arete[i].second_point()))) {

#ifdef DEBUG
      cout << "P n'est pas entierement devant Q " << endl;
#endif
      return false;
    }
  }

#ifdef DEBUG
    cout << "P est entierement devant Q " << endl;
#endif
  return true;
}

bool PS_facet_3::is_z_overlap(PS_facet_3& Q) {
  
  if ((zmax() < Q.zmin()) || (zmin() > Q.zmax())) {
#ifdef DEBUG
    cout<<"les z ne se chevauchent pas"<<endl;
#endif
 return false;
  }

#ifdef DEBUG
  cout<<"les z se chevauchent"<<endl;
#endif
  return true;
}

bool PS_facet_3::is_x_overlap(PS_facet_3& Q) {
  
  if ((xmax() < Q.xmin()) || (xmin() > Q.xmax())) {
#ifdef DEBUG
    cout<<"les x ne se chevauchent pas"<<endl;
#endif
    return false;
  }
#ifdef DEBUG
  cout<<"les x se chevauchent"<<endl;
#endif
  return true;
}

bool PS_facet_3::is_y_overlap(PS_facet_3& Q) {
  
  if ((ymax() < Q.ymin()) || (ymin() > Q.ymax())) {
#ifdef DEBUG
    cout<<"les y ne se chevauchent pas"<<endl;
#endif
    return false;
  }
#ifdef DEBUG
  cout<<"les y se chevauchent"<<endl;
#endif
  return true;
}



bool PS_facet_3::has_a_cycle(PS_facet_3& Q){
  return ((_mark==true) && (Q.mark()==true));
}


bool PS_facet_3::two_test(PS_facet_3& Q){
   
  //A REVOIR
  if (is_back_side(Q)) {return true;}
  else if (Q.is_front_side(*this)) {return true;}
  else return false;
}

//Check if P(this) can be draw before Q
//Return true if one of the five test is true

bool PS_facet_3::five_test(PS_facet_3& Q){
   
  //A REVOIR
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
 
 if((a==0) && (c==0)) {_projection=ZX;
#ifdef DEBUG
 cout << "Plan de projection choisi : ZX" << endl;
#endif
}
 else if(v.z()==0) {_projection=YZ;
#ifdef DEBUG
cout << "Plan de projection choisi : YZ" << endl;
#endif
}
 else {_projection=XY;
#ifdef DEBUG
cout << "Plan de projection choisi : XY" << endl;
#endif
 }
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
 
#ifdef DEBUG
  cout << "\n############ ARRANGEMENT CREATION BEGIN ###########\n" << endl;
  cout << "Number of faces before arrangement creation " << (*this).arr.number_of_faces() << endl; 
  cout << "Number of vertices before arrangement creation : " << (*this).arr.number_of_vertices() << endl;
  cout << "The face is : " << (*this) << endl;
#endif

  check_projection();
 
  Arr_Point2 p1 = projection(l.first_point());
  Arr_Point2 q1 = projection(l.second_point());

#ifdef DEBUG  
  cout << "Inserted Point : " << p1 << " | " << q1 << endl;
#endif
  Curve *c1 = new Curve(p1,q1);

#ifdef DEBUG
  cout << "Adress of c1 : *c1 " << *c1 << " , &c1 " << &c1 << endl; 
#endif 

  Arr_2::Curve_iterator cit = arr.insert(*c1);
  pair<Curve_node *,PS_edge_3 *> pair1( &(*cit),&l);

#ifdef DEBUG
   cout << "Map Size " << mapp.size() << endl;
#endif
  
   mapp.insert(pair1);
   
#ifdef DEBUG 
   cout << "I insert in the Map : " <<  &(*cit) << " | " << &l << " | l= " << l << endl;
#endif
   
   Arr_Point2 initial = projection(Face_Arete[0].first_point());
   Arr_Point2 first = (*&initial);
   int i;
   for(i=0;i<_number_of_edge-1;i++) {

#ifdef DEBUG
     cout << "Etape : " << i << endl;
     cout << "Arete  " << i << " : " << Face_Arete[i].first_point() << " ";
     cout << Face_Arete[i].second_point() << endl; 
#endif

    Arr_Point2 second = projection(Face_Arete[i].second_point());

#ifdef DEBUG
    cout << "p2 = : " << first << " q2 = " << second << endl;
    cout << "Inserted Point : " << first << " | " << second << "\n" << endl;
#endif
    
    Curve *c2 = new Curve(first,second);
    cit = (*this).arr.insert(*c2);
    pair<Curve_node *,PS_edge_3 *> pair2(&(*cit),&(Face_Arete[i]));
    mapp.insert(pair2);
    
    first = second;
   }

#ifdef DEBUG
   cout << "Step : " << i << endl;
   cout << "Edge  " << i << " : " << Face_Arete[i].first_point() << " ";
   cout << Face_Arete[i].second_point() << endl;	       
#endif
 
   Curve *cf = new Curve(first,initial);

#ifdef DEBUG
    cout << "Inserted Point : " << first << " | " << initial << endl;
#endif

   cit = (*this).arr.insert(*cf);

#ifdef DEBUG
   cout << "TEST Last edge !!!! : " << &(*cit) << " | " << &(Face_Arete[_number_of_edge]) << endl; 
#endif
   
   pair<Curve_node *,PS_edge_3 *> pairf(&(*cit),&(Face_Arete[_number_of_edge-1]));
   mapp.insert(pairf); 

#ifdef DEBUG		  
  cout << "\nArrangement Display \n" << endl;
  for(cit=(*this).arr.curve_node_begin(); cit!=(*this).arr.curve_node_end(); ++cit){
    cout << "Curve level: " << cit->curve() << " | Pointeur : " << &(cit->curve()) << endl;
  }
  cit = (*this).arr.curve_node_begin();
  
  cout << "\nMap Display \n" << endl;
  std::map<Curve_node *,PS_edge_3 *>::iterator mi;
  cout << "Map size " << mapp.size() << "\n" <<endl;
  for (mi = mapp.begin();mi!=mapp.end();++mi) {
    cout << "Curve " << (*mi->first).curve() << " " << mi->first << endl;
    cout << "PS_edge  " << mi->second << " -> " << *(mi->second) << endl;
    cout << endl;
  }
  
  cout << "\nNumber of faces in the arrangement : " << arr.number_of_faces() << endl;
  cout << "Number of vertices in the arrangement : " << (*this).arr.number_of_vertices() << "\n" << endl;

  Arr_2::Vertex_iterator vit=(*this).arr.vertices_begin();
  int u=1;
  for(vit=(*this).arr.vertices_begin(); vit!=(*this).arr.vertices_end(); ++vit){
    cout << "Vertex " << u << " : " <<
      CGAL::to_double(vit->point().x()) << ","<<
      CGAL::to_double(vit->point().y()) << ","<< " addr : " << &*vit << endl;
    u++;
  }
  cout << "&(*this): " << &(*this) << endl;
  cout << "\n############ ARRANGEMENT CREATION END ###########\n" << endl;
#endif 
  if (arr.number_of_faces() >2) {return 1;}
  else return 0;
}

//To recover the facet of the arrangement
vector<PS_facet_3> PS_facet_3::create_facets(PS_edge_3 &l) {
  
  int nb_face=1;
  vector<PS_facet_3> facets;
  Arr_2::Face_iterator fit;

#ifdef DEBUG
  cout << "############ FACET CREATION BEGIN ###########\n" << endl;
  cout << "The face used is : " << *this << endl;
  cout << "Number of faces before facets creation : " << arr.number_of_faces() << endl;
  cout << "Number of vertices before facets creation : " << arr.number_of_vertices() << endl;
#endif 

  if (arr.number_of_faces() >2) {
    std::map<Curve_node *,PS_edge_3 *>::iterator mi;
    for(fit=arr.faces_begin();fit!=arr.faces_end();++fit) {
      
#ifdef DEBUG
      cout << "NEW FACE " << endl;
      cout << "FACE NUMBER " << nb_face << endl;
#endif  
      //the face is bounded and there is no hole inside the face
      if ( (!(*fit).is_unbounded()) && ( (*fit).holes_begin()== (*fit).holes_end()) ) {
	Arr_2::Ccb_halfedge_circulator initial=fit->outer_ccb(); 
	Arr_2::Ccb_halfedge_circulator e=fit->outer_ccb(); 
	//We have to check if the edge is not an edge create by the
	//cutting of the facet, because it is not possible to determine
	//if the facet is inner or outer.So if it is the case ,we take
	//another edge
	mi = mapp.find(&( * (e->edge_node()->curve_node())));

#ifdef DEBUG
     cout << "1St Time : I compare mi->second : " << mi->second << "and the cutting line :" << &l << endl;
#endif     

     while (mi->second == &l) {

#ifdef DEBUG
       cout << "I am in the WHILE" << endl;
       cout << "I have found that I compare mi->second : " <<
	 mi->second << " | and the cutting line :" << &l << endl;
#endif 
       e++;initial++;	
       mi = mapp.find(&( * (e->edge_node()->curve_node())));
     }
     
     //Here I am sure that the edge that I check is not provide by an
     //outer face or a hole,and moreover is not on the curve of the
     //cutting line.
#ifdef DEBUG      
     cout <<
       "\n#############################################################" << endl;
#endif
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

#ifdef DEBUG
     cout << "\nEdge : " ;
     cout << "Source : " << CGAL::to_double(source.x()) << " " << CGAL::to_double(source.y());
     cout << "| Target : " << CGAL::to_double(target.x()) << " " <<  CGAL::to_double(target.y()) << endl;
     cout << "mi->first : " << mi->first << endl;
     cout <<"*(mi->first) : " << (*mi->first).curve() << endl;
     cout << "mi->second : " << mi->second << endl; 
     cout << "*(mi->second) : " << *(mi->second) << "\n" << endl;    
     cout << "Je compare mi->second : " << mi->second << " | et la ligne de coupure :" << &l << endl;
#endif

     if( mi->second == &l) {
#ifdef DEBUG
       cout << "I have found that I compare mi->second : " <<
	 mi->second << " | and the cutting line : -> CUTTING LINE " << &l << endl;
       ///////A REVOIR
#endif
     
     }
     Vector2 v1(target_double-source_double);
     Vector2 v2(curve_target_double-curve_source_double);

#ifdef DEBUG
     cout << "Scalar Product between : " << endl;
     cout << "V1 " << v1 << " and " << endl;
     cout << "V2 " << v2 << endl;
     cout << "Result : v1*v2=" << v1*v2 << endl;
#endif
     
     //It is necessary to test if the halfedge and the PS_edge_3
     //corresponding got the same orientation. If not,the face is an
     //outer face and must be reject.To check this we make the scalar
     //product between the two vectors (the vector of the edge that we 
     //check and the curve_node of the arrangement which the edge
     //comes from).
     //If >0 the two edges got the same orientation so we keep the facet
     //else they don't have the same orientation so we reject the facet
     if (v1*v2>0) {
#ifdef DEBUG
       cout << "Same direction : the facet is keep\n" << endl;
#endif
       //It is necessary to reinitialize the edge iterator.if not It
       //miss the first edge use to know if the facet got a good orientation
       e=initial;
       std::vector<PS_edge_3> v_edge;	 
       
       do {
#ifdef DEBUG
	 cout << "I search : " << &( * (e->edge_node()->curve_node())) << endl;
#endif	 
	 mi = mapp.find(&( * (e->edge_node()->curve_node())));
       
#ifdef DEBUG
	 cout << "addr : " << &( * (e->edge_node()->curve_node())) << endl;  
	 cout << "PS_edge corresponding : 1-2 " << mi->first << "   " << 
	   mi->second  << "  " << *(mi->second) << endl;
	 /////////////////////////////////////////
	 //Edge creation
	 cout << "*e->source() " << &(*(e->source())) << " | *e->target() "
	      << &(*(e->target())) << endl;
#endif

	 Arr_Point2 p1 = (*(e->source())).point();
	 Arr_Point2 p2 = (*(e->target())).point();
	 Point3 p1_3D = lift(p1);
	 Point3 p2_3D = lift(p2);

#ifdef DEBUG	 
	 cout << "p1=" << CGAL::to_double(p1.x()) << "," <<
	   CGAL::to_double(p1.y()) << " |p2=" <<
	   CGAL::to_double(p2.x()) << "," << CGAL::to_double(p2.y());
	 
	 cout << "\nlift : " << lift(p1) << "   " << lift(p2) << endl;
	 cout << "EDGE Display : " << (*(mi->second)) << endl;
	 cout << "Color : " << (*(mi->second)).color();
	 cout << " | VISIBILITY : " << (*(mi->second)).visibility() << endl;
#endif

	 PS_edge_3 edge(p1_3D,p2_3D, (*(mi->second)).color(),(*(mi->second)).visibility());

#ifdef DEBUG
	 cout << "Creation of the edge : " << edge << "\n" << endl;
#endif
	 v_edge.push_back(edge);
	 //////////////////////////////////////
	 e++;
       }
     while (e != initial);
       //When the new facet is create we give it the number of the
       //mother facet.
       PS_facet_3 f(v_edge,color(),filling(),number_of_the_facet());
       facets.push_back(f);
       cout << endl;
     }
     else {
#ifdef DEBUG
       cout << "The face is rejected : negative scalar product !!! \n" << endl;
#endif
     }
      }
      else {
#ifdef DEBUG
	cout << "This face is an hole or is an outer face !!! \n " << endl;
#endif
      }
      nb_face++;
    }    
#ifdef DEBUG 
    cout << "############ FACET CREATION END ###########\n" << endl; 
#endif 
    return facets;
  }
  else {
#ifdef DEBUG
    cout << "Number of face in the arrangement < 2 " << endl;
#endif
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
#ifdef DEBUG
	  cout << "The segment cut themselves REALLY !!!" << endl;
#endif
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
  
  //cout << "JE SUIS APT TEST DU UNBOUNDED " << endl;
  //cout << "du coup je renvoit faux "<< endl;

#ifdef DEBUG
 cout << "The xy don't overlap !!! " << endl;
#endif
 
 return false;
}

void PS_facet_3::transformation(Transformation &t) {
  for (int i=0;i<_number_of_edge;i++) {
    Face_Arete[i].transformation(t);
  }
}

CGAL_END_NAMESPACE

#endif
