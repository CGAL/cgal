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
// file          : src/PS_Stream_3.C
// package       : PS_Stream
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : 
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================
 
#include <CGAL/basic.h>
#include <CGAL/IO/PS_Stream_3.h>

CGAL_BEGIN_NAMESPACE

//Return the maximum in x,y,z of the BBox3
coord_type PS_Stream_3::search_max(PS_BBox3& b) {
  
  coord_type zmax= b.zmax();
  coord_type ymax= b.ymax();
  coord_type xmax= b.xmax();
 
  vector<coord_type> v;
  v.push_back(xmax);
  v.push_back(ymax);
  v.push_back(zmax);

  sort(v.begin(),v.end());

  return v[2];
}
//Return the minimum in x,y,z of the BBox3
coord_type PS_Stream_3::search_min(PS_BBox3& b) {
  coord_type zmin = b.zmin();
  coord_type ymin = b.ymin();
  coord_type xmin = b.xmin();
 
  vector<coord_type> v;
  v.push_back(xmin);
  v.push_back(ymin);
  v.push_back(zmin);

  sort(v.begin(),v.end());
  return (v[0]);
} 


PS_Stream_3::PS_Stream_3(const PS_BBox3& bb3,const Direction& d,
			 const Direction &l,ostream& os, OutputMode mode):
 PS_Stream(os,mode),_bbox3(bb3),_point_of_view(d),_light(l)
{
double xn=d.dx()/norme(d.dx(),d.dy(),d.dz());
double yn=d.dy()/norme(d.dx(),d.dy(),d.dz());
double zn=d.dz()/norme(d.dx(),d.dy(),d.dz());

_transformation=Transformation (yn/den(zn),-xn/den(zn),0,zn*xn/den(zn),zn*yn/den(zn),-den(zn),xn,yn,zn);

init();

for ( int i=0;i<=7;i++)
  transform(_transformation,v[i]);

_bbox=PS_BBox(CGAL::to_double(search_xmin(v)),CGAL::to_double(search_ymin(v)),CGAL::to_double(search_xmax(v)),CGAL::to_double(search_ymax(v)));
 set_scale(_bbox);

}

PS_Stream_3::PS_Stream_3(const PS_BBox3& bb3, const Direction& d, const Direction &l,const
			 char* fname, OutputMode mode):
 PS_Stream(fname,mode),_bbox3(bb3),_point_of_view(d),_light(l)

{
double xn=d.dx()/norme(d.dx(),d.dy(),d.dz());
double yn=d.dy()/norme(d.dx(),d.dy(),d.dz());
double zn=d.dz()/norme(d.dx(),d.dy(),d.dz());

_transformation=Transformation (yn/den(zn),-xn/den(zn),0,zn*xn/den(zn),zn*yn/den(zn),-den(zn),xn,yn,zn);
init();
for ( int i=0;i<=7;i++)
  transform(_transformation,v[i]);
_bbox=PS_BBox(CGAL::to_double(search_xmin(v)),CGAL::to_double(search_ymin(v)),CGAL::to_double(search_xmax(v)),CGAL::to_double(search_ymax(v)));

 set_scale(_bbox);

}

PS_Stream_3::PS_Stream_3(const PS_BBox3& bb3, const Direction& d,
			 const Direction &l,float H, const
	    char* fname, OutputMode mode):
 PS_Stream(H,fname,mode),_bbox3(bb3),_point_of_view(d),_light(l)
{

double xn=d.dx()/norme(d.dx(),d.dy(),d.dz());
double yn=d.dy()/norme(d.dx(),d.dy(),d.dz());
double zn=d.dz()/norme(d.dx(),d.dy(),d.dz());
double deno=sqrt(xn*xn+yn*yn);

 if (deno != 0) {
_transformation=Transformation (yn/deno,-xn/deno,0,zn*xn/deno,zn*yn/deno,-deno,xn,yn,zn);
 }
 else{
   if (zn>0) {
     _transformation=Transformation (zn,0,0,0,zn,0,0,0,zn);}
   else  _transformation=Transformation (zn,0,0,0,-zn,0,0,0,zn);
 

 }

init();

for ( int i=0;i<=7;i++)
  transform(_transformation,v[i]);

_bbox=PS_BBox(CGAL::to_double(search_xmin(v)),CGAL::to_double(search_ymin(v)),CGAL::to_double(search_xmax(v)),CGAL::to_double(search_ymax(v)));
set_window(_bbox,H); 
set_scale(_bbox);
insert_catalogue(); 

 _current_filling=CGAL::NO_FILL;
_number_of_face=0;

}

PS_Stream_3::PS_Stream_3(const PS_BBox3& bb3,const Direction& d, const Direction &l,float H, ostream& os,
	    OutputMode mode ) :
 PS_Stream(H,os,mode),_bbox3(bb3),_point_of_view(d),_light(l)
{
double xn=d.dx()/norme(d.dx(),d.dy(),d.dz());
double yn=d.dy()/norme(d.dx(),d.dy(),d.dz());
double zn=d.dz()/norme(d.dx(),d.dy(),d.dz());

_transformation=Transformation (yn/den(zn),-xn/den(zn),0,zn*xn/den(zn),zn*yn/den(zn),-den(zn),xn,yn,zn);

init();

for ( int i=0;i<=7;i++)
  transform(_transformation,v[i]);

_bbox=PS_BBox(CGAL::to_double(search_xmin(v)),CGAL::to_double(search_ymin(v)),CGAL::to_double(search_xmax(v)),CGAL::to_double(search_ymax(v)));
set_window(_bbox,H); 
set_scale(_bbox);
insert_catalogue();
}
/////////////////////////////// END OF THE CONSTRUCTORS ////////////////////////////////////////////////




Point2  PS_Stream_3::transform(Transformation  t,const Point3& p)
{
  Point3 r(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
  Point3  q=t(r);
  return Point2 (q.x(),q.y());
}

coord_type PS_Stream_3::search_xmin(Point3  v[])
{
coord_type xmin=v[0].x();
for(int i=1;i<8;i++)
  {if (xmin>v[i].x())
    xmin=v[i].x();
  }
return xmin;
}

coord_type PS_Stream_3::search_ymin(Point3  v[])
{
coord_type ymin=v[0].y();
for(int i=1;i<8;i++)
  {if (ymin>v[i].y())
    ymin=v[i].y();
  }
return ymin;
}

coord_type PS_Stream_3::search_xmax(Point3  v[])
{
coord_type xmax=v[0].x();
for(int i=1;i<8;i++)
  {if (xmax<v[i].x())
    xmax=v[i].x();
  }
return xmax;
}

coord_type PS_Stream_3::search_ymax(Point3  v[])
{
coord_type ymax=v[0].y();
for(int i=1;i<8;i++)
  {if (ymax<v[i].y())
    ymax=v[i].y();
  }
return ymax;
}

Point2 PS_Stream_3::transform( const Point3 & p){
Point2  q(p.x(),p.y());
return q;
}

Point_2<D> PS_Stream_3::convert_to_double(Point2& p) {
  Point_2<D> q(CGAL::to_double(p.x()),CGAL::to_double(p.y()));
  return q;
}


void PS_Stream_3::transformation() {

  for(int i=0;i<_number_of_face;i++) {
    _scene[i].transformation(_transformation);
    _scene[i].check_projection();
  } 
  //After the transformation the point of view become 0,0,1 !! And the 
  //light had to be transform also.
  _point_of_view=_transformation(_point_of_view); 
  _light=_transformation(_light);
}


//Make the sort according the increasing z
vector<PS_facet_3> PS_Stream_3::zmin_sort() {

  for(int i=_number_of_face;i>0;i--) {
    for(int j=1;j<i;j++) {
      if ((_scene[j-1].zmin()) > (_scene[j].zmin())) {
	PS_facet_3 temp(_scene[j-1]);
	_scene[j-1]=_scene[j];
	_scene[j]=temp;
      }
      else if ((_scene[j-1].zmin()) == (_scene[j].zmin())) {
	if ((_scene[j-1].zmax()) < (_scene[j].zmax())) {
	  PS_facet_3 temp(_scene[j-1]);
	  _scene[j-1]=_scene[j];
	  _scene[j]=temp; 
	}
	else {
	  PS_facet_3 temp(_scene[j]);
	  _scene[j]=_scene[j-1];
	  _scene[j-1]=temp; 
	}
      }
    }
  } 
  return _scene;
}


//make the normal sort to determine visible surfaces
vector<PS_facet_3> PS_Stream_3::normal_sort() {

Direction l=_light;
Direction pov=_point_of_view;

//Vector with only the visible PS_facet_3  
vector<PS_facet_3> scene_visible;
  
 for(int i=0;i<_number_of_face;i++) {

   Plane3 P1(_scene[i][0].first_point(),_scene[i][0].second_point(),_scene[i][_scene[i].number_of_edge()-1].first_point());
   //Creation of the normal vector of the plane
   Vector3 v_ortho(P1.orthogonal_vector());
   Vector3
     normal_normalise(v_ortho.x(),v_ortho.y(),v_ortho.z(),norme(v_ortho.x(),v_ortho.y(),v_ortho.z()));
   //Creation of the direction of pov vector  
   Vector3 v_pov(pov.vector());
   Vector3 pov_normalise(v_pov.x(),v_pov.y(),v_pov.z(),norme(v_pov.x(),v_pov.y(),v_pov.z()));

   //Creation of the direction of light vector  
   Vector3 v_light(l.vector());
   Vector3 light_normalise(v_light.x(),v_light.y(),v_light.z(),norme(v_light.x(),v_light.y(),v_light.z()));
   
   double angle = CGAL::to_double(normal_normalise*pov_normalise);
   double angle_light = light_normalise*normal_normalise;
   //The face is visible if angle>0
   if (angle>0) {//The facet is visible
     _scene[i].set_gray_level(angle_light);
     scene_visible.push_back(_scene[i]);
   }
 }
 _scene=scene_visible;
 _number_of_face=scene_visible.size();

 return scene_visible;
}

void PS_Stream_3::display() {

  transformation(); 
  normal_sort();
  zmin_sort();
  depth_sort();

  for (int i=0;i<_number_of_face;i++) {
    (*this) << _scene[i];
  }
}


int PS_Stream_3::cutting(vector<PS_facet_3>& liste_face,int i,int j) {

  //The index i give the facet which is going to be cut and the j the face which intersect facet i
  Line3 l;Plane3 p;CGAL::Object obj;  
  Plane3 P1(liste_face[i][0].first_point(),liste_face[i][0].second_point(),liste_face[i][liste_face[i].number_of_edge()-1].first_point());
 Plane3
   P2(liste_face[j][0].first_point(),liste_face[j][0].second_point(),liste_face[j][liste_face[j].number_of_edge()-1].first_point());

 //POSSIBLE PROBLEM IF THE INTERSECTION IS A PLAN
 //REDONE THE CREATION OF THE INTERSECT LINE(NOT VERY SURE !!!)
 obj = intersection(P1,P2);
 assign(l, obj); 

 int valmin = (int)(CGAL::to_double(liste_face[i].min()) -1);
 int valmax = (int)(CGAL::to_double(liste_face[i].max()) +1);
 Point3 la = l.point(valmin); 
 Point3 lb = l.point(valmax);

 CGAL::PS_edge_3 line_edge(la,lb,CGAL::GREEN,CGAL::INVISIBLE);
 
 vector<CGAL::PS_facet_3> facet1;
 if (liste_face[i].make_arrangement(line_edge)) {
   facet1 = liste_face[i].create_facets(line_edge);
   liste_face.erase(&liste_face[i]);
  
   for(int k=0;k<(int)liste_face.size();k++) {
     facet1.push_back(liste_face[k]);
   }
   liste_face=facet1;
   return 1;
 }
 else {
   facet1.push_back(liste_face[i]);
   liste_face.erase(&liste_face[i]);
   
   for(int k=0;k<(int)liste_face.size();k++) {
     facet1.push_back(liste_face[k]);
   }

   liste_face=facet1;
   return 0;
 }
}

//Check if the face i is under the face j
int PS_Stream_3::dessous(PS_facet_3& facei,PS_facet_3& facej ) {
  
  Plane3 Pi(facei[0].first_point(),facei[0].second_point(),facei[facei.number_of_edge()-1].first_point());
  double a= Pi.a();double b= Pi.b();double c=Pi.c();double d=Pi.d();

  vector<double> v;
  double result;
  for (int i=0;i<facej.number_of_edge();i++) {
    result = a*facej[i].first_point().x()+b*facej[i].first_point().y()+c*facej[i].first_point().z()+d;
    v.push_back(result);
  }
  //Sort by increasing values
  sort(v.begin(),v.end());
  
  int indice=0;
  for(int i=0;i<(int)v.size();i++) {
    if (CGAL_NTS abs(v[indice]) < CGAL_NTS abs(v[i])) {indice=i;}
  }
  if (v[indice]>=0) return 1;
  else return (-1);
}


//Precondition at least 2 faces !!!!
void PS_Stream_3::depth_sort() {

  if (_scene.size() > 1) {
    vector<PS_facet_3> liste_face(_scene);
    vector<PS_facet_3> liste_face_triee;
    int i=0;int j=1;

    while( i<(int)liste_face.size()-1) {
      i=0;j=1;
      while(j<(int)liste_face.size()) { 
	if (liste_face[i].number_of_the_facet() == liste_face[j].number_of_the_facet()) {
	  //Coplanar faces so I don't cut , I try with the next facet.
	  j++;}
	else {
	  if (liste_face[i].five_test(liste_face[j])) {j++;}
	  else {
	    if (liste_face[j].has_a_cycle(liste_face[i])) {
	      cutting(liste_face,j,i);
	      i=0;j=1;
	    }
	    else {
	      if(liste_face[j].two_test(liste_face[i])) {//Swap P and Q
		PS_facet_3 temp=liste_face[i];
		liste_face[i]=liste_face[j];
		liste_face[j]=temp;
		liste_face[i].set_mark(true);
		i=0;j=1;
	      }
	      else {  
		if (cutting(liste_face,i,j)==0) {
		  if (dessous(liste_face[i],liste_face[j])) {j++;}
		  else  {//Swap P and Q
		    PS_facet_3 temp=liste_face[i];
		    liste_face[i]=liste_face[j];
		    liste_face[j]=temp;
		    liste_face[i].set_mark(true);
		    i=0;j=1;
		  }
		}
		else {i=0;j=1;}
	      }
	    }
	  }
	}
      }
      liste_face_triee.push_back(liste_face[i]);
      liste_face.erase(&liste_face[i]);
    }
    //Add the last element
    liste_face_triee.push_back(liste_face[i]);
    liste_face.erase(&liste_face[i]);
    _scene=liste_face_triee; 
    _number_of_face=liste_face_triee.size();
  }
}

PS_Stream_3& operator<< (PS_Stream_3& ps, const PS_edge_3& ar) {

if (ar.visibility() == VISIBLE) {

  if (ps.is_readable()) {
  ps.os() << "\n%***************************************************%" << endl;
  ps.os() << "%CGAL% - PS_EDGE_3  " << endl;
  ps.os() << "%CGAL% - [(" << ar.first_point() << ") ,(" <<
    ar.second_point() << ")]" << endl;
  }
  //We regain the color of the edge
  ps.os() << ar.color() << " setrgbcolor" << endl;  

  //Projection of the points
  Point2  pa = ps.transform(ar.first_point());
  Point2  pb = ps.transform(ar.second_point());
  
  //Transformation to pass in postscript
  ps.os() << ps.x2ps(CGAL::to_double(pa.x())) << " " << ps.y2ps(CGAL::to_double(pa.y())) << " " << ps.x2ps(CGAL::to_double(pb.x()))<< "  "<< ps.y2ps(CGAL::to_double(pb.y())) << " arete" << endl;
  }
 if (ps.is_readable()) {
  ps.os() << "\n%***************************************************%" << endl;
 }
return ps;
}

PS_Stream_3& operator << (PS_Stream_3& ps, const PS_facet_3& face) {

  if(ps.is_readable()){
  ps.os() << "\n%***************************************************%" << endl;
  ps.os() << "%CGAL% - PS_FACET_3 " << endl;
  ps.os() << "%CGAL% - Number of edges : " << face.number_of_edge() << " " << endl;  
  }
  switch (face.filling()) {
    
  case NO_FILL :
    { 
      ps.os() << endl;
      break;
    }
    
  case WIRED_CULLBACK_FACING :
    { 
      if(ps.is_readable()){
	ps.os() << "%CGAL% - WIRED_CULLBACK_FACING" << endl;  
	ps.os() <<
	  "%***************************************************%" <<
	  endl;
      }
      //Put the color of the backgroung (white by default)
      //Need to be modified if the background color change
      ps.os() << "1 1 1 setrgbcolor" << endl;  
   
      for (int i=0;i<face.number_of_edge();i++) {
	Point2  pb = ps.transform(face[i].second_point());
	//Transformation to pass in postscript coordinates
	ps.os () << ps.x2ps(CGAL::to_double(pb.x()))<< " " << ps.y2ps(CGAL::to_double(pb.y())) << " ";
      }
      ps.os() << face.number_of_edge() << " face" << endl;
      ps.os() << "fill" << endl;
      ps.os() << endl;
      break;
    }

  case UNIFORM_FILL :
    {
       if(ps.is_readable()){
	 ps.os() << "%CGAL% - UNIFORM FILL - Color :" << face.color() << endl;  
	 ps.os() <<
	   "%***************************************************%" <<
	   endl;
       }
       //We got the color of the facet
       ps.os() << face.color() << " setrgbcolor" << endl;  
       
       for (int i=0;i<face.number_of_edge();i++) {
	 Point2  pb = ps.transform(face[i].second_point());
	 //Transformation to pass in postscript coordinates
	 ps.os () << ps.x2ps(CGAL::to_double(pb.x()))<< " " << ps.y2ps(CGAL::to_double(pb.y())) << " ";
       }
       ps.os() << face.number_of_edge() << " face" << endl;
       ps.os() << "fill" << endl;
       ps.os() << endl;
       break;
    }
    
  case NORMAL_FILL :
    {
       if(ps.is_readable()){
      ps.os() << "%CGAL% - NORMAL FILL " << endl; 
      ps.os() <<
	"%***************************************************%" << endl;
       }
      for (int i=0;i<face.number_of_edge();i++) {
	Point2  pb = ps.transform(face[i].second_point());
	ps.os () << ps.x2ps(CGAL::to_double(pb.x())) << "  "<< ps.y2ps(CGAL::to_double(pb.y())) << " " << endl;
      }
      ps.os() << face.number_of_edge() << " face" << endl;
      ps.os() << face.gray_level() <<" setgray" <<endl;			    
      ps.os() << "fill" <<endl;
      ps.os() << endl;
      break;
    }

    
 default :

    if(ps.is_readable()){
      ps.os() << "%***************************************************%\n\n" << endl;
    }
   return ps;
} 

//In all the cases we paint the edges with their color
  for (int i=0;i<face.number_of_edge();i++) {
    ps << face[i];
  }
   if(ps.is_readable()){
     ps.os() << "%****************************************************%\n\n" << endl;
   }
  return ps;
}



//Convert a triangle in PS_facet_3
PS_facet_3 PS_Stream_3::triangle_to_facet(Triangle3& T) {
  
  Color edge_color = context().get_border_color();
  Color filling_color = context().get_fill_color();
  FILLING filling_type = _current_filling;

  Point3 a(T[0].x(),T[0].y(),T[0].z());
  Point3 b(T[1].x(),T[1].y(),T[1].z());
  Point3 c(T[2].x(),T[2].y(),T[2].z());
  
  vector<Point3> triangle;
  triangle.push_back(a);
  triangle.push_back(b);
  triangle.push_back(c);
  
  PS_facet_3 face(triangle,edge_color,filling_color,filling_type,_scene_initiale.size()+1);
  return face;
}

//Allows to add a triangle 3D in the PS_Stream_3
void PS_Stream_3::add_triangle(Triangle3 &triangle) {
PS_facet_3 tri = triangle_to_facet(triangle);
add_facet(tri);
}



PS_Stream_3& operator << (PS_Stream_3& ps, Triangle3& t) {
ps.add_triangle(t);
return ps;
}


//Convert a tetrahedron in a vector of PS_facet_3
vector<PS_facet_3> PS_Stream_3::tetrahedron_to_facet(Tetrahedron& T1) {
  
Color edge_color = context().get_border_color();
Color filling_color = context().get_fill_color();
FILLING filling_type = _current_filling;

   Tetrahedron T=T1; 

   Point3 a(T[0].x(),T[0].y(),T[0].z());
   Point3 b(T[1].x(),T[1].y(),T[1].z());
   Point3 c(T[2].x(),T[2].y(),T[2].z());
   Point3 d(T[3].x(),T[3].y(),T[3].z());

   vector<Point3> f1;
   f1.push_back(a);f1.push_back(b);f1.push_back(c);
   PS_facet_3 face1(f1,edge_color,filling_color,filling_type,_scene_initiale.size()+1);
   
   vector<Point3> f2;
   f2.push_back(d);f2.push_back(a);f2.push_back(c);
   PS_facet_3 face2(f2,edge_color,filling_color,filling_type,_scene_initiale.size()+1);

   vector<Point3> f3;
   f3.push_back(d);f3.push_back(c);f3.push_back(b);
   PS_facet_3 face3(f3,edge_color,filling_color,filling_type,_scene_initiale.size()+1);
   
   vector<Point3> f4;
   f4.push_back(a);f4.push_back(d);f4.push_back(b);
   PS_facet_3 face4(f4,edge_color,filling_color,filling_type,_scene_initiale.size()+1);

   vector<PS_facet_3> tetrahedre;
   tetrahedre.push_back(face1); tetrahedre.push_back(face2);
   tetrahedre.push_back(face3); tetrahedre.push_back(face4);
   return tetrahedre;
}

//Allows to add a tetrahedron in the PS_Stream_3
void PS_Stream_3::add_tetrahedron(Tetrahedron &t) {
vector<PS_facet_3> tetrahedron = tetrahedron_to_facet(t);
add_vfacet(tetrahedron);
}


PS_Stream_3& operator << (PS_Stream_3& ps, Tetrahedron& t) {
ps.add_tetrahedron(t);
return ps;
}

//Allow to add a point 3D in the scene (to begin approximate by a cube)
vector<PS_facet_3> PS_Stream_3::point_to_facet(Point3& point){
  
  vector<PS_facet_3> cube;
  unsigned int mid_size = context().get_dot_size()/2;
  Color edge_color = context().get_border_color();
  Color filling_color = context().get_fill_color();
  FILLING filling_type = _current_filling;
 
  Point3 p1(point.x() - mid_size,point.y() - mid_size,point.z() +
	    mid_size);
  Point3 p2(point.x() + mid_size,point.y() - mid_size,point.z() +
	    mid_size);
  Point3 p3(point.x() + mid_size,point.y() + mid_size,point.z() +
	    mid_size);
  Point3 p4(point.x() - mid_size,point.y() + mid_size,point.z() +
	    mid_size); 
  Point3 p5(point.x() - mid_size,point.y() - mid_size,point.z() -
	    mid_size);
  Point3 p6(point.x() + mid_size,point.y() - mid_size,point.z() -
	    mid_size);
  Point3 p7(point.x() + mid_size,point.y() + mid_size,point.z() -
	    mid_size);
  Point3 p8(point.x() - mid_size,point.y() + mid_size,point.z() -
	    mid_size); 
  
  vector<Point3> face1;
  face1.push_back(p1);face1.push_back(p2);face1.push_back(p3);face1.push_back(p4);
  vector<Point3> face2;
  face2.push_back(p2);face2.push_back(p6);face2.push_back(p7);face2.push_back(p3);
  vector<Point3> face3;
  face3.push_back(p6);face3.push_back(p5);face3.push_back(p8);face3.push_back(p7);
  vector<Point3> face4;
  face4.push_back(p5);face4.push_back(p1);face4.push_back(p4);face4.push_back(p7);
  vector<Point3> face5;
  face5.push_back(p5);face5.push_back(p6);face5.push_back(p2);face5.push_back(p1);
  vector<Point3> face6;
  face6.push_back(p4);face6.push_back(p3);face6.push_back(p7);face6.push_back(p8);

  PS_facet_3 facet1(face1,edge_color,filling_color,filling_type,_scene_initiale.size()+1);
  PS_facet_3 facet2(face2,edge_color,filling_color,filling_type,_scene_initiale.size()+2);
  PS_facet_3 facet3(face3,edge_color,filling_color,filling_type,_scene_initiale.size()+3);
  PS_facet_3 facet4(face4,edge_color,filling_color,filling_type,_scene_initiale.size()+4);
  PS_facet_3 facet5(face5,edge_color,filling_color,filling_type,_scene_initiale.size()+5);
  PS_facet_3 facet6(face6,edge_color,filling_color,filling_type,_scene_initiale.size()+6);

  cube.push_back(facet1);cube.push_back(facet2);cube.push_back(facet3);cube.push_back(facet4);
  cube.push_back(facet5);cube.push_back(facet6);

return cube;
}

void PS_Stream_3::add_point(Point3 &point){

vector<PS_facet_3> cube = point_to_facet(point);
add_vfacet(cube);
}

PS_Stream_3& operator <<(PS_Stream_3& ps,  Point3& p) {

ps.add_point(p);
return ps;
}





vector<PS_facet_3> PS_Stream_3::segment_to_facet(Segment3& segment){
  
  vector<PS_facet_3> para;
  unsigned int mid_size = context().get_thickness()/2;
  mid_size=15/2;

  Color edge_color = context().get_border_color();
  Color filling_color = context().get_fill_color();
  FILLING filling_type = _current_filling;
  Point3 source = segment.source();
  Point3 target = segment.target();
 
  Point3 p1(source.x() - mid_size,source.y() - mid_size,source.z());
  Point3 p2(source.x() + mid_size,source.y() - mid_size,source.z());
  Point3 p3(source.x() + mid_size,source.y() + mid_size,source.z());
  Point3 p4(source.x() - mid_size,source.y() + mid_size,source.z());
 
  Point3 p5(target.x() - mid_size,target.y() - mid_size,target.z());
  Point3 p6(target.x() + mid_size,target.y() - mid_size,target.z());
  Point3 p7(target.x() + mid_size,target.y() + mid_size,target.z());
  Point3 p8(target.x() - mid_size,target.y() + mid_size,target.z()); 
 
  vector<Point3> face1;
  face1.push_back(p1);face1.push_back(p2);face1.push_back(p3);face1.push_back(p4);
  vector<Point3> face2;
  face2.push_back(p2);face2.push_back(p6);face2.push_back(p7);face2.push_back(p3);
   vector<Point3> face3;
   face3.push_back(p6);face3.push_back(p5);face3.push_back(p8);face3.push_back(p7);
  vector<Point3> face4;
  face4.push_back(p5);face4.push_back(p1);face4.push_back(p4);face4.push_back(p7);
  vector<Point3> face5;
  face5.push_back(p5);face5.push_back(p6);face5.push_back(p2);face5.push_back(p1);
     vector<Point3> face6;
   face6.push_back(p4);face6.push_back(p3);face6.push_back(p7);face6.push_back(p8);

  PS_facet_3 facet1(face1,edge_color,filling_color,filling_type,_scene_initiale.size()+1);
  PS_facet_3 facet2(face2,edge_color,filling_color,filling_type,_scene_initiale.size()+2);
  PS_facet_3 facet3(face3,edge_color,filling_color,filling_type,_scene_initiale.size()+3);
  PS_facet_3 facet4(face4,edge_color,filling_color,filling_type,_scene_initiale.size()+4);
  PS_facet_3 facet5(face5,edge_color,filling_color,filling_type,_scene_initiale.size()+5);
  PS_facet_3 facet6(face6,edge_color,filling_color,filling_type,_scene_initiale.size()+6);
  
  para.push_back(facet1);para.push_back(facet2);
  para.push_back(facet3);para.push_back(facet4);
  para.push_back(facet5);
  para.push_back(facet6);

return para;
}

void PS_Stream_3::add_segment(Segment3 &segment){

vector<PS_facet_3> para = segment_to_facet(segment);
add_vfacet(para);
}

//Allows to add a Segment_3 in the PS_Stream_3
PS_Stream_3& operator << (PS_Stream_3& ps,  Segment3& 
			  segment) {

ps.add_segment(segment);
return ps;
}

vector<PS_facet_3> PS_Stream_3::line_to_facet(Line3& l){
//UTILISER SEGMENT AVEC UNE TAILLE SUFFISSANTE QUI SORT DE LA BBOX

Point3 initiale = l.point(0);
Direction d = l.direction();
Vector3 vect = d.vector();

PS_BBox3 bbox3 = _bbox3;

Point3 source = initiale + search_min(bbox3)*vect;
Point3 destination = initiale + search_max( bbox3)*vect;

// For the line I use the segment with two points which go out the BBox 
Segment3 seg(source,destination);
vector<PS_facet_3> line = segment_to_facet(seg);
return line;
}

void PS_Stream_3::add_line(Line3 &l){

vector<PS_facet_3> line = line_to_facet(l);
add_vfacet(line);
}

//Allows to add a new Line_3 in the PS_Stream_3 with the operator <<
PS_Stream_3& operator << (PS_Stream_3& ps, Line3& 
			  line) {

ps.add_line(line);
return ps;
}



PS_Stream_3& operator << (PS_Stream_3& ps, FILLING& filling) {
ps.set_current_filling(filling);
return ps;
}




template <class Traits, class HDS>
vector<CGAL::PS_facet_3> PS_Stream_3::polyhedron_to_facet(Polyhedron_3<Traits,HDS>& p) {

  typedef Polyhedron_3<Traits,HDS>::Facet_const_iterator FCI;
  typedef Polyhedron_3<Traits,HDS>::Halfedge_around_facet_const_circulator HFCC;
  typedef Polyhedron_3<Traits,HDS>::Point Point;
  
  Color edge_color = context().get_border_color();
  Color filling_color = context().get_fill_color();
  FILLING filling_type = _current_filling;
  vector<CGAL::PS_facet_3> polyh;
  
  FCI fci;
  for(fci = p.facets_begin(); fci != p.facets_end();++fci)
    {
      vector<Point3> polygone;
      HFCC hfcc = fci->facet_begin();
      HFCC hfcc_end = hfcc;
      do {
	polygone.push_back((Point3)(hfcc->vertex())->point());
	++hfcc;
      }
      while(hfcc != hfcc_end);
      
      CGAL::PS_facet_3 facet(polygone,edge_color,filling_color,filling_type,_scene_initiale.size()+1);
      polyh.push_back(facet); 
    }
  return polyh;
} 

template <class Traits, class HDS>
void PS_Stream_3::add_polyhedron(Polyhedron_3<Traits,HDS>& p){
vector<PS_facet_3> polyhedron = polyhedron_to_facet(p);
add_vfacet(polyhedron);
}

template <class Traits, class HDS>
PS_Stream_3& operator <<(PS_Stream_3& ps,Polyhedron_3<Traits,HDS>& p) {
ps.add_polyhedron(p);
return ps;
}


//Allows to add a new face in the scene
void PS_Stream_3::add_facet(PS_facet_3 &face){
  _scene.push_back(face);
  _scene_initiale.push_back(face);
  _number_of_face++;
}

//Allows to add a vector of facet in the scene
void PS_Stream_3::add_vfacet(vector<PS_facet_3> &vfacet){

  for (int i=0;i<(int)vfacet.size();i++) {
  _scene.push_back(vfacet[i]);
  _scene_initiale.push_back(vfacet[i]);
  _number_of_face++;
  }
}

template <class triangulation_2>
vector<PS_facet_3> PS_Stream_3::triangulation_2_to_facet(triangulation_2& tr) {

  typedef typename triangulation_2::Face_iterator Face_iterator;
  typedef typename triangulation_2::Point Point_tr;

  //typedef triangulation_2::Point Point_tr;
  Face_iterator fit = tr.faces_begin();
  vector<CGAL::PS_facet_3> vfacet;
  Color edge_color = context().get_border_color();
  Color filling_color = context().get_fill_color();
  FILLING filling_type = _current_filling;
  int nbface=0;
  
  for (fit =  tr.faces_begin(); fit !=  tr.faces_end(); fit ++){
    nbface++;
    Point_tr p0 = (fit)->vertex(0)->point();
    Point_tr p1 = (fit)->vertex(1)->point();
    Point_tr p2 = (fit)->vertex(2)->point();

    Point3 p0_3D(p0.x(),p0.y(),p0.z());
    Point3 p1_3D(p1.x(),p1.y(),p1.z());
    Point3 p2_3D(p2.x(),p2.y(),p2.z());
    
    vector<Point3> vpoint1;
    vpoint1.push_back(p2_3D);vpoint1.push_back(p1_3D);vpoint1.push_back(p0_3D);  
    vector<Point3> vpoint2;
    vpoint2.push_back(p0_3D);vpoint2.push_back(p1_3D);vpoint2.push_back(p2_3D);
      
    PS_facet_3 facet1(vpoint1,edge_color,filling_color,filling_type,_scene_initiale.size()+1);
    PS_facet_3 facet2(vpoint2,edge_color,filling_color,filling_type,_scene_initiale.size()+1);
   
    vfacet.push_back(facet1); vfacet.push_back(facet2);
  }
return vfacet;
}

template <class triangulation_2>
void PS_Stream_3::add_triangulation(triangulation_2& tr){

vector<PS_facet_3> triangulation = triangulation_2_to_facet(tr);
add_vfacet(triangulation);
}


template <class Gt,class Tds>
PS_Stream_3& operator << (PS_Stream_3& ps, Triangulation_2<Gt,Tds> &t){
  ps.add_triangulation(t);
  return ps;
}

template < class Gt, class Tds >
PS_Stream_3& operator << (PS_Stream_3& ps, Delaunay_triangulation_2<Gt,Tds> &t){
  ps.add_triangulation(t); 
  return ps; 
}

template < class Gt, class Tds>
PS_Stream_3& operator<<(PS_Stream_3& ps, Constrained_triangulation_2<Gt,Tds> &t){
  ps.add_triangulation(t);  
  return ps;
}

template < class Gt, class Tds >
PS_Stream_3& operator << (PS_Stream_3& ps, Regular_triangulation_2<Gt,Tds> &t){
  ps.add_triangulation(t);
  return ps;
}

template < class Gt, class Tds >
PS_Stream_3& operator << (PS_Stream_3& ps, Constrained_Delaunay_triangulation_2<Gt,Tds> &t){
 ps.add_triangulation(t);
 return ps;
}

CGAL_END_NAMESPACE
