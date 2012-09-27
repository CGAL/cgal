#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_on_sphere_2.h>
#include <CGAL/Triangulation_sphere_traits_2.h>
#include <CGAL/enum.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/number_utils.h>

#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Triangulation_sphere_traits_2<K>          Gt;
typedef CGAL::Triangulation_on_sphere_2<Gt>               TOS;


typedef TOS::Point             Point;
typedef TOS::Face_handle      Face_handle;
typedef TOS::Vertex_handle Vertex_handle;
typedef TOS::Faces_iterator Face_iterator;
typedef TOS::Locate_type  Locate_type;
typedef TOS::Edge       Edge;
typedef TOS::Vertices_iterator                           Vertex_iterator;

template <class Stream>
Stream &write_points(Stream &out,std::vector<Point> &t);

template <class Stream,  class Triangulation>
Stream &write_vertices(Stream &out, Triangulation &t);

template <class Stream>
Stream &write_faces_to_off(Stream &out, std::vector<Face_handle> t);

template <class Stream, class Triangulation>
Stream &write_triangulation_to_off(Stream &out, Triangulation &t);

template <class Stream, class Triangulation>
Stream &write_triangulation_to_off_2(Stream &out,Stream &out2, Triangulation &t);


int main(){
  TOS tr_2_h;
  TOS tr_2_s;
  TOS tr;
  TOS tr2_1;
  TOS tr2_2;
  TOS tr2_3;
  TOS tr1_1;
  TOS tr0_1;
  Locate_type lt;
  int li;
  Face_handle Loc_result;
  Face_handle start;

  std::ofstream verts("verts.off");
  std::ofstream pos("pos.off");
  std::ofstream neg("neg.off");
  //std::ofstream pos2("pos2.off");
  //std::ofstream neg2("neg2.off");
  std::ofstream hpts("hpts.off");
  //std::ofstream hpts2("hpts2.off");
 

  /*--------------------------*/
  //       TEST INPUTS        //
  /*--------------------------*/
  //             SPHERES
  Point S1=Point(0,0,0);

 
  //             POINTS
  Point p1(1,0,0);
  Point p2(0,1,0);
  Point p3(0,0,1);
  Point p4(-1/sqrt(3.),-1/sqrt(3.),-1/sqrt(3.));  

  //on a circle z=0
  Point P0=Point(1,0,0);
  Point P1=Point(1/sqrt(2.),1/sqrt(2.),0);
  Point P2=Point(0,1,0);
  Point P3=Point(-1/sqrt(2.),1/sqrt(2.),0);
  Point P4=Point(-1,0,0);
  Point P5=Point(-1/sqrt(2.),-1/sqrt(2.),0);
  Point P6=Point(1/sqrt(2.),-1/sqrt(2.),0);
  Point P30=Point(0.5,0,0);
  Point P31=Point(0.25,0.25,0);
  Point P32=Point(0.5,0.25,0);
  Point P33=Point(2,-0.25,0);
  //on a line
  Point P34=Point(1,1,0);
  Point P35=Point(1,-1,0);
  Point P36=Point(1,0.5,0);
  Point P37=Point(2,2,0);

  //top tetrahedron
  Point P15=Point(1/sqrt(2.),0,1/sqrt(2.));
  Point P9=Point(-0.5,-0.5,1/sqrt(2.));
  Point P10=Point(0,1/sqrt(2.),1/sqrt(2.));
  Point P11=Point(0,0,1);
  Point P12=Point(-0.5,0.5,1/sqrt(2.));
  Point P13=Point(0.5,-0.5,1/sqrt(2.));
  Point P14=Point(0,-0.5,1/sqrt(2.));//edge P9,P13
  //Point P15=Point(0.25/sqrt(2.),0.25/sqrt(2.),0.5/sqrt(2.));//under edge
  Point P16=Point(-0.5,0.5*(1+sqrt(2.)),2/sqrt(2.));//over edge
  Point P17=Point(0.5,0,0.5);
  Point P18=Point(-1/sqrt(2.),0,-1/sqrt(2.));
  Point P19=Point(0,0,-1);
  Point P40=Point(-0.25,-0.25,0.5/sqrt(2.));
  Point P41=Point(-1,-1,2/sqrt(2.));
  Point P42=Point(0,0,0.5);
  Point P43=Point(2,-0.5,1);

  Point P20=Point(0.5,0,0);  
  Point P21=Point(1,0,0);   
  Point P22=Point(2,0,0);
  Vertex_handle v20;
  Vertex_handle v21;
  Vertex_handle v22;

  Point P23=Point(0,0.25,0.25);  
  Point P24=Point(0,0.5,0.5);     
  Vertex_handle v23;
  Vertex_handle v24;

  Point P25=Point(0.75,0.25,0.1);  
  Point P26=Point(0.5,0.5,0.5);    
  Point P27=Point(1,1,1);          
  Vertex_handle v25;
  Vertex_handle v26;
  Vertex_handle v27;
  Vertex_handle v0;
  Vertex_handle vv1;
  Vertex_handle vv0;
  Vertex_handle v1;
  Vertex_handle v2;  
  Vertex_handle v3; 
  Vertex_handle v4;
  Vertex_handle v01;
  Vertex_handle vv2;  
  Vertex_handle vv3; 
  Vertex_handle vv4;
  Vertex_handle v5;
  Vertex_handle v6;
  Vertex_handle v7;
  Vertex_handle v8;
  Vertex_handle v9;
  Vertex_handle v10;
  Vertex_handle v11;
  Vertex_handle v12;
  Vertex_handle v13;
  Vertex_handle v14;
  Vertex_handle v15;
  Vertex_handle v33;
  Vertex_handle v16;
  Vertex_handle v31;
  Vertex_handle v34;
  Vertex_handle v32;
  Vertex_handle v19;
  Vertex_handle v42;
  Vertex_handle v37;
  Vertex_handle vv33;
  Vertex_handle vv34;
 Vertex_handle vv31;



  /*--------------------------*/
  //       TEST INPUTS        //
  /*--------------------------*/

  //             SPHERES

    //dim-1
  std::cout<<"test insertion dim-1"<<std::endl;
  assert(tr.dimension()==-2);
  v01=tr.insert(P1);
  tr.insert(P1);
  assert(tr.number_of_vertices()==1);assert(tr.dimension()==-1);
   
  //dim0
  std::cout<<"test insertion dim0"<<std::endl;
  tr0_1.insert(P30);
  tr0_1.insert(P32);
  assert(tr0_1.number_of_vertices()==2);assert(tr0_1.dimension()==0);
  tr0_1.insert(P30);//vertex
  tr0_1.insert(P32);//vertex
  assert(tr0_1.number_of_vertices()==2);assert(tr0_1.dimension()==0);

  //opposed vertices case
  tr0_1.clear();
  tr0_1.insert(P0);
  tr0_1.insert(P4);
  assert(tr0_1.dimension()==0);
  tr0_1.insert(P0);
  tr0_1.insert(P4);
  tr0_1.insert(P30);
  assert(tr0_1.dimension()==0);
  tr0_1.insert(P2);
  assert(tr0_1.dimension()==1);
  
  //dim1
  std::cout<<"test insertion dim1"<<std::endl;
  tr0_1.clear();

  tr1_1.insert(P42);//outside affine hull 
  assert(tr1_1.dimension()==2);assert(tr1_1.number_of_vertices()==4);
  
  //dim1->dim2
  std::cout<<"dim1 to dim2"<<std::endl;
  //with full circle
  tr2_1.clear();
  vv0=tr2_1.insert(P0);
  tr2_1.insert(P3);
  v5=tr2_1.insert(P5);
  v11=tr2_1.insert(P11);
  assert(tr2_1.dimension()==2);
  assert(tr2_1.is_valid());
  
  //with half circle
  tr2_2.insert(P0);
  tr2_2.insert(P6);
  tr2_2.insert(P1);
  assert(tr2_2.is_valid());
  tr2_2.insert(P2);
  assert(tr2_2.is_valid());
  tr2_2.insert(P11);
  assert(tr2_2.dimension()==2);
  assert(tr2_2.is_valid());
  
  //from non coplanar points with _sphere circle
  tr2_3.insert(P15);
  v9=tr2_3.insert(P9);
  tr2_3.insert(P10);
  tr2_3.insert(P11);
  assert(tr2_3.dimension()==2);
  assert(tr2_3.is_valid());

  //dim2
  std::cout<<"insertion dim 2"<<std::endl;
  //with full circle
  v1=tr2_1.insert(P1);//edge boundary
  write_triangulation_to_off_2(pos,neg,tr2_1);
  write_vertices(verts,tr2_1);
  v3=tr2_1.insert(P2);//
  v4=tr2_1.insert(P4);//
  v6=tr2_1.insert(P6);//
  v3=tr2_1.insert(P3);//vertex of boundary
  assert(tr2_1.is_valid());
  assert(tr2_1.number_of_vertices()==8);
  v19=tr2_1.insert(P19);//outside convex_hull
  v42=tr2_1.insert(P42);
  assert(tr2_1.is_valid());


  //with half circle
  tr2_2.insert(P36); //edge 
  tr2_2.insert(P5); //outside
  tr2_2.insert(P19);//full sphere
  assert(tr2_2.is_valid());
  
  //general case
  vv1=tr_2_s.insert(p1);
  vv2=tr_2_s.insert(p2);
  vv3=tr_2_s.insert(p3);
  vv4=tr_2_s.insert(p4);
  assert(tr_2_s.is_valid());
  //VERTEX    
  v20 =tr_2_s.insert(P20); 
  assert(tr_2_s.number_of_vertices()==4); 
  //EDGE
  v24 = tr_2_s.insert(P24);  //on edge
  assert(tr_2_s.number_of_vertices()==5);
  //FACE
  v25 = tr_2_s.insert(P25); 
  assert(tr_2_s.number_of_vertices()==6);
  v26 = tr_2_s.insert(P26); //on face
  assert(tr_2_s.number_of_vertices()==7);
  assert(tr_2_s.is_valid());

  /*--------------------------*/
  //          LOCATION        //
  /*--------------------------*/
  
  std::cout<<"\tLOCATION TEST :\n"<<std::endl;
  /*
  // DIM 0
  // DIM 1
  // DIM 2
  //March locate half sphere

  std::vector<Point> points;
 

  tr_2_h.clear();
  tr_2_h.insert_four_init_vertices_half_sphere();

  std::cout<<"Location half sphere1"<<std::endl;
 
  start=tr_2_h.faces_begin()->neighbor(2);
  tr_2_h.show_face(start);
  

  //start from a negative face
  Loc_result=tr_2_h.locate(P12,lt,li,start); //outside
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  Loc_result=tr_2_h.locate(P31,lt,li,start); //on positive face
  assert(lt==TOS::FACE);
  assert(tr_2_h.orientation(Loc_result)==1);
  Loc_result=tr_2_h.locate(P32,lt,li,start); //outside
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  assert(tr_2_h.orientation(Loc_result)==-1);
  Loc_result=tr_2_h.locate(P34,lt,li,start); //edge boundary
  std::vector<Face_handle> loc; loc.push_back(Loc_result);
  std::ofstream locc("loc.off");
  write_faces_to_off(locc,loc);
  std::cout<<P34<<std::endl;
  std::ofstream pos("pos.off");
  std::ofstream neg("neg.off");
  write_triangulation_to_off_2(pos,neg,tr_2_h); 
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  Loc_result=tr_2_h.locate(P12,lt,li,start); //outside opposite
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  assert(tr_2_h.orientation(Loc_result)==-1);
  Loc_result=tr_2_h.locate(P8,lt,li,start); //edge boundary
  assert(lt==TOS::EDGE);
  Loc_result=tr_2_h.locate(P9,lt,li,start); //outside
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  Loc_result=tr_2_h.locate(P10,lt,li,start); //point boundary
  assert(lt==TOS::VERTEX); 
  Loc_result=tr_2_h.locate(P37,lt,li,start); //on positive face
  assert(lt==TOS::FACE);
  assert(tr_2_h.orientation(Loc_result)==1);
  Loc_result=tr_2_h.locate(P11,lt,li,start); //vertex positive face
  assert(lt==TOS::VERTEX);
  Loc_result=tr_2_h.locate(P38,lt,li,start); //on positive face
  assert(lt==TOS::FACE);
  assert(tr_2_h.orientation(Loc_result)==1);
  Loc_result=tr_2_h.locate(P40,lt,li,start);//edge positive/positive
  assert(lt==TOS::EDGE);
  Loc_result=tr_2_h.locate(P36,lt,li,start); //outside coplanar with an edge of interior //most tricky case i saw
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  
  std::cout<<"Location half sphere2"<<std::endl;
  start=tr_2_h.faces_begin();
  tr_2_h.show_face(start);

  //start from a positive face 
  Loc_result=tr_2_h.locate(P12,lt,li,start); //outside
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  Loc_result=tr_2_h.locate(P31,lt,li,start); //on positive face
  assert(lt==TOS::FACE);
  assert(tr_2_h.orientation(Loc_result)==1);
  Loc_result=tr_2_h.locate(P32,lt,li,start); //outside
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  assert(tr_2_h.orientation(Loc_result)==-1);
  Loc_result=tr_2_h.locate(P12,lt,li,start); //outside opposite
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  assert(tr_2_h.orientation(Loc_result)==-1);
  Loc_result=tr_2_h.locate(P34,lt,li,start); 
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  Loc_result=tr_2_h.locate(P36,lt,li,start); //outside coplanar with an edge of interior
  std::cout<<lt<<std::endl;
  std::cout<<P36<<std::endl;
 
  assert(lt==TOS::OUTSIDE_CONVEX_HULL);
  Loc_result=tr_2_h.locate(P37,lt,li,start); //on positive face
  assert(lt==TOS::FACE);
  assert(tr_2_h.orientation(Loc_result)==1);
  Loc_result=tr_2_h.locate(P38,lt,li,start); //on positive face
  assert(lt==TOS::FACE);
  assert(tr_2_h.orientation(Loc_result)==1);
  Loc_result=tr_2_h.locate(P11,lt,li,start); //vertex positive face
  assert(lt==TOS::VERTEX);
  Loc_result=tr_2_h.locate(P40,lt,li,start);//edge positive/positive
  assert(lt==TOS::EDGE);
  
  // March_locate full sphere
  std::cout<<"Location full sphere :\n"<<std::endl;
 
  std::vector<Point> pts_to_locate;
  std::vector<Vertex_handle> init;

  tr_2_s.clear();
  init=tr_2_s.insert_four_init_vertices();
  assert(tr_2_s.is_valid());

  //general cases
  start=tr_2_s.faces_begin();
  Loc_result=tr_2_s.locate(P0,lt,li,start); //on point
  assert(lt==TOS::VERTEX);
  assert(tr_2_s.xy_equal(Loc_result->vertex(li)->point(),P0));
  Loc_result=tr_2_s.locate(P200,lt,li,start);                //coradial
  assert(lt==TOS::VERTEX);
  assert(tr_2_s.xy_equal(Loc_result->vertex(li)->point(),P200));
  tr_2_s.locate(P300,lt,li,start); //on edge
  assert(lt==TOS::EDGE);
  tr_2_s.locate(P500,lt,li,start);               //sector of edge
  assert(lt==TOS::EDGE);
  tr_2_s.locate(P1100,lt,li,start);           //sector of face
  assert(lt==TOS::FACE);

  //cases P is at the opposite of an element of the starting face
  start=init[2]->face()->neighbor(2) ;
  tr_2_s.show_face(start);

  tr_2_s.locate(P4,lt,li,start); //face opposite of a vertex
  assert(lt==TOS::FACE);
  tr_2_s.locate(P700,lt,li,start); //face opposite of a vertex
  assert(lt==TOS::FACE);
  tr_2_s.locate(P12,lt,li,start);//face opposite of a vertex
  assert(lt==TOS::FACE);
  tr_2_s.locate(P400,lt,li,start); //Edge opposite of an edge
  assert(lt==TOS::EDGE);
  tr_2_s.locate(P14,lt,li,start);//Edge opposite of an edge
  assert(lt==TOS::EDGE);
  tr_2_s.locate(P15,lt,li,start); //edge opposite of an edge
  assert(lt==TOS::EDGE);
  Loc_result=tr_2_s.locate(P13,lt,li,start); //vertex opposite of the starting face
  assert(lt==TOS::VERTEX);
  assert(tr_2_s.xy_equal(Loc_result->vertex(li)->point(),P13));


  /*--------------------------*/
  //           REMOVAL        //
  /*--------------------------*/
  // DIM 2 
  // DIM 0
  // DIM 1

}


template <class Stream, class Triangulation>
Stream &write_triangulation_to_off_2(Stream &out,Stream &out2, Triangulation &t){

  // Points of triangulation
  for (Face_iterator it = t.faces_begin(); it != t.faces_end(); it++) {  
    if(!it->is_negative()
       /*(t.orientation(it->vertex(0)->point(),it->vertex(1)->point(),it->vertex(2)->point())==1)*/
	 ){
	for (int i=0 ; i<3 ; i++) {
	  if(it->vertex(i)!=Vertex_handle())
	  {
	    Point p = it->vertex(i)->point();
	    out << p.x() << " " 
		<< p.y() << " " 
		<< p.z() << std::endl;
	  }
	}
      }
  else{
	for (int i=0 ; i<3 ; i++){
	  if(it->vertex(i)!=Vertex_handle())
	    {
	      Point p = it->vertex(i)->point();
	      out2 << p.x() << " " 
	           << p.y() << " " 
	           << p.z() << std::endl;
	    }	  
	  }
      }
    }
 

  return out;
}

template <class Stream, class Triangulation>
Stream &write_triangulation_to_off(Stream &out, Triangulation &t){
  std::cout<<"drawing begins"<<std::cout;
  // Points of triangulation
  for (Face_iterator it = t.faces_begin_2(); it != t.faces_end(); it++) {
    for (int i=0 ; i<3 ; i++) {
	if(it->vertex(i)!=Vertex_handle())
	  {
	Point p = it->vertex(i)->point();
	out << p.x() << " " 
	    << p.y() << " " 
	    << p.z() << std::endl;
	  }
      }
  }

  return out;
}




template <class Stream>
Stream &write_faces_to_off(Stream &out, std::vector<Face_handle> t){

  for(std::vector<Face_handle>::iterator it= t.begin(); it!= t.end(); ++it)
  {
    for (int i=0 ; i<3 ; i++) {
	Point p = (*it)->vertex(i)->point();
	out << p.x() << " " 
	    << p.y() << " " 
	    << p.z() << std::endl;
      }
  }

  return out;
}


template <class Stream>
Stream &write_points(Stream &out,std::vector<Point> &t)
{
 for(std::vector<Point>::iterator it= t.begin(); it!= t.end(); ++it)
  {
    out << (*it).x() << " " 
	<< (*it).y() << " " 
	<< (*it).z() << std::endl;
    
  }

}

template <class Stream,  class Triangulation>
Stream &write_vertices(Stream &out, Triangulation &t)
{
  
  for(Vertex_iterator it= t.vertices_begin(); it!= t.vertices_end(); ++it)
  {
    Point p=(it)->point();
    out << p.x() << " " 
	<< p.y() << " " 
	<< p.z() << std::endl;
    
  }
}
