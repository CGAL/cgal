#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_sphere_traits_2.h>
#include <CGAL/Regular_triangulation_on_sphere_2.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/point_generators_2.h>
#include <fstream>
#include <CGAL/Timer.h>
#include <CGAL/number_utils.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Regular_triangulation_sphere_traits_2<K>             Gt;
typedef CGAL::Regular_triangulation_on_sphere_2<Gt>              RTOS;
typedef RTOS::Vertex_handle                             Vertex_handle;
typedef RTOS::Face_handle                                 Face_handle;
typedef RTOS::Point                                             Point;
typedef RTOS::Faces_iterator                            Face_iterator;
typedef RTOS::Vertices_iterator                           Vertex_iterator;
typedef RTOS::Locate_type                                 Locate_type;
typedef RTOS::Edge                                               Edge;
typedef RTOS::Hidden_vertices_iterator           Hidden_points_iterator;

template <class Stream>
Stream &write_vertices(Stream &out,std::vector<Vertex_handle> &t);

template <class Stream,  class Triangulation>
Stream &write_vertices(Stream &out, Triangulation &t);

template <class Stream>
Stream &write_points(Stream &out,std::vector<Point> &t);

template <class Stream>
Stream &write_face_to_off(Stream &out, Face_handle loc);

template <class Stream, class Triangulation>
Stream &write_triangulation_to_off_2(Stream &out,Stream &out2, Triangulation &t);

template <class Stream, class Triangulation>
Stream &write_hpts_to_off(Stream &out, Triangulation &t);

template <class Stream, class Triangulation>
Stream &write_triangulation_to_off(Stream &out, Triangulation &t);


int main()
{
	RTOS tr_2_h;
	RTOS tr_2_s;
	RTOS tr;
	RTOS tr2_1;
	RTOS tr2_2;
	RTOS tr2_3;
	RTOS tr1_1;
	RTOS tr0_1;
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
	
	std::list<Point> l; l.push_back(P15);
	l.push_back(P9); l.push_back(P10); l.push_back(P19);
	l.push_back(P0); l.push_back(P1); l.push_back(P2);
	
	std::vector<Point> v; v.push_back(P15);
	v.push_back(P9); v.push_back(P10); v.push_back(P19);
	v.push_back(P0); v.push_back(P1); v.push_back(P2);
	
	/*------------------------------*/
	/*         INSERTION            */
	/*------------------------------*/
	
	
	//dim-1
	std::cout<<"test insertion dim-1"<<std::endl;
	
	RTOS T0_0;
	assert(T0_0.dimension()==-2);
	v01=T0_0.insert(P1);
	T0_0.insert(P1);
	assert(T0_0.number_of_vertices()==1);assert(T0_0.number_of_hidden_vertices()==0);assert(T0_0.dimension()==-1);
	vv31=T0_0.insert(P31);//hidden
	assert(T0_0.number_of_vertices()==1);assert(T0_0.number_of_hidden_vertices()==1);assert(T0_0.dimension()==-1);
	vv34=T0_0.insert(P34);//hidding
	assert(T0_0.number_of_vertices()==1);assert(T0_0.number_of_hidden_vertices()==2);assert(T0_0.dimension()==-1);
	// tr.insert(P31);//already hidden            here ve store 2 versions of the hidden vertex
	//assert(tr.number_of_vertices()==1);//assert(tr0_1.number_of_hidden_vertices()==2);
	//assert(tr.dimension()==-1);
	//assert(tr0_1.is_valid());
	
	//dim0
	std::cout<<"test insertion dim0"<<std::endl;
	tr0_1.insert(P30);
	tr0_1.insert(P32);
	assert(tr0_1.number_of_vertices()==2);assert(tr0_1.number_of_hidden_vertices()==0);assert(tr0_1.dimension()==0);
	tr0_1.insert(P30);//vertex
	tr0_1.insert(P32);//vertex
	assert(tr0_1.number_of_vertices()==2);assert(tr0_1.number_of_hidden_vertices()==0);assert(tr0_1.dimension()==0);
	tr0_1.insert(P0);//outside affine hull hidding
	assert(tr0_1.number_of_vertices()==2);assert(tr0_1.number_of_hidden_vertices()==1);assert(tr0_1.dimension()==0);
	tr0_1.insert(P1);//outside affine hull hidding
	std::cout<<tr0_1.dimension()<<std::endl;
	assert(tr0_1.number_of_vertices()==2);assert(tr0_1.number_of_hidden_vertices()==2);assert(tr0_1.dimension()==0);
	tr0_1.insert(P31);//outside affine hull hidden
	assert(tr0_1.number_of_vertices()==2);assert(tr0_1.number_of_hidden_vertices()==3);assert(tr0_1.dimension()==0);
	tr0_1.insert(P33);//hidding
	assert(tr0_1.number_of_vertices()==2);assert(tr0_1.number_of_hidden_vertices()==4);assert(tr0_1.dimension()==0);
	tr0_1.insert(P32);//allready hidden
	assert(tr0_1.number_of_vertices()==2);//assert(tr0_1.number_of_hidden_vertices()==4);
	assert(tr0_1.dimension()==0);
	//assert(tr0_1.is_valid());
	
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
	v0=tr0_1.insert(P0);
	tr0_1.insert(P3);
	tr0_1.insert(P1);//dim 1 half circle
	assert(tr0_1.dimension()==1); assert(tr0_1.number_of_vertices()==3);
	assert(tr0_1.is_valid());
	tr0_1.insert(P0);
	tr0_1.insert(P3);
	tr0_1.insert(P1);
	tr0_1.insert(P2);//in edge
	tr0_1.insert(P30);//vertex hidden
	assert(tr0_1.number_of_vertices()==4);assert(tr0_1.number_of_hidden_vertices()==1);
	assert(tr0_1.is_valid());
	tr0_1.insert(P32);//edge hidden
	vv33=tr0_1.insert(P33);//outside hidding P0
	assert(tr0_1.number_of_hidden_vertices()==3); assert(tr0_1.number_of_vertices()==4);
	v37=tr0_1.insert(P37);//vertex hidding P2 and P1
	assert(tr0_1.number_of_hidden_vertices()==5); assert(tr0_1.number_of_vertices()==3);
	assert(tr0_1.is_valid());
	//full circle todo
	
	//opposed vertices case
	RTOS tr1_2;
	tr1_2.insert(P0);
	tr1_2.insert(P2);
	tr1_2.insert(P4);//half circle
	assert(tr1_2.dimension()==1);
	tr1_2.insert(P1);
	tr1_2.insert(P3);;
	assert(tr1_2.number_of_vertices()==5);
	tr1_2.insert(P5);
	assert(tr1_2.is_valid());
	
	//non-coplanar with _sphere
	tr1_1.insert(P15);
	tr1_1.insert(P9);
	tr1_1.insert(P10);
	assert(tr1_1.dimension()==1);
	tr1_1.insert(P40);//vertex hidden
	tr1_1.insert(P41);//vertex_hidding
	assert(tr1_1.dimension()==1);assert(tr1_1.number_of_vertices()==3);assert(tr1_1.number_of_hidden_vertices()==2);
	assert(tr1_1.is_valid());
	tr1_1.insert(P42);//outside affine hull hidden
	assert(tr1_1.dimension()==1);assert(tr1_1.number_of_hidden_vertices()==3);
	tr1_1.insert(P43);//outside affine hull hidding case 1lt42n
	// write_hpts_to_off(hpts,tr1_1);
	//write_vertices(verts,tr1_1);
	//write_triangulation_to_off_2(pos,neg,tr1_1);
	assert(tr1_1.dimension()==1);assert(tr1_1.number_of_hidden_vertices()==4);
	assert(tr1_1.is_valid());
	
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
	assert(tr2_2.number_of_hidden_vertices()==0);
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
	//must write a test for excessive hidding corner point
	std::cout<<"insertion dim 2"<<std::endl;
	//with full circle
	v1=tr2_1.insert(P1);//edge boundary
	v3=tr2_1.insert(P2);//
	v4=tr2_1.insert(P4);//
	v6=tr2_1.insert(P6);//
	assert(tr2_1.is_valid());
	v3=tr2_1.insert(P3);//vertex of boundary
	v31=tr2_1.insert(P31);//vertex of boundary hidden
	v34=tr2_1.insert(P34);//vertex of boundary hidding
	assert(tr2_1.is_valid());
	assert(tr2_1.number_of_hidden_vertices()==2);assert(tr2_1.number_of_vertices()==8);
	v33=tr2_1.insert(P33);//edge of boundary hidding
	v32=tr2_1.insert(P32);//edge of boundary hidden
	assert(tr2_1.number_of_hidden_vertices()==4);assert(tr2_1.number_of_vertices()==8);
	v19=tr2_1.insert(P19);//outside convex_hullx
	v42=tr2_1.insert(P42);
	assert(tr2_1.is_valid());
	
	//with half circle
	tr2_2.insert(P36); //edge 
	tr2_2.insert(P5); //outside
	tr2_2.insert(P19);//full sphere
	assert(tr2_2.is_valid());
	
	
	//general case
	std::cout<<"general case"<<std::endl;
	
	vv1=tr_2_s.insert(p1);
	vv2=tr_2_s.insert(p2);
	vv3=tr_2_s.insert(p3);
	vv4=tr_2_s.insert(p4);
	assert(tr_2_s.number_of_vertices()==4);
	assert(tr_2_s.is_valid());
	//write_triangulation_to_off_2(pos,neg,tr2_3);
	// write_vertices(verts,tr2_3);
	// write_hpts_to_off(hpts,tr2_3);
	//VERTEX    
	v20 =tr_2_s.insert(P20);//hidden 
	assert(tr_2_s.number_of_vertices()==4); assert(tr_2_s.number_of_hidden_vertices()==1);
	v21 = tr_2_s.insert(P21);//exactly on point
	assert(tr_2_s.number_of_vertices()==4); assert(tr_2_s.number_of_hidden_vertices()==1);
	v22 = tr_2_s.insert(P22);//hidding init vertex (1,1,1)
	assert(tr_2_s.number_of_vertices()==4); assert(tr_2_s.number_of_hidden_vertices()==2);
	//EDGE
	v23 = tr_2_s.insert(P23);  //hidden
	assert(tr_2_s.number_of_vertices()==4); assert(tr_2_s.number_of_hidden_vertices()==3);
	v24 = tr_2_s.insert(P24);  //on edge
	assert(tr_2_s.number_of_vertices()==5); assert(tr_2_s.number_of_hidden_vertices()==3);
	//FACE
	v25 = tr_2_s.insert(P25); //hidden
	assert(tr_2_s.number_of_vertices()==5); assert(tr_2_s.number_of_hidden_vertices()==4);
	v26 = tr_2_s.insert(P26); //on face, not hidding points
	assert(tr_2_s.number_of_vertices()==6); assert(tr_2_s.number_of_hidden_vertices()==4);
	v27 = tr_2_s.insert(P27); //on face, hidding point
	assert(tr_2_s.number_of_vertices()==6); assert(tr_2_s.number_of_hidden_vertices()==5);
	assert(tr_2_s.is_valid());
	
	/*
	 std::cout << "test copy constructor" << std::endl;
	 
	 // test copy_constructor with non-empty 1-triangulation
	 RTOS T1_5_1( tr0_1 );
	 assert( T1_5_1.dimension() == tr0_1.dimension() );
	 assert( T1_5_1.number_of_vertices() == tr0_1.number_of_vertices());
	 assert( T1_5_1.number_of_hidden_vertices()== 
	 tr0_1.number_of_hidden_vertices() );
	 assert( T1_5_1.is_valid() );
	 
	 // Test assignment operator
	 RTOS T1_5_2 = tr0_1;
	 assert( T1_5_1.dimension() == tr0_1.dimension() );
	 assert( T1_5_1.number_of_vertices() == tr0_1.number_of_vertices());
	 assert( T1_5_1.number_of_hidden_vertices()== 
	 tr0_1.number_of_hidden_vertices() );
	 assert( T1_5_1.is_valid() );
	 
	 // test copy_constructor with non-empty 2-triangulation
	 RTOS T2_3_1( tr2_2 );
	 assert( T2_3_1.dimension() == tr2_2.dimension());
	 assert( T2_3_1.number_of_vertices() == tr2_2.number_of_vertices());
	 assert( T2_3_1.number_of_hidden_vertices()== 
	 tr2_2.number_of_hidden_vertices() ); 
	 assert( T2_3_1.is_valid() );
	 
	 // test assignment operator
	 RTOS T2_3_4 =  tr2_2;
	 assert( T2_3_1.dimension() == tr2_2.dimension());
	 assert( T2_3_1.number_of_vertices() == tr2_2.number_of_vertices());
	 assert( T2_3_1.number_of_hidden_vertices()== 
	 tr2_2.number_of_hidden_vertices() ); 
	 assert( T2_3_1.is_valid() );
	 
	 // make sure push_back exists and does the same thing as insert
	 assert( tr_2_s.push_back(P27) == v27 );
	 assert( tr_2_s.number_of_vertices() == 6 );
	 
	 // test list iterator insert
	 RTOS T2_5;
	 assert( T2_5.insert(l.begin(), l.end()) == 7 );
	 assert( T2_5.dimension() == 2 );
	 assert( T2_5.number_of_vertices() == 7 );
	 assert( T2_5.is_valid() );
	 
	 // test vector iterator insert
	 RTOS T2_6;
	 assert( T2_6.insert(v.begin(), v.end()) == 7 );
	 assert( T2_6.dimension() == 2 );
	 assert( T2_6.number_of_vertices() == 7 );
	 assert( T2_6.is_valid() );
	 /*
	 
	 // test flip
     std::cout << "    test flip " << std::endl;
     Triangul T2_8;
     T2_8.insert(Point(0,0,1));
     T2_8.insert(Point(1,0,1));
     T2_8.insert(Point(1,1,1));
     T2_8.insert(Point(0,1,1));
     ff = T2_8.locate(Point(1,1,2),lt,li);
     assert(lt == Triangul::EDGE);
     assert(!T2_8.is_infinite(ff));
     Face_handle f2 = ff->neighbor(li);
     assert(!T2_8.is_infinite(f2));
     T2_8.flip(ff,0);
     assert( T2_8.is_valid() );
     /*
	 
	 /********************/
	/******** I/O *******/
	/*
	 std::cout << "test output to a file" << std::endl;
	 std::ofstream of0_0("T00.triangulation", std::ios::out);
	 CGAL::set_ascii_mode(of0_0); 
	 of0_0 << T0_0; of0_0.close();
	 std::ofstream of0_1("T01.triangulation");
	 CGAL::set_ascii_mode(of0_1); 
	 of0_1 << tr0_1; of0_1.close();
	 std::ofstream of1_2("T12.triangulation");
	 CGAL::set_ascii_mode(of1_2); 
	 of1_2 << tr1_2; of1_2.close();
	 std::ofstream of1_5("T11.triangulation");
	 CGAL::set_ascii_mode(of1_5); 
	 of1_5 << tr1_1; of1_5.close();
	 std::ofstream of1_6("T21.triangulation");
	 CGAL::set_ascii_mode(of1_6); 
	 of1_6 << tr2_1; of1_6.close();
	 std::ofstream of2_1("T22.triangulation");
	 CGAL::set_ascii_mode(of2_1); 
	 of2_1 << tr2_2; of2_1.close();
	 std::ofstream of2_3("T23.triangulation");
	 CGAL::set_ascii_mode(of2_3); 
	 of2_3 << tr2_3; of2_3.close();
	 std::ofstream of2_5("T25.triangulation");
	 CGAL::set_ascii_mode(of2_5); 
	 of2_5 << T2_5; of2_5.close();
	 std::ofstream of2_6("T26.triangulation");
	 CGAL::set_ascii_mode(of2_6); 
	 of2_6 << T2_6; of2_6.close();
	 /*
	 std::cout << "test input from a file" << std::endl;
	 std::ifstream if0_0("T00.triangulation"); CGAL::set_ascii_mode(if0_0);
	 RTOS T0_0_copy;   if0_0 >> T0_0_copy;
	 assert( T0_0_copy.is_valid() &&
	 T0_0_copy.number_of_vertices() == T0_0.number_of_vertices() );
	 std::ifstream if0_1("T01.triangulation"); CGAL::set_ascii_mode(if0_1);
	 RTOS T0_1_copy; if0_1 >> T0_1_copy;
	 assert( T0_1_copy.is_valid() &&
	 T0_1_copy.number_of_vertices() == tr0_1.number_of_vertices() );
	 std::ifstream if1_2("T12.triangulation"); CGAL::set_ascii_mode(if1_2); 
	 RTOS T1_2_copy; if1_2 >> T1_2_copy;
	 assert( T1_2_copy.is_valid() &&
	 T1_2_copy.number_of_vertices() == tr1_2.number_of_vertices() );
	 std::ifstream if1_5("T11.triangulation"); CGAL::set_ascii_mode(if1_5); 
	 RTOS T1_1_copy; if1_5 >> T1_2_copy;
	 assert( T1_1_copy.is_valid() &&
	 T1_1_copy.number_of_vertices() == tr1_1.number_of_vertices() );
	 std::ifstream if1_6("T21.triangulation"); CGAL::set_ascii_mode(if1_6);
	 RTOS T2_1_copy; if1_6 >> T2_1_copy;
	 assert( T2_1_copy.is_valid() &&
	 T2_1_copy.number_of_vertices() == tr2_1.number_of_vertices() );
	 std::ifstream if2_1("T22.triangulation"); CGAL::set_ascii_mode(if2_1);
	 RTOS T2_2_copy; if2_1 >> T2_2_copy;
	 assert( T2_2_copy.is_valid() &&
	 T2_2_copy.number_of_vertices() == tr2_2.number_of_vertices() );
	 std::ifstream if2_3("T23.triangulation"); CGAL::set_ascii_mode(if2_3);
	 RTOS T2_3_copy; if2_3 >> T2_3_copy;
	 assert( T2_3_copy.is_valid() &&
	 T2_3_copy.number_of_vertices() == tr2_3.number_of_vertices() );
	 std::ifstream if2_5("T25.triangulation"); CGAL::set_ascii_mode(if2_5); 
	 RTOS T2_5_copy; if2_5 >> T2_5_copy;
	 assert( T2_5_copy.is_valid() &&
	 T2_5_copy.number_of_vertices() == T2_5.number_of_vertices() );
	 std::ifstream if2_6("T26.triangulation"); CGAL::set_ascii_mode(if2_6);
	 RTOS T2_6_copy; if2_6 >> T2_6_copy;
	 assert( T2_6_copy.is_valid() &&
	 T2_6_copy.number_of_vertices() == T2_6.number_of_vertices() );
	 */
	/*
	 std::cout<<"test_remove"<<std::endl;
	 //remove dim2
	 tr_2_s.remove(v27); 
	 assert(tr_2_s.number_of_hidden_vertices()==4);assert(tr_2_s.number_of_vertices()==6);
	 tr_2_s.remove(v26); 
	 assert(tr_2_s.number_of_hidden_vertices()==4);assert(tr_2_s.number_of_vertices()==5);//a bug occurs sometimes here
	 tr_2_s.remove(v25); 
	 assert(tr_2_s.number_of_hidden_vertices()==3);assert(tr_2_s.number_of_vertices()==5);
	 tr_2_s.remove(v24);  
	 assert(tr_2_s.number_of_hidden_vertices()==3);assert(tr_2_s.number_of_vertices()==4);
	 tr_2_s.remove(v23);  
	 assert(tr_2_s.number_of_hidden_vertices()==2);assert(tr_2_s.number_of_vertices()==4);
	 assert(tr_2_s.is_valid());
	 
	 //full->half
	 tr2_1.remove(v11);//->half sphere but reinsertion of an hidden vertex restores it
	 tr2_1.remove(v42); //->half here 
	 assert(tr2_1.number_of_vertices()==8);assert(tr2_1.number_of_hidden_vertices()==4);
	 tr2_1.remove(v33);
	 
	 //dim2-dim1
	 tr2_1.remove(v19);
	 assert(tr2_1.dimension()==1);assert(tr2_1.is_valid());
	 
	 tr2_3.remove(v9);
	 assert(tr2_1.dimension()==1);assert(tr2_3.is_valid());
	 tr_2_s.remove(v22); //->dim1 but reinsertion of hidden vertex restores dim2
	 assert(tr_2_s.number_of_hidden_vertices()==1);assert(tr_2_s.number_of_vertices()==4);assert(tr_2_s.dimension()==2);
	 tr_2_s.remove(vv2); //->dim1
	 assert(tr_2_s.number_of_hidden_vertices()==1);assert(tr_2_s.number_of_vertices()==3);assert(tr_2_s.dimension()==1);
	 assert(tr_2_s.is_valid());
	 //  write_triangulation_to_off_2(pos,neg,tr2_1);
	 //  write_vertices(verts,tr2_1);
	 //  write_hpts_to_off(hpts,tr2_1);
	 
	 //dim1
	 tr2_1.remove(v31);//hidden
	 assert(tr2_1.number_of_hidden_vertices()==2);
	 tr2_1.remove(v34);//unhidding
	 assert(tr2_1.number_of_hidden_vertices()==1);
	 assert(tr2_1.is_valid());
	 tr2_1.remove(v3);
	 tr2_1.remove(v4);
	 tr2_1.remove(v5);
	 assert(tr2_1.is_valid());
	 
	 tr0_1.remove(v37);                                                 //assertion fault due to removing
	 //of faces_begin_2 that was cheating to acces faces in dim1
	 //right solition is to use edges_begin
	 //I implemented it but it causes this bug
	 //I don't know why
	 //maybee comparing with previous version with face will help
	 assert(tr0_1.is_valid());
	 //tr0_1.remove(vv33); // bug here because an hidden vertex (1,0,0) occurs twice in a list (also in 2D package)
	 //remove 1 of them first :
	 tr0_1.remove(v0);
	 tr0_1.remove(vv33);
	 assert(tr0_1.is_valid());
	 
	 //dim1->dim0
	 tr2_1.remove(v6);
	 tr2_1.remove(vv0);
	 assert(tr2_1.dimension()==1);
	 tr2_1.remove(v32);
	 assert(tr2_1.dimension()==0);
	 
	 //dim0
	 Point P50=Point(0.25,0.5,0);
	 Vertex_handle v50=tr2_1.insert(P31);
	 tr2_1.remove(v1);
	 assert(tr2_1.dimension()==0);
	 tr2_1.remove(v50);
	 assert(tr2_1.dimension()==-1);
	 
	 //dim-1
	 T0_0.remove(vv34);
	 assert(T0_0.dimension()==-2);assert(T0_0.number_of_vertices()==0);assert(T0_0.number_of_hidden_vertices()==0);
	 //test if hidden vertices reinserted in dim -1, currently they are deleted
	 //assert(tr.number_of_hidden_vertices()==1);
	 //tr.remove(vv31);
	 //assert(tr.number_of_hidden_vertices()==0);assert(tr.dimension()==-1);
	 //tr.remove(v01);
	 //assert(tr.dimension()==-2);
	 
	 
	 std::cout<<"revomal test passed successfully"<<std::endl;
	 
	 /**********************/
	/******* MOVE *********/
	/*  
	 std::cout << "    moves" << std::endl;
	 
	 std::cout << "    degenerate cases: " << std::endl;
	 
	 Triangul TM_0, TM_1;
	 Vertex_handle tmv1 = TM_0.insert(Point(0,0));
	 Vertex_handle tmv2 = TM_0.insert(Point(1,0));
	 Vertex_handle tmv3 = TM_0.insert(Point(2,0));
	 Vertex_handle tmv4 = TM_0.insert(Point(1,1));
	 assert(TM_0.dimension() == 2);
	 TM_0.move(tmv4, Point(2, 1));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 2);
	 TM_0.move(tmv4, Point(3, 0));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 1);
	 TM_0.move(tmv3, Point(1, 1));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 2);
	 TM_0.move(tmv3, Point(-1, 0));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 1);
	 TM_0.move(tmv2, Point(-1, 0, 2));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 1);
	 TM_0.move(tmv2, Point(-1, 0, 2));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 1);
	 TM_0.move(tmv2, Point(-1, 0, 4));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 1);
	 TM_0.move(tmv2, Point(-1, 0, 2));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 1);
	 TM_0.move(tmv2, Point(-1, 1, 2));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 2);
	 TM_0.move(tmv1, Point(0, 2));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 2);
	 TM_0.move(tmv1, Point(0, 1));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 2);
	 TM_0.move(tmv1, Point(0, 0));
	 assert(TM_0.tds().is_valid());
	 assert(TM_0.is_valid());
	 assert(TM_0.dimension() == 2);
	 assert(!TM_0.move(tmv1, Point(3, 0)));
	 
	 std::cout << "    non-degenerate cases: " << std::endl;
	 // non-degenerate cases
	 CGAL::Random_points_in_square_2<Point> g(50);
	 std::list<Point> points;
	 for(int count=0; count<50; count++) points.push_back(*g++);
	 TM_1.insert(points.begin(), points.end());
	 for(int i=0; i<50; i++) {
	 for(Finite_vertices_iterator fvi = TM_1.finite_vertices_begin();
	 fvi != TM_1.finite_vertices_end(); fvi++) {
	 Point p = *g++;
	 TM_1.move(fvi, p);
	 assert(TM_1.tds().is_valid());
	 assert(TM_1.is_valid());
	 }
	 */
	
	
	
	
	
	//Test random insertion in or on the sphere on a restricted area of the sphere
	//good way to find bugs because points are close
	//at the moment, in sphere we have a bug,
	//this bug must be related with a case documented in my notes
	//this case is not tested, we talk about it with manuel
	//the source is in the first face of the find conflict
	//it has been fixed with get conflict but it may be the source of another bug
	//when the center is outside the sphere we saw vertives inside the lower envelope faces not hidden
	//bug occurs at nb = 107
	tr_2_s.clear();
	assert(tr_2_s.is_valid());
	// tr_2_s.insert_four_init_vertices();
	assert(tr_2_s.is_valid());
	
	CGAL::Random random(7);
	typedef CGAL::Creator_uniform_3<double, Point> Creator;
	int nb_of_pts=107;
	std::cout<<"\tINSERTION OF "<<nb_of_pts<<" RANDOM POINTS\n"<<std::endl;
	
	std::ofstream fout("test.off");
	
	//Without hidden points
	CGAL::Random_points_in_sphere_3<Point, Creator> on_sphere(1, random);
	for (int count=0; count<nb_of_pts; count++) {
		
		Point p = *on_sphere; on_sphere ++;
		p = Point(p.x(),p.y(),p.z());
		
		if(p.z()>0.1 && p.y()>0.1 && p.x()>0.1){
			std::cout<<count<<std::endl;
			tr_2_s.insert(p);
		}
		
		assert(tr_2_s.is_valid());
		
		std::ofstream pos2("pos2.off");
		std::ofstream neg2("neg2.off");
		write_vertices(verts,tr_2_s);
		write_hpts_to_off(hpts,tr_2_s);
		write_triangulation_to_off_2(pos2,neg2,tr_2_s);
		
	}
	
    write_vertices(verts,tr_2_s);
    write_hpts_to_off(hpts,tr_2_s);
    write_triangulation_to_off_2(pos,neg,tr_2_s);
	
	
    assert(tr_2_s.number_of_hidden_vertices()==0);
	assert(tr_2_s.is_valid());
	
	
	return 0;
}


template <class Stream, class Triangulation>
Stream &write_triangulation_to_off(Stream &out, Triangulation &t) {
	std::vector<Face_handle> faces;
	
	// Points of triangulation
	for (Face_iterator it = t.faces_begin(); it != t.faces_end(); it++) {
		for (int i=0 ; i<3 ; i++) {
			Point p = it->vertex(i)->point();
			out << p.x() << " " 
			<< p.y() << " " 
			<< p.z() << std::endl;
		}
    }
	
	return out;
}


template <class Stream, class Triangulation>
Stream &write_hpts_to_off(Stream &out, Triangulation &t) {
	
	
	// Points of triangulation
	for (Hidden_points_iterator it = t.hidden_vertices_begin(); it != t.hidden_vertices_end(); it++) {
		out << (it)->point().x() << " " 
		<< (it)->point().y() << " " 
		<< (it)->point().z() << std::endl;
		
	}
	
	return out;
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

template <class Stream>
Stream &write_face_to_off(Stream &out, Face_handle loc){
	
	for (int i=0 ; i<3 ; i++) {
		Point p = loc->vertex(i)->point();
		out << p.x() << " " 
	    << p.y() << " " 
	    << p.z() << std::endl;
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

template <class Stream>
Stream &write_vertices(Stream &out,std::vector<Vertex_handle> &t)
{
	for(std::vector<Vertex_handle>::iterator it= t.begin(); it!= t.end(); ++it)
	{
		Point p=(*it)->point();
		out << p.x() << " " 
		<< p.y() << " " 
		<< p.z() << std::endl;
		
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