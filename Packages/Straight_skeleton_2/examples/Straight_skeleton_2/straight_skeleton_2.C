
#include<CGAL/Simple_cartesian.h>
#include<CGAL/Straight_skeleton_builder_2.h>

typedef CGAL::Simple_cartesian<double>                   Kernel;
typedef CGAL::Straight_skeleton_2<Kernel>                Ssds;
typedef CGAL::Straight_skeleton_builder_traits_2<Kernel> BuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<BuilderTraits,Ssds> Builder;
typedef Ssds::Halfedge_iterator          Halfedge_iterator;
typedef Ssds::Vertex_handle              Vertex_handle;

typedef Kernel::Point_2                   Point_2;


int main()
{
  /*
  Point_2 outer[] = {  Point_2(0,0)
                     , Point_2(50,30)
		     , Point_2(100,0)
		     , Point_2(120,50)
		     , Point_2(90,100)
		     , Point_2(0,100)
		    } ;
  */	  
  Point_2 outer[] = {  Point_2(0,0)
                     , Point_2(20,0)
		     , Point_2(20,10)
		     , Point_2(0,10)
		    } ;  
 
  Point_2 hole[] = { Point_2(20,70)
                    ,Point_2(50,90)
		    ,Point_2(70,50)
		   } ;
		   
  // Straight Skeleton builder
  Ssds ss;

  Ssds ss2 = Builder().insert_CCB(outer,outer+4)
    //                       .insert_CCB(hole,hole+3)
   .proceed() ;
  ss = ss2;
  for ( Halfedge_iterator i = ss.halfedges_begin(); i != ss.halfedges_end(); ++i ){
    std::string edge_type = (i->is_bisector())? " bisector" : " border";
    std::cout << i->segment() << edge_type << std::endl;
  }
  return 0;
} 
