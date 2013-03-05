#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>

#include <CGAL/Triangulation_2_filtered_projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <iostream>
#include <algorithm>
#include <vector>

/* If you want to use a viewer, you can use one of the following file
 * depending if you use vtk or qglviewer. */
#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#else 
#ifdef CGAL_LCC_USE_VTK
#include "linear_cell_complex_3_viewer_vtk.h"
#endif
#endif

typedef CGAL::Linear_cell_complex<3> LCC_3;
typedef LCC_3::Dart_handle           Dart_handle;
typedef LCC_3::Point                 Point;
typedef LCC_3::FT                    FT;
typedef LCC_3::Traits Traits;

typedef CGAL::Triangulation_2_filtered_projection_traits_3<Traits> P_traits;
typedef CGAL::Triangulation_vertex_base_with_info_2<Dart_handle,P_traits> Vb;

struct Face_info {
  bool exist_edge[3];
  bool is_external;
};

typedef CGAL::Triangulation_face_base_with_info_2<Face_info,P_traits> Fb1;

typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fb1>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                   TDS;
typedef CGAL::No_intersection_tag                                     Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS,
                                                   Itag>              CDT;

bool is_external(CDT::Face_handle fh)
{
  return fh->info().is_external;
}

int number_of_existing_edge(CDT::Face_handle fh)
{
  unsigned res=0;
  for (int i=0; i<3; ++i)
    if (fh->info().exist_edge[i]) ++res;
  return res;
}

int get_free_edge(CDT::Face_handle fh)
{
  CGAL_assertion( number_of_existing_edge(fh)==2 );
  for (int i=0; i<3; ++i)
    if (!fh->info().exist_edge[i]) return i;
  
  CGAL_assertion(false);
  return -1;
}

void constrained_delaunay_triangulation(LCC_3 &lcc, Dart_handle d1)
{
  CGAL::set_ascii_mode(std::cout);
  std::cout<<"Vertices: ";
  for (LCC_3::Vertex_attribute_const_range::iterator
         v=lcc.vertex_attributes().begin(),
         vend=lcc.vertex_attributes().end(); v!=vend; ++v)
    std::cout << v->point() << "; ";
  std::cout<<std::endl;
 
  LCC_3::Vector normal = CGAL::compute_normal_of_cell_2(lcc,d1);
  P_traits cdt_traits(normal);
  CDT cdt(cdt_traits); 
    
  //inserting the constraints edge by edge
  LCC_3::Dart_of_orbit_range<1>::iterator
    it(lcc.darts_of_orbit<1>(d1).begin());

   CDT::Vertex_handle previous=NULL, first=NULL, vh=NULL;

   for (LCC_3::Dart_of_orbit_range<1>::iterator
          itend(lcc.darts_of_orbit<1>(d1).end()); it!=itend; ++it)
   {     
     vh = cdt.insert(LCC_3::point(it));
     vh->info()=it;
     if( first==NULL ){
       first=vh;
     }
     if( previous!=NULL){
       CGAL_assertion( previous !=vh );
       cdt.insert_constraint(previous,vh);
     }

     previous=vh;
   }
   cdt.insert_constraint(previous,first);
   CGAL_assertion(cdt.is_valid());
   
   // sets mark is_external
   for( CDT::All_faces_iterator fit = cdt.all_faces_begin(),
          fitend = cdt.all_faces_end(); fit != fitend; ++fit)
   {
     fit->info().is_external = false;
     fit->info().exist_edge[0]=false;
     fit->info().exist_edge[1]=false;
     fit->info().exist_edge[2]=false;
   }
	
   std::queue<CDT::Face_handle> face_queue;
      
   face_queue.push(cdt.infinite_vertex()->face());
   while(! face_queue.empty() )
   {
     CDT::Face_handle fh = face_queue.front();
     face_queue.pop();
     if(!fh->info().is_external)
     {
       fh->info().is_external = true;
       for(int i = 0; i <3; ++i)
       {
         if(!cdt.is_constrained(std::make_pair(fh, i)))
         {
           face_queue.push(fh->neighbor(i));
         }
       }
     }
   }

   for( CDT::Finite_edges_iterator eit = cdt.finite_edges_begin(),
          eitend = cdt.finite_edges_end(); eit != eitend; ++eit)
   {
     CDT::Face_handle fh = eit->first;
     int index = eit->second;
     CDT::Face_handle opposite_fh = fh->neighbor(index);
     if(cdt.is_constrained(std::make_pair(fh, index)))
     {
       fh->info().exist_edge[index]=true;
       opposite_fh->info().exist_edge[cdt.mirror_index(fh,index)]=true;
       
       if ( !fh->info().is_external && number_of_existing_edge(fh)==2 )
         face_queue.push(fh);
       if ( !opposite_fh->info().is_external &&
            number_of_existing_edge(opposite_fh)==2 )
         face_queue.push(opposite_fh);
     }
   }
   
   while( !face_queue.empty() )
   {
     CDT::Face_handle fh = face_queue.front();
     face_queue.pop();
     CGAL_assertion( number_of_existing_edge(fh)>=2 ); // i.e. ==2 or ==3
     CGAL_assertion( !fh->info().is_external );
     
     if (number_of_existing_edge(fh)==2)
     {
       int index = get_free_edge(fh);
       CDT::Face_handle opposite_fh = fh->neighbor(index);

       CGAL_assertion( !fh->info().exist_edge[index] );
       CGAL_assertion( !opposite_fh->info().
                       exist_edge[cdt.mirror_index(fh,index)] );
       
       const CDT::Vertex_handle va = fh->vertex(cdt. cw(index));
       const CDT::Vertex_handle vb = fh->vertex(cdt.ccw(index));
       
       Dart_handle ndart=
         CGAL::insert_cell_1_in_cell_2(lcc,va->info(),vb->info());         
       va->info()=ndart->beta(2);

       fh->info().exist_edge[index]=true;
       opposite_fh->info().exist_edge[cdt.mirror_index(fh,index)]=true;
       
       if ( !opposite_fh->info().is_external &&
            number_of_existing_edge(opposite_fh)==2 )
         face_queue.push(opposite_fh);
     }
   }   
}

Dart_handle make_facet(LCC_3& lcc,const std::vector<Point>& points)
{
  Dart_handle d =
    CGAL::make_combinatorial_polygon<LCC_3>(lcc,(unsigned int)points.size());
  for (unsigned int i=0; i<points.size(); ++i)
  {
    lcc.set_vertex_attribute_of_dart(d, lcc.create_vertex_attribute(points[i]));
    d=d->beta(1);
  }
  return d;
}


int main()
{
  LCC_3 lcc;

  // Create one tetrahedra.
  Dart_handle d1 = lcc.make_tetrahedron(Point(-1, 0, 0), Point(0, 2, 0),
                                        Point(1, 0, 0), Point(1, 1, 2));

  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid()<<std::endl;
  constrained_delaunay_triangulation(lcc,d1);
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid()<<std::endl;
  lcc.clear();
  std::cout<<std::endl
           <<"###################################################### \n"
           <<std::endl;

  // Create one hexahedron.
  d1 = lcc.make_hexahedron(Point(0,0,0), Point(1,0,0), Point(1,1,0),
                           Point(0,1,0), Point(0,1,1), Point(0,0,1),
                           Point(1,0,1), Point(1,1,1));

  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid()<<std::endl;
  
  constrained_delaunay_triangulation(lcc,d1);
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid()<<std::endl;  
  
  constrained_delaunay_triangulation(lcc,d1->beta(2));
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid()<<std::endl;
  lcc.clear();
  std::cout<<std::endl
           <<"###################################################### \n"
           <<std::endl;
  
  std::vector<Point> points;
  points.push_back(Point(0,0,0));
  points.push_back(Point(5,15,0));
  points.push_back(Point(8,18,0));
  points.push_back(Point(12,5,0));
  points.push_back(Point(8,3,0));
  points.push_back(Point(8,-9,0));
  points.push_back(Point(5,0,0));
  points.push_back(Point(2,-3,2));
  d1=make_facet(lcc,points);

  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid()<<std::endl;
#ifdef CGAL_LCC_USE_VIEWER
  display_lcc(lcc);
#endif // CGAL_LCC_USE_VIEWER
  
  constrained_delaunay_triangulation(lcc,d1);
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid()<<std::endl;
#ifdef CGAL_LCC_USE_VIEWER
  display_lcc(lcc);
#endif // CGAL_LCC_USE_VIEWER

  lcc.clear();
  std::cout<<std::endl
           <<"###################################################### \n"
           <<std::endl;

  return EXIT_SUCCESS;
}

