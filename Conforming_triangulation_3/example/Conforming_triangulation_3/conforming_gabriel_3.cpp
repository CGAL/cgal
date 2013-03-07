#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Conforming_Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>


typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;

typedef CGAL::Triangulation_vertex_base_with_info_3<int,K,CGAL::Conforming_triangulation_vertex_base_3<K> > Vb;
typedef CGAL::Conforming_triangulation_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3 < Vb, Cb > Tds;

typedef CGAL::Conforming_Delaunay_triangulation_3< K, Tds > CDT3;


template <class Triangulation, class InputIterator>
std::ptrdiff_t insert_with_info(Triangulation& t, InputIterator first,InputIterator last)
{
  typedef typename Triangulation::Geom_traits Geom_traits;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  
  std::size_t n = t.number_of_vertices();
  std::vector<std::ptrdiff_t> indices;
  std::vector<Point_3> points;
  std::vector<typename Triangulation::Vertex::Info> infos;
  std::ptrdiff_t index=0;
  for (InputIterator it=first;it!=last;++it){
    std::pair<Point_3,int> value=*it;
    points.push_back( value.first  );
    infos.push_back ( value.second );
    indices.push_back(index++);
  }

  typedef CGAL::Spatial_sort_traits_adapter_3<Geom_traits,Point_3*> Search_traits;
  
  CGAL::spatial_sort(indices.begin(),indices.end(),Search_traits(&(points[0]),Geom_traits()));

  Vertex_handle hint;
  for (typename std::vector<std::ptrdiff_t>::const_iterator
    it = indices.begin(), end = indices.end();
    it != end; ++it){
    hint = t.insert(points[*it], hint);
    if (hint!=Vertex_handle()) hint->info()=infos[*it];
  }

  return t.number_of_vertices() - n;
}

int main()
{
  
  std::vector< std::pair<Point_3,int> > points;
  points.push_back( std::make_pair(Point_3(0,0,0),0) );
  points.push_back( std::make_pair(Point_3(1,0,0),1) );
  points.push_back( std::make_pair(Point_3(0,1,0),2) );
  points.push_back( std::make_pair(Point_3(0.2,0.2,100),3) );
  points.push_back( std::make_pair(Point_3(0.2,0.2,-100),4) );
  
  
  std::vector< std::pair<int, int> > constraints;
  constraints.push_back( std::make_pair(0,1) );
  constraints.push_back( std::make_pair(0,2) );
  constraints.push_back( std::make_pair(1,2) );
  constraints.push_back( std::make_pair(3,4) );
  
  CDT3 cdt3;
  insert_with_info( cdt3, points.begin(), points.end() );
  
  
  std::vector< CDT3::Vertex_handle > vertices( points.size() );

  for (CDT3::Finite_vertices_iterator vit=cdt3.finite_vertices_begin(), 
                                      vit_end=cdt3.finite_vertices_end();vit!=vit_end;++vit)
  {
    vertices[ vit->info() ]=vit;
  }

  std::cout << cdt3.number_of_vertices() << std::endl;
  
  std::cout << "Making Delaunay conform" << std::endl;
  for (std::vector< std::pair<int,int> >::iterator cst_it=constraints.begin(),
                                                  cst_end=constraints.end();
       cst_it!=cst_end; ++cst_it)
  {
    cdt3.insert_conforming( std::make_pair(vertices[cst_it->first], vertices[cst_it->second]) );
  }


  for(CDT3::Finite_edges_iterator it=cdt3.finite_edges_begin(); it!=cdt3.finite_edges_end(); ++it)
    if ( it->first->vertex(it->second)->steiner() || it->first->vertex(it->third)->steiner() )
      std::cout << "edge with steiner endpoint" << std::endl;

  std::cout << "Making Gabriel conform" << std::endl;
  for (std::vector< std::pair<int,int> >::iterator cst_it=constraints.begin(),
                                                  cst_end=constraints.end();
       cst_it!=cst_end; ++cst_it)
  {
    cdt3.insert_conforming_Gabriel( vertices[cst_it->first], vertices[cst_it->second] );
  }

  for(CDT3::Finite_edges_iterator it=cdt3.finite_edges_begin(); it!=cdt3.finite_edges_end(); ++it)
    if ( it->first->vertex(it->second)->steiner() || it->first->vertex(it->third)->steiner() )
      std::cout << "edge with steiner endpoint" << std::endl;
}

