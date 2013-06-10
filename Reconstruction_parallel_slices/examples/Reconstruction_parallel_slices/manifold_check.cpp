#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/tuple.h>
#include <fstream>
#include <set>
#include <iostream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;

template <class PointOutputIterator>
void read_points(std::size_t n,std::ifstream& input,PointOutputIterator out){
  Kernel::Point_3 p;
  for (;n>0;--n){
    input >> p;
    *out++=p;
  }
}

template <class TripleOutputIterator>
void read_facets(std::size_t n,std::ifstream& input,TripleOutputIterator out){
  int three,i,j,k;
  for (;n>0;--n){
    input >> three >> i >> j >> k;
    CGAL_assertion(three==3);
    *out++=CGAL::cpp0x::make_tuple(i,j,k);
  }
}

std::pair<int,int> make_sorted_pair(int i,int j){
  if (i < j) return std::make_pair(i,j);
  return std::make_pair(j,i);
}

typedef CGAL::cpp0x::tuple<int,int,int> Facet;
typedef std::list<Facet> Facets;
typedef std::map<std::pair<int,int>,int> Edge_map;

bool keep_facets_edge_criterium(int i,int j,Edge_map& edge_map){
  Edge_map::iterator it = edge_map.insert( std::make_pair( make_sorted_pair(i,j),0 ) ).first;
  if (it->second==2) return false;
  it->second+=1;
  return true;
}

int main(int,char** argv)
{
  std::ifstream input(argv[1]);
  
  /*
  std::string OFF;
  input >> OFF;
  CGAL_assertion(OFF=="OFF");
  std::size_t nbv,nbf,nbt;
  input >> nbv >> nbf >> nbt;
  
  
  std::vector<Kernel::Point_3> points;
  points.reserve(nbv);
  
  Facets facets;
  
  
  read_points(nbv,input,std::back_inserter(points));
  read_facets(nbf,input,std::back_inserter(facets));
 
  Edge_map edge_map;
  Facets removed_facets;

//remove facets incident to an edge with more than two incident facets
  for (Facets::iterator it=facets.begin();it!=facets.end();){
    int indices[3]={CGAL::cpp0x::get<0>(*it),CGAL::cpp0x::get<1>(*it),CGAL::cpp0x::get<2>(*it)};
    int i=0;
    for (;i<3;++i){
      if (! keep_facets_edge_criterium( indices[(i+1)%3],indices[(i+2)%3],edge_map ) ){
        Facets::iterator oit=it;
        ++it;
        removed_facets.splice(removed_facets.end(),facets,oit);
        break;
      }
    }
    if (i==3) ++it;
  }

  typedef std::map<int,std::list<std::pair<Facets::iterator,int> > > Central_vertex_to_edges;
  Central_vertex_to_edges central_vertex_to_edges;
  
//remove facets that are incident to a common vertex but a the link of the facets is not connected
  bool facets_removed=true;
  while (facets_removed){
    facets_removed=false;
    //collect edges of the link of each vertex
    for (Facets::iterator it=facets.begin();it!=facets.end();++it){
      int indices[3]={CGAL::cpp0x::get<0>(*it),CGAL::cpp0x::get<1>(*it),CGAL::cpp0x::get<2>(*it)};
      for (int i=0;i<3;++i){
        Central_vertex_to_edges::iterator mit=
          central_vertex_to_edges.insert( std::make_pair(indices[i],std::list<std::pair<Facets::iterator,int> >() ) ).first;
        mit->second.push_back(std::make_pair(it,i));
      }
    }
  
    for (Central_vertex_to_edges::iterator it=central_vertex_to_edges.begin(),end=central_vertex_to_edges.end();it!=end;++it){
      int index=it->first;
      std::list<std::pair<Facets::iterator,int> >& edges=it->second;
      typedef std::map<int,std::list<Facets::iterator> > Contributing_facets;
      Contributing_facets contributing_facets;
      for (std::list<std::pair<Facets::iterator,int> >::iterator eit=edges.begin();eit!=edges.end();++eit)
      {
        int indices[3]={CGAL::cpp0x::get<0>(*(eit->first)),CGAL::cpp0x::get<1>(*(eit->first)),CGAL::cpp0x::get<2>(*(eit->first))};
        CGAL_assertion(index==indices[eit->second]);
        
        contributing_facets
          .insert(std::make_pair(indices[(eit->second+1)%3],std::list<Facets::iterator>()))
              .first->second.push_back(eit->first);
        contributing_facets
          .insert(std::make_pair(indices[(eit->second+2)%3],std::list<Facets::iterator>()))
              .first->second.push_back(eit->first);
      }
      
      int nb_one=0;
      for (Contributing_facets::iterator cit=contributing_facets.begin();cit!=contributing_facets.end();++cit){
        std::size_t size_of_list=cit->second.size();
        if (size_of_list==1) ++nb_one;
        CGAL_assertion(size_of_list==2 || size_of_list==1);
      }
      
      if (nb_one > 2){
        for( std::list<std::pair<Facets::iterator,int> >::iterator eit=edges.begin();eit!=edges.end();++eit )
          removed_facets.splice(removed_facets.end(),facets,eit->first);
        facets_removed=true;
        std::cout << "Remove facets around vertex " << index << std::endl;
        break;
      }
    }
  }
  
  input.close();
  input.open(argv[1]);
  */
  
  Polyhedron_3 p;
  CGAL::scan_OFF( input, p,true);
  std::cout << p.size_of_vertices() << std::endl;
  if ( !p.is_closed () ) std::cout << "WARNING the polyhedron is not closed!!!!" << std::endl;
  
  //check that there is not degenerate triangles
  std::cout << "Checking for degenerate triangles\n";
  int degen_triangle=0;
  for (Polyhedron_3::Facet_iterator it=p.facets_begin();it!=p.facets_end();++it)
  {
    std::set<Kernel::Point_3> points;
    points.insert( it->halfedge()->vertex()->point() );
    points.insert( it->halfedge()->next()->vertex()->point() );
    points.insert( it->halfedge()->opposite()->vertex()->point() );
    if (points.size()!=3) ++degen_triangle;
  }
  if (degen_triangle!=0) std::cout << degen_triangle << " degenerate triangles" << std::endl;
  
}
  

