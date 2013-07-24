#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Timer.h>

#include <cassert>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;


// usage example CDT_insert_on_constraints < norway.cin
int main()
{
  std::size_t nbs;
  char s;
  CDT::Point p1, p2;
  std::cin >> nbs;

  CDT cdt;

  do
  {
    std::cin >> s >> p1 >> p2;
    cdt.insert_constraint(p1,p2);
  }
  while ( --nbs!=0 );

  std::cout << "Triangulation built "<< cdt.number_of_vertices() << "\n";

  {
    CDT cdt2(cdt);
    CGAL::Timer time;
    time.start();
    std::vector< std::pair<CDT::Vertex_handle, CDT::Vertex_handle> >csts;
    for ( CDT::Finite_edges_iterator eit=cdt2.finite_edges_begin(),
                                     eit_end=cdt2.finite_edges_end();
                                     eit!=eit_end; ++eit)
    {
      if ( cdt2.is_constrained(*eit) )
      {
        CDT::Vertex_handle v1 = eit->first->vertex( CDT::cw(eit->second) );
        CDT::Vertex_handle v2 = eit->first->vertex( CDT::ccw(eit->second) );
        csts.push_back( std::make_pair(v1,v2) );
      }
    }

    std::cout << csts.size() << " edges to process\n";
    while ( !csts.empty() )
    {
      CDT::Face_handle face;
      int index=0;
      CDT::Vertex_handle v1=csts.back().first;
      CDT::Vertex_handle v2=csts.back().second;
      csts.pop_back();
      cdt2.is_edge( v1, v2, face, index );
      cdt2.remove_constrained_edge(face, index);
      CDT::Vertex_handle vn =
        cdt2.insert( CGAL::midpoint(v1->point(), v2->point()) );
      cdt2.insert_constraint(vn,v1);
      cdt2.insert_constraint(vn,v2);
    }

    time.stop();
    std::cout << "Unconstraining first " << cdt2.number_of_vertices() << " " << time.time() << std::endl;
  }

  {
    CDT cdt2(cdt);
    CGAL::Timer time;
    time.start();
    std::vector< std::pair<CDT::Vertex_handle, CDT::Vertex_handle> >csts;
    for ( CDT::Finite_edges_iterator eit=cdt2.finite_edges_begin(),
                                     eit_end=cdt2.finite_edges_end();
                                     eit!=eit_end; ++eit)
    {
      if ( cdt2.is_constrained(*eit) )
      {
        CDT::Vertex_handle v1 = eit->first->vertex( CDT::cw(eit->second) );
        CDT::Vertex_handle v2 = eit->first->vertex( CDT::ccw(eit->second) );
        csts.push_back( std::make_pair(v1,v2) );
      }
    }

    std::cout << csts.size() << " edges to process\n";
    while ( !csts.empty() )
    {
      CDT::Face_handle face;
      int index=0;
      CDT::Vertex_handle v1=csts.back().first;
      CDT::Vertex_handle v2=csts.back().second;
      csts.pop_back();
      cdt2.is_edge( v1, v2, face, index );
      cdt2.remove_constrained_edge(face, index);
      CDT::Vertex_handle vn =
        cdt2.insert( CGAL::midpoint(v1->point(), v2->point()), v1->face() );
      cdt2.insert_constraint(vn,v1);
      cdt2.insert_constraint(vn,v2);
    }

    time.stop();
    std::cout << "Unconstraining first +hint " << cdt2.number_of_vertices() << " " << time.time() << std::endl;
  }

  {
    CDT cdt2(cdt);
    CGAL::Timer time;
    time.start();
    std::vector< std::pair<CDT::Vertex_handle, CDT::Vertex_handle> >csts;
    for ( CDT::Finite_edges_iterator eit=cdt2.finite_edges_begin(),
                                     eit_end=cdt2.finite_edges_end();
                                     eit!=eit_end; ++eit)
    {
      if ( cdt2.is_constrained(*eit) )
      {
        CDT::Vertex_handle v1 = eit->first->vertex( CDT::cw(eit->second) );
        CDT::Vertex_handle v2 = eit->first->vertex( CDT::ccw(eit->second) );
        csts.push_back( std::make_pair(v1,v2) );
      }
    }

    std::cout << csts.size() << " edges to process\n";
    while ( !csts.empty() )
    {
      CDT::Vertex_handle v1=csts.back().first;
      CDT::Vertex_handle v2=csts.back().second;
      csts.pop_back();
      CDT::Vertex_handle vn =
        cdt2.insert( CGAL::midpoint(v1->point(), v2->point()) );
    }

    time.stop();
    std::cout << "Not unconstraining first " << cdt2.number_of_vertices() << " " << time.time() << std::endl;
  }

  {
    CDT cdt2(cdt);
    CGAL::Timer time;
    time.start();
    std::vector< std::pair<CDT::Vertex_handle, CDT::Vertex_handle> >csts;
    for ( CDT::Finite_edges_iterator eit=cdt2.finite_edges_begin(),
                                     eit_end=cdt2.finite_edges_end();
                                     eit!=eit_end; ++eit)
    {
      if ( cdt2.is_constrained(*eit) )
      {
        CDT::Vertex_handle v1 = eit->first->vertex( CDT::cw(eit->second) );
        CDT::Vertex_handle v2 = eit->first->vertex( CDT::ccw(eit->second) );
        csts.push_back( std::make_pair(v1,v2) );
      }
    }

    std::cout << csts.size() << " edges to process\n";
    while ( !csts.empty() )
    {
      CDT::Face_handle face;
      int index;
      CDT::Vertex_handle v1=csts.back().first;
      CDT::Vertex_handle v2=csts.back().second;
      csts.pop_back();
      cdt2.is_edge( v1, v2, face, index );
      CDT::Vertex_handle vn =
        cdt2.insert( CGAL::midpoint(v1->point(), v2->point()), face );
    }

    time.stop();
    std::cout << "Not unconstraining + hint " << cdt2.number_of_vertices() << " " << time.time() << std::endl;
  }

  {
    CDT cdt2(cdt);
    CGAL::Timer time;
    time.start();
    std::vector< std::pair<CDT::Vertex_handle, CDT::Vertex_handle> >csts;
    for ( CDT::Finite_edges_iterator eit=cdt2.finite_edges_begin(),
                                     eit_end=cdt2.finite_edges_end();
                                     eit!=eit_end; ++eit)
    {
      if ( cdt2.is_constrained(*eit) )
      {
        CDT::Vertex_handle v1 = eit->first->vertex( CDT::cw(eit->second) );
        CDT::Vertex_handle v2 = eit->first->vertex( CDT::ccw(eit->second) );
        csts.push_back( std::make_pair(v1,v2) );
      }
    }

    std::cout << csts.size() << " edges to process\n";
    while ( !csts.empty() )
    {
      CDT::Face_handle face;
      int index=0;
      CDT::Vertex_handle v1=csts.back().first;
      CDT::Vertex_handle v2=csts.back().second;
      csts.pop_back();
      cdt2.is_edge( v1, v2, face, index );
      CDT::Vertex_handle vn =
        cdt2.insert( CGAL::midpoint(v1->point(), v2->point()), CDT::EDGE, face, index );
    }

    time.stop();
    std::cout << "Not unconstraining + locate " << cdt2.number_of_vertices() << " " << time.time() << std::endl;
  }
}
