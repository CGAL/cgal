// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : /test/Triangulation/test_triangulation_2.C
// package       : Triangulation
// source        : $URL$
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Herve Bronnimann, Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
// ============================================================================

#include <utility>
#include <list>

#include <CGAL/_test_types.h>

#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_2.h>

#include <CGAL/_test_traits.h>
#include <CGAL/_test_cls_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

// Explicit instantiation of the whole class :
template class CGAL::Triangulation_2<TestK>;

template<typename Vb>
struct Custom_vertex_base : public Vb
{
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef Custom_vertex_base<Vb2>           Other;
  };
  
  std::ostream& write_data(std::ostream& os,
                           const CGAL::Unique_hash_map<typename Vb::Vertex_handle, int > &
                           )const 
  {
    os << this->info() << std::endl;
    return os;
  }
  
  std::istream& read_data(std::istream& is,
                          std::vector<typename Vb::Face_handle>&,
                          std::vector< typename Vb::Vertex_handle>&)
  {
    is >> this->info();
    return is;
  }
};

template<typename Fb>
struct Custom_face_base : public Fb
{
  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Fb::template Rebind_TDS<TDS2>::Other       Fb2;
    typedef Custom_face_base<Fb2>  Other;
  };
  Custom_face_base():Fb(){}
  Custom_face_base(typename Fb::Vertex_handle v0, typename Fb::Vertex_handle v1,
                   typename Fb::Vertex_handle v2):Fb(v0, v1, v2){}
  Custom_face_base(typename Fb::Vertex_handle v0, typename Fb::Vertex_handle v1, 
                   typename Fb::Vertex_handle v2, typename Fb::Face_handle n0, 
                   typename Fb::Face_handle n1, typename Fb::Face_handle n2)
    :Fb(v0,v1,v2,n0,n1,n2){}
  
  std::ostream& write_data(std::ostream& os,
                           const CGAL::Unique_hash_map<typename Fb::Vertex_handle, 
                           int >&)const
  {
    os << this->info() << std::endl;
    return os;
  }
  
  std::istream& read_data(std::istream& is,
                          std::vector<typename Fb::Face_handle>&,
                          std::vector<typename Fb::Vertex_handle>&)
  {
    is >> this->info();
    return is;
  }
};

int main()
{

  std::cout << "Testing Triangulation_2 " << std::endl; 
  std::cout << " with Cartesian : " << std::endl ;
  typedef Test_rep_cartesian Gt1;
  typedef CGAL::Triangulation_vertex_base_2<Gt1>                     Vb1;
  typedef CGAL::Triangulation_face_base_2<Gt1>                       Fb1;
  typedef CGAL::Triangulation_data_structure_2<Vb1,Fb1> Tds1;
  typedef CGAL::Triangulation_2<Gt1,Tds1>                            Cls1;
   _test_cls_triangulation_2( Cls1() );


  std::cout << std::endl << "Testing Triangulation_2" << std::endl;
  std::cout << " using Homogeneous  Kernel traits : " << std::endl;
  std::cout << " and defaults setting " << std::endl;
  typedef CGAL::Homogeneous<Rtype>      Gt6;
  typedef CGAL::Triangulation_2<Gt6>    Cls6;
  _test_cls_triangulation_2( Cls6() );
  
  typedef Custom_vertex_base<CGAL::Triangulation_vertex_base_with_info_2<int, Gt6> > Vertex_base;
  typedef Custom_face_base<CGAL::Triangulation_face_base_with_info_2<int, Gt6> > Cell_base;
  
  typedef CGAL::Triangulation_2<Gt6, CGAL::Triangulation_data_structure_2<Vertex_base, Cell_base > > Tr;
  typedef typename Tr::Point                Point;
  typedef typename Tr::Vertex_handle Vertex_handle;

  Tr TM_0, TM_1;
  Vertex_handle v0 = TM_0.insert(Point(0,0));
  Vertex_handle v1 = TM_0.insert(Point(1,0));
  Vertex_handle v2 = TM_0.insert(Point(2,0));
  Vertex_handle v3 = TM_0.insert(Point(1,1));
  
  v0->info() = 0;
  v1->info() = 1;
  v2->info() = 2;
  v3->info() = 3;
  
  int i=0;
  for(auto it = TM_0.tds().faces_begin();
      it != TM_0.tds().faces_end(); ++it)
  {
    it->info() = i++;
  }
  
  assert(TM_0.dimension() == 2);
  std::ofstream ofs("triangulation_output");
  ofs << TM_0;
  ofs.close();
  std::ifstream ifs("triangulation_output");
  ifs >> TM_1;
  ifs.close();
  auto it2 = TM_1.tds().faces_begin();
  
  for(auto it = TM_0.tds().faces_begin();
      it != TM_0.tds().faces_end(); ++it)
  {
    assert(it->info() == it2->info());
    ++it2;
  }
  
  
  return 0;
}
