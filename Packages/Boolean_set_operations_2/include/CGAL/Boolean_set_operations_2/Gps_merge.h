
#ifndef GPS_MERGE_H
#define GPS_MERGE_h

#include <CGAL/Boolean_set_operations_2/Gps_agg_op.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_join_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_xor_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_bfs_intersection_visitor.h>
#include <CGAL/Boolean_set_operations_2/Gps_agg_op.h>
#include <vector>

CGAL_BEGIN_NAMESPACE

template <class Arrangement>
class Join_merge
{
public:
   void operator()(unsigned int i,
                   unsigned int j,
                   unsigned int jump,
                   std::vector<Arrangement*>& arr_vec)                     
  {
    if(i==j)
      return;

    typename Arrangement::Traits_2*  tr = arr_vec[i]->get_traits();
    Arrangement* res = new Arrangement(tr);
    Gps_agg_op<Arrangement, Gps_bfs_join_visitor<Arrangement> >
      agg_op(*res, *tr);
    agg_op.sweep_arrangements(i, j, jump, arr_vec);

    for(unsigned int count=i; count<=j; count+=jump)
    {
      delete arr_vec[count];
    }
    
    arr_vec[i] = res;
  }

};


template <class Arrangement>
class Intersection_merge
{
public:
   void operator()(unsigned int i,
                   unsigned int j,
                   unsigned int jump,
                   std::vector<Arrangement*>& arr_vec)                     
  {
    if(i==j)
      return;

    typename Arrangement::Traits_2*  tr = arr_vec[i]->get_traits();
    Arrangement* res = new Arrangement(tr);
    Gps_agg_op<Arrangement, Gps_bfs_intersection_visitor<Arrangement> >
      agg_op(*res, *tr);
    agg_op.sweep_arrangements(i, j, jump, arr_vec);

    for(unsigned int count=i; count<=j; count+=jump)
    {
      delete arr_vec[count];
    }
    
    arr_vec[i] = res;
  }
};

template <class Arrangement>
class Xor_merge
{
public:
   void operator()(unsigned int i,
                   unsigned int j,
                   unsigned int jump,
                   std::vector<Arrangement*>& arr_vec)                     
  {
    if(i==j)
      return;

    typename Arrangement::Traits_2*  tr = arr_vec[i]->get_traits();
    Arrangement* res = new Arrangement(tr);
    Gps_agg_op<Arrangement, Gps_bfs_xor_visitor<Arrangement> >
      agg_op(*res, *tr);
    agg_op.sweep_arrangements(i, j, jump, arr_vec);

    for(unsigned int count=i; count<=j; count+=jump)
    {
      delete arr_vec[count];
    }
    
    arr_vec[i] = res;
  }
};

CGAL_END_NAMESPACE

#endif
