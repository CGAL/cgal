#include <CGAL/Combinatorial_map.h>

struct Map_3_dart_items_3: public CGAL::Generic_map_min_items
{
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif
};

using Map3=CGAL::Combinatorial_map<3, Map_3_dart_items_3>;
Map3 m;

template<typename Range>
bool test_empty_it(const Range& r, const std::string& txt)
{
  bool res=true;
  if(r.size()!=0)
  {
    std::cout<<"[ERROR "<<txt<<"] non zero size with range: "<<r.size()<<std::endl;
    res=false;
  }
  std::size_t nb=0;
  for(auto it=r.begin(), itend=r.end(); it!=itend; ++it)
  { ++nb; }
  if(nb!=0)
  {
    std::cout<<"[ERROR "<<txt<<"] non zero size with it: "<<nb<<std::endl;
    res=false;
  }

  return res;
}
int main()
{
  bool res=true;

  res=res && test_empty_it<Map3::Dart_range>(m.darts(), "Dart_range");
  res=res && test_empty_it<Map3::Dart_const_range>
    (const_cast<const Map3&>(m).darts(), "Dart_const_range");

  res=res && test_empty_it<Map3::One_dart_per_cell_range<3>>
    (m.one_dart_per_cell<3>(), "One_dart_per_cell_range<0>");

  res=res && test_empty_it<Map3::One_dart_per_cell_const_range<3>>
    (const_cast<const Map3&>(m).one_dart_per_cell<3>(),
     "One_dart_per_cell_const_range<0>");

  if(!res)
  { return(EXIT_FAILURE); }

  std::cout<<"ALL SUCCESS."<<std::endl;
  return EXIT_SUCCESS;
}
