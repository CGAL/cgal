#include <CGAL/Combinatorial_map.h>

struct Map_3_dart_items_3: public CGAL::Generic_map_min_items
{
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif
};

using Map3=CGAL::Combinatorial_map<3, Map_3_dart_items_3>;

template<typename Range>
bool test_empty_it(Range& r, const std::string& txt)
{
  bool res=true;
  Map3 m;
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

}
