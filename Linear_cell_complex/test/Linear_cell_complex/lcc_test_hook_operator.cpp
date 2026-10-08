#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>

struct Items_index: public CGAL::Linear_cell_complex_min_items
{ using Use_index=CGAL::Tag_true; };

template<typename LCC>
bool test_lcc()
{
  LCC lcc;
  typename LCC::Dart_descriptor dd;
  for(int i=0; i<10; ++i)
  {
    typename LCC::Dart_descriptor cur=lcc.create_dart();
    if(lcc.dart_descriptor(lcc.darts()[lcc.darts().index(cur)])!=cur)
    { return false; }
  }

  return true;
}

int main()
{
  using LCC2=CGAL::Linear_cell_complex_for_combinatorial_map<2,2>;
  using LCC2_INDEX=CGAL::Linear_cell_complex_for_combinatorial_map
    <2,2,CGAL::Linear_cell_complex_traits<2>, Items_index>;

  if(!test_lcc<LCC2>() || !test_lcc<LCC2_INDEX>())
  {
    std::cout<<"ERROR lcc_test_hook_operator."<<std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
