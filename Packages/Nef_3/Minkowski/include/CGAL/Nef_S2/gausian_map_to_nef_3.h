
#include <CGAL/Modifier_base.h>

CGAL_BEGIN_NAMESPACE

template<typename Kernel_, typename Items_, typename Mark_>
  class gausian_map_to_nef_3 : public Modifier_base<SNC_structure<Kernel_,Items_,Mark_> {

  typedef Kernel_                                 Kernel;
  typedef Items_                                  Items;
  typedef Mark_                                   Mark;
  typedef CGAL::SNC_structure<Kernel,Items,Mark>  SNC_structure;
  typedef CGAL::Gausian_map<Kernel>               Gausian_map;

  Gausian_map& G;

 public:
  gausian_map_to_nef_3(Gausian_map& Gin) : G(Gin) {}

  void operator()(SNC_structure& snc) {
    
    polyhedron_3_to_nef_3<Gausian_map, SNC_structure>(G, snc);    
  }






};

CGAL_END_NAMESPACE
