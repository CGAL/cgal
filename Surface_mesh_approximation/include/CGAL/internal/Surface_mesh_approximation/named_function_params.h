#ifndef CGAL_VSA_BGL_NAMED_FUNCTION_PARAMS_H
#define CGAL_VSA_BGL_NAMED_FUNCTION_PARAMS_H

#include <CGAL/boost/graph/named_function_params.h>

namespace CGAL{
namespace internal_np{

// define enum types and values for new named parameters
#define CGAL_add_named_parameter(X, Y, Z)            \
  enum X { Y };
#include <CGAL/Surface_mesh_approximation/internal/parameters_interface.h>
#undef CGAL_add_named_parameter

}//internal_np

  template <typename T, typename Tag, typename Base = boost::no_property>
  struct vsa_bgl_named_params
    : CGAL::cgal_bgl_named_params<T, Tag, Base>
  {
    typedef CGAL::cgal_bgl_named_params<T, Tag, Base> base;
    typedef vsa_bgl_named_params self;

    vsa_bgl_named_params(T v = T()) : base(v) {}
    vsa_bgl_named_params(T v, const Base& b) : base(v, b) {}

// define the functions for new named parameters and the one imported from BGL and boost
// used to concatenate several parameters
#define CGAL_add_named_parameter(X, Y, Z)              \
  template<typename K>                               \
  vsa_bgl_named_params<K, internal_np::X, self>                   \
  Z(const K& k) const                                \
  {                                                  \
    typedef vsa_bgl_named_params<K, internal_np::X, self> Params; \
    return Params(k, *this);                         \
  }
#include <CGAL/Surface_mesh_approximation/internal/parameters_interface.h>
#include <CGAL/boost/graph/parameters_interface.h>
#include <CGAL/boost/graph/boost_parameters_interface.h>
#undef CGAL_add_named_parameter
  };

namespace VSA{

namespace parameters{

// define free functions for new named parameters and the one imported from BGL and boost
#define CGAL_add_named_parameter(X, Y, Z)          \
  template<typename K>                           \
  vsa_bgl_named_params<K, internal_np::X>                     \
  Z(const K& k)                                  \
  {                                              \
    typedef vsa_bgl_named_params<K, internal_np::X> Params;   \
    return Params(k);                            \
  }
#include <CGAL/boost/graph/boost_parameters_interface.h>
#include <CGAL/boost/graph/parameters_interface.h>
#include <CGAL/Surface_mesh_approximation/internal/parameters_interface.h>
#undef CGAL_add_named_parameter
} //namespace parameters
} //namespace VSA

} //namespace CGAL

// partial specializations hate inheritance and we need to repeat
// those here. this is rather fragile.
namespace boost {
#if BOOST_VERSION < 105100
  template <class Tag1, class Tag2, class T1, class Base>
  inline
  typename property_value< CGAL::vsa_bgl_named_params<T1,Tag1,Base>, Tag2>::type
  get_param(const CGAL::vsa_bgl_named_params<T1,Tag1,Base>& p, Tag2 tag2)
  {
    enum { match = detail::same_property<Tag1,Tag2>::value };
    typedef typename
      boost::property_value< CGAL::vsa_bgl_named_params<T1,Tag1,Base>, Tag2>::type T2;
    T2* t2 = 0;
    typedef detail::property_value_dispatch<match> Dispatcher;
    return Dispatcher::const_get_value(p, t2, tag2);
  }
#endif

  template <typename T, typename Tag, typename Base, typename Def>
  struct lookup_named_param_def<Tag, CGAL::vsa_bgl_named_params<T, Tag, Base>, Def> {
    typedef T type;
    static const type& get(const bgl_named_params<T, Tag, Base>& p, const Def&) {
      return p.m_value;
    }
  };

  template <typename Tag1, typename T, typename Tag, typename Base, typename Def>
  struct lookup_named_param_def<Tag1, CGAL::vsa_bgl_named_params<T, Tag, Base>, Def> {
    typedef typename lookup_named_param_def<Tag1, Base, Def>::type type;
    static const type& get(const bgl_named_params<T, Tag, Base>& p, const Def& def) {
      return lookup_named_param_def<Tag1, Base, Def>::get(p.m_base, def);
    }
  };

} // boost

#endif //CGAL_VSA_BGL_NAMED_FUNCTION_PARAMS_H
