#ifndef CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H
#define CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H

#include <vector>

#include <CGAL/Cartesian.h>
//#include <CGAL/Interval_arithmetic.h>
//#include <CGAL/Arithmetic_filter.h>
#include <CGAL/Filtered_exact.h>

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Segment_Voronoi_diagram_site_2.h>

#include <CGAL/predicates/Segment_Voronoi_diagram_predicates_ftC2.h>

CGAL_BEGIN_NAMESPACE

static unsigned int num_failures_side_of_bisector = 0;
static unsigned int num_failures_vertex_conflict = 0;
static unsigned int num_failures_finite_edge_conflict = 0;
static unsigned int num_failures_infinite_edge_conflict = 0;

template < class CT, class ET, bool Protected, class Cache, class Method_tag >
inline
Comparison_result
svd_compare_distanceC2(const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& qx,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& qy,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& sx,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& sy,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& tx,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& ty,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& px,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& py,
		       Method_tag method_tag)
{
  //  do_not_compile();

  try {
    Protect_FPU_rounding<Protected> Protection;

    Interval_nt_advanced qx_it = qx.interval();
    Interval_nt_advanced qy_it = qy.interval();
    Interval_nt_advanced sx_it = sx.interval();
    Interval_nt_advanced sy_it = sy.interval();
    Interval_nt_advanced tx_it = tx.interval();
    Interval_nt_advanced ty_it = ty.interval();
    Interval_nt_advanced px_it = px.interval();
    Interval_nt_advanced py_it = py.interval();

    return svd_compare_distanceC2(qx_it, qy_it, sx_it, sy_it,
				  tx_it, ty_it, px_it, py_it,
				  method_tag);
  }
  catch (Interval_nt_advanced::unsafe_comparison) {
    Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

    num_failures_side_of_bisector++;

    //    std::cerr << "inside svd_compare_distanceC2 #1 after exception"
    //    	      << std::endl;

    ET qx_et = qx.exact();
    ET qy_et = qy.exact();
    ET sx_et = sx.exact();
    ET sy_et = sy.exact();
    ET tx_et = tx.exact();
    ET ty_et = ty.exact();
    ET px_et = px.exact();
    ET py_et = py.exact();

    return svd_compare_distanceC2(qx_et, qy_et, sx_et, sy_et,
				  tx_et, ty_et, px_et, py_et,
				  method_tag);

  }
}



template < class CT, class ET, bool Protected, class Cache, class Method_tag >
inline
Comparison_result
svd_compare_distanceC2(const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& qx,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& qy,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& s1x,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& s1y,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& t1x,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& t1y,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& s2x,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& s2y,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& t2x,
		       const Filtered_exact<CT,ET,Dynamic,Protected,Cache>& t2y,
		       Method_tag method_tag)
{
  //  do_not_compile();

  try {
    Protect_FPU_rounding<Protected> Protection;

    Interval_nt_advanced qx_it = qx.interval();
    Interval_nt_advanced qy_it = qy.interval();
    Interval_nt_advanced s1x_it = s1x.interval();
    Interval_nt_advanced s1y_it = s1y.interval();
    Interval_nt_advanced t1x_it = t1x.interval();
    Interval_nt_advanced t1y_it = t1y.interval();
    Interval_nt_advanced s2x_it = s2x.interval();
    Interval_nt_advanced s2y_it = s2y.interval();
    Interval_nt_advanced t2x_it = t2x.interval();
    Interval_nt_advanced t2y_it = t2y.interval();

    return svd_compare_distanceC2(qx_it, qy_it,
				  s1x_it, s1y_it, t1x_it, t1y_it,
				  s2x_it, s2y_it, t2x_it, t2y_it,
				  method_tag);
  }
  catch (Interval_nt_advanced::unsafe_comparison) {
    Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

    //    std::cerr << "inside svd_compare_distanceC2 #2 after exception"
    //    	      << std::endl;

    num_failures_side_of_bisector++;

    ET qx_et = qx.exact();
    ET qy_et = qy.exact();
    ET s1x_et = s1x.exact();
    ET s1y_et = s1y.exact();
    ET t1x_et = t1x.exact();
    ET t1y_et = t1y.exact();
    ET s2x_et = s2x.exact();
    ET s2y_et = s2y.exact();
    ET t2x_et = t2x.exact();
    ET t2y_et = t2y.exact();

#if 0
    std::cerr << "q: " << CGAL::to_double(qx_et) << " "
	      << CGAL::to_double(qy_et) << std::endl;
    std::cerr << "S1: " << CGAL::to_double(s1x_et) << " "
	      << CGAL::to_double(s1y_et) << ", "
	      << CGAL::to_double(t1x_et) << " "
	      << CGAL::to_double(t1y_et) << std::endl;
    std::cerr << "S2: " << CGAL::to_double(s2x_et) << " "
	      << CGAL::to_double(s2y_et) << ", "
	      << CGAL::to_double(t2x_et) << " "
	      << CGAL::to_double(t2y_et) << std::endl;
#endif

    return svd_compare_distanceC2(qx_et, qy_et,
				  s1x_et, s1y_et, t1x_et, t1y_et,
				  s2x_et, s2y_et, t2x_et, t2y_et,
				  method_tag);

  }
}

//----------------------------------------------------------------------------

template < class CT, class ET, bool Protected, class Cache, class Method_tag >
inline
Sign
svd_vertex_conflict_ftC2(const std::vector< 
			 Filtered_exact<CT,ET,Dynamic,Protected,Cache> >& v,
			 char site_types[], int num_sites,
			 Method_tag method_tag)
{
  try {
    Protect_FPU_rounding<Protected> Protection;

    std::vector<Interval_nt_advanced> v_IT(v.size());

    for (unsigned int i = 0; i < v.size(); i++) {
      v_IT[i] = v[i].interval();
    }
    return svd_vertex_conflict_ftC2(v_IT, site_types, num_sites,
				    method_tag);
  }
  catch (Interval_nt_advanced::unsafe_comparison) {
    Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

    num_failures_vertex_conflict++;

    //    std::cerr << "inside svd_vertex_conflict_ftC2 after exception"
    //	      << std::endl;

#if 0
    {
      int k = 0;
      for (int i = 0; i < num_sites; i++) {
	if ( site_types[i] == 'p' ) {
	  std::cerr << "p " << CGAL::to_double(v[k].exact())
		    << " " << CGAL::to_double(v[k+1].exact())
		    << std::endl;
	  k += 2;
	} else {
	  std::cerr << "s " << CGAL::to_double(v[k].exact())
		    << " " << CGAL::to_double(v[k+1].exact())
		    << " " << CGAL::to_double(v[k+2].exact())
		    << " " << CGAL::to_double(v[k+3].exact())
		    << std::endl;
	  k += 4;
	}
      }
    }
#endif

    std::vector<ET> v_ET(v.size());

    for (unsigned int i = 0; i < v.size(); i++) {
      v_ET[i] = v[i].exact();
    }

    return svd_vertex_conflict_ftC2(v_ET, site_types, num_sites,
				    method_tag);
  }
}


//----------------------------------------------------------------------------


template < class CT, class ET, bool Protected, class Cache, class Method_tag >
inline
bool
svd_finite_edge_conflict_ftC2(const std::vector< 
			      Filtered_exact<CT,ET,Dynamic,
			      Protected,Cache> >& v,
			      Sign sgn, char site_types[], int num_sites,
			      Method_tag method_tag)
{
  try {
    Protect_FPU_rounding<Protected> Protection;

    std::vector<Interval_nt_advanced> v_IT(v.size());

    for (unsigned int i = 0; i < v.size(); i++) {
      v_IT[i] = v[i].interval();
    }
    return svd_finite_edge_conflict_ftC2(v_IT, sgn, site_types,
					 num_sites, method_tag);
  }
  catch (Interval_nt_advanced::unsafe_comparison) {
    Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

    num_failures_finite_edge_conflict++;

    std::vector<ET> v_ET(v.size());

    for (unsigned int i = 0; i < v.size(); i++) {
      v_ET[i] = v[i].exact();
    }

    return svd_finite_edge_conflict_ftC2(v_ET, sgn, site_types,
					 num_sites, method_tag);
  }
}


//----------------------------------------------------------------------------


template < class CT, class ET, bool Protected, class Cache, class Method_tag >
inline
bool
svd_infinite_edge_conflict_ftC2(const std::vector< 
				Filtered_exact<CT,ET,Dynamic,
				Protected,Cache> >& v,
				Sign sgn, char site_types[], int num_sites,
				Method_tag method_tag)
{
  try {
    Protect_FPU_rounding<Protected> Protection;

    std::vector<Interval_nt_advanced> v_IT(v.size());

    for (unsigned int i = 0; i < v.size(); i++) {
      v_IT[i] = v[i].interval();
    }
    return svd_infinite_edge_conflict_ftC2(v_IT, sgn, site_types,
					   num_sites, method_tag);
  }
  catch (Interval_nt_advanced::unsafe_comparison) {
    Protect_FPU_rounding<!Protected> Protection(CGAL_FE_TONEAREST);

    num_failures_infinite_edge_conflict++;

    std::vector<ET> v_ET(v.size());

    for (unsigned int i = 0; i < v.size(); i++) {
      v_ET[i] = v[i].exact();
    }

    return svd_infinite_edge_conflict_ftC2(v_ET, sgn, site_types,
					   num_sites, method_tag);
  }
}


//----------------------------------------------------------------------------




CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H
