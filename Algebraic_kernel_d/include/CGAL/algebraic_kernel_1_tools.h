// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
// 
// Author(s)     : Michael Hemmer <hemmer@mpi-ing.mpg.de>
//
// ============================================================================
//

// The tools inside this file are not documented. 
// However, they may find their way into the kernel. 

// TODO: is there a more elaborated way to get the first non negative root ?
// TODO: is there a more elabotated way to compare the first non negative roots ?

#ifndef CGAL_ALGEBRAIC_KERNEL_1_TOOLS_H
#define CGAL_ALGEBRAIC_KERNEL_1_TOOLS_H

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_1.h>
#include <CGAL/Polynomial_traits_d.h>
#include <vector>

#include <boost/optional.hpp>

CGAL_BEGIN_NAMESPACE

// returns the first nonnegative root of p
// precondition: p has at least one non negative root
template <class Algebraic_kernel_1>
boost::optional< typename Algebraic_kernel_1::Algebraic_real_1 >
compute_smallest_nonnegative_root(
    const Algebraic_kernel_1& ak,
    const typename Algebraic_kernel_1::Polynomial_1& p){
  
  typedef Algebraic_kernel_1 AK;
  typedef typename AK::Algebraic_real_1 Root;
  typedef boost::optional< Root > Root_option; 
  
  typename AK::Solve_1 solve_1 = ak.solve_1_object();
  std::vector<Root> roots;
  
  solve_1(p, std::back_inserter(roots));
  typename std::vector<Root>::const_iterator it = roots.begin();
  for (; it != roots.end(); ++it)
  {
    if(CGAL::is_negative(*it) == false)
    {
      break;
    }
  }

  if (it == roots.end())
    return Root_option();
  else 
    return Root_option(*it);
}

template <class Algebraic_kernel_1>
CGAL::Comparison_result
compare_smallest_nonnegative_roots(
    const Algebraic_kernel_1& ak,
    const typename Algebraic_kernel_1::Polynomial_1& p1,
    const typename Algebraic_kernel_1::Polynomial_1& p2){
  typedef  typename Algebraic_kernel_1::Algebraic_real_1 Root; 
  typedef boost::optional<Root> Root_option; 
  Root_option r1  = compute_smallest_nonnegative_root(ak,p1);
  Root_option r2  = compute_smallest_nonnegative_root(ak,p2);
  if(!r1 && !r2 ) return CGAL::EQUAL;
  if(!r1)         return CGAL::LARGER;
  if(!r2)         return CGAL::SMALLER; 
  return                 CGAL::compare(*r1,*r2);
}

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_1_TOOLS_H
