//  (C) Copyright Gennadiy Rozental 2001-2004.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : defines algoirthms for comparing 2 floating point values
// ***************************************************************************

#ifndef BOOST_FLOATING_POINT_COMPARISON_HPP_071894GER
#define BOOST_FLOATING_POINT_COMPARISON_HPP_071894GER

#include <boost/limits.hpp>  // for std::numeric_limits

#include <boost/test/detail/class_properties.hpp>

namespace boost {

namespace test_tools {

using unit_test::readonly_property;

// ************************************************************************** //
// **************        floating_point_comparison_type        ************** //
// ************************************************************************** //

enum floating_point_comparison_type { FPC_STRONG, FPC_WEAK };

// ************************************************************************** //
// **************                    details                   ************** //
// ************************************************************************** //

namespace tt_detail {

template<typename FPT>
inline FPT
fpt_abs( FPT arg ) 
{
    return arg < 0 ? -arg : arg;
}

//____________________________________________________________________________//

// both f1 and f2 are unsigned here
template<typename FPT>
inline FPT 
safe_fpt_division( FPT f1, FPT f2 )
{
    return  (f2 < 1 && f1 > f2 * (std::numeric_limits<FPT>::max)())               ? (std::numeric_limits<FPT>::max)()
            : ((f2 > 1 && f1 < f2 * (std::numeric_limits<FPT>::min)() || f1 == 0) ? 0
                                                                                  : f1/f2 );
}

//____________________________________________________________________________//

} // namespace tt_detail

// ************************************************************************** //
// **************             close_at_tolerance               ************** //
// ************************************************************************** //

template<typename FPT, typename PersentType = FPT >
class close_at_tolerance {
public:
    explicit    close_at_tolerance( PersentType percentage_tolerance, floating_point_comparison_type fpc_type = FPC_STRONG ) 
    : p_fraction_tolerance( static_cast<FPT>(0.01)*percentage_tolerance ), p_strong_or_weak( fpc_type ==  FPC_STRONG ) {}

    bool        operator()( FPT left, FPT right ) const
    {
        FPT diff = tt_detail::fpt_abs( left - right );
        FPT d1   = tt_detail::safe_fpt_division( diff, tt_detail::fpt_abs( right ) );
        FPT d2   = tt_detail::safe_fpt_division( diff, tt_detail::fpt_abs( left ) );
        
        return p_strong_or_weak ? (d1 <= p_fraction_tolerance.get() && d2 <= p_fraction_tolerance.get()) 
                                : (d1 <= p_fraction_tolerance.get() || d2 <= p_fraction_tolerance.get());
    }

    // Public properties
    readonly_property<FPT>  p_fraction_tolerance;
    readonly_property<bool> p_strong_or_weak;
};

//____________________________________________________________________________//

// ************************************************************************** //
// **************               check_is_close                 ************** //
// ************************************************************************** //

template<typename FPT, typename PersentType>
inline bool
check_is_close( FPT left, FPT right, PersentType percentage_tolerance, floating_point_comparison_type fpc_type = FPC_STRONG )
{
    close_at_tolerance<FPT,PersentType> pred( percentage_tolerance, fpc_type );

    return pred( left, right );
}

//____________________________________________________________________________//

template<typename FPT>
inline FPT
compute_tolerance( FPT percentage_tolerance )
{
    close_at_tolerance<FPT> pred( percentage_tolerance );

    return pred.p_fraction_tolerance.get();
}

//____________________________________________________________________________//

} // namespace test_tools
} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1.1.2  2004/11/20 10:52:13  spion
//  Import of Boost v. 1.32.0
//
//  Revision 1.18  2004/07/19 12:14:09  rogeeff
//  guard rename
//  tolerance parameter renamed for clarity
//
//  Revision 1.17  2004/06/07 07:33:49  rogeeff
//  detail namespace renamed
//
//  Revision 1.16  2004/05/21 06:19:35  rogeeff
//  licence update
//
//  Revision 1.15  2004/05/11 11:00:34  rogeeff
//  basic_cstring introduced and used everywhere
//  class properties reworked
//
//  Revision 1.14  2004/02/26 18:26:57  eric_niebler
//  remove minmax hack from win32.hpp and fix all places that could be affected by the minmax macros
//
//  Revision 1.13  2003/12/01 00:41:56  rogeeff
//  prerelease cleaning
//
// ***************************************************************************

#endif // BOOST_FLOATING_POINT_COMAPARISON_HPP_071894GER
