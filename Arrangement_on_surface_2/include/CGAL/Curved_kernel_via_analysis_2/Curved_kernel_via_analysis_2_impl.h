// Copyright (c) 2007,2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_IMPL_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_IMPL_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*!\file include/CGAL/Curved_kernel_via_analysis_2.h
 * \brief defines class \c Curved_kernel_via_analysis_2
 *  
 * Defines points and arcs supported by curves that can be analyzed.
 */

#include <CGAL/config.h>

#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Curved_kernel_via_analysis_2/Point_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Arc_2.h>

#include <CGAL/Curved_kernel_via_analysis_2/Non_x_monotone_arc_2.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curve_interval_arcno_cache.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h>

namespace CGAL {

namespace internal {

/*!\brief
 * Provides basic types for Curved_kernel_via_analysis_2
 */
template < class NewCKvA, class BaseCKvA, class CurveKernel_2,
        template <class, class> class FunctorBase =
             Curved_kernel_via_analysis_2_functors >
class Curved_kernel_via_analysis_2_base 
       : public FunctorBase< NewCKvA, BaseCKvA > {
    
public:
    //!\name Global types
    //!@{
    
    //! type of curve kernel
    typedef CurveKernel_2 Curve_kernel_2;

    //! rebinds functor base to a new CKvA
    template <class X>
    struct rebind {

        typedef FunctorBase<X, typename
            BaseCKvA::Curved_kernel_via_analysis_2> Functor_base;      
    };

    //!\name Embedded types to fulfill \c ArrangementTraits_2 concept

    //! type of curve that can be analyzed
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_2;

    //! type of non x-monotone arc on a curve that can be analyzed
    typedef internal::Non_x_monotone_arc_2<NewCKvA>
             Non_x_monotone_arc_2;

    //! the multiplicity type
    typedef unsigned int Multiplicity;
    
    //!@}

public:

    //!\name Tags
    //!@{
    
    //! tag specifies that "to the left of" comparisons are supported
    typedef CGAL::Tag_true Has_left_category;

    //! tag specifies that merge and split functors are supported
    typedef CGAL::Tag_true Has_merge_category; 

    typedef CGAL::Tag_false Has_do_intersect_category;

    typedef Arr_open_side_tag Left_side_category;
    typedef Arr_open_side_tag Bottom_side_category;
    typedef Arr_open_side_tag Top_side_category;
    typedef Arr_open_side_tag Right_side_category;

    //!@}

public:
    //!\name Caching
    
    //!@{
    
    //! type of inverval arcno cache
    typedef internal::Curve_interval_arcno_cache< Curve_kernel_2 > 
        Curve_interval_arcno_cache;

    //!@}
    
public:
    //! \name Constructors
    //!@{

    //! default constructor
    Curved_kernel_via_analysis_2_base() :
        _m_kernel(Curve_kernel_2()) {
    }

    //! construct using specific Curve_kernel_2 instance \c kernel
    Curved_kernel_via_analysis_2_base(const Curve_kernel_2& kernel) :
        _m_kernel(kernel) {
    }
    
    //!@}

    //!\name underlying curve kernel + caching
    //!@{
    
    /*!\brief
     * access to static Curve_interval_arcno_cache
     */
    const Curve_interval_arcno_cache& interval_arcno_cache() const {
        return _m_interval_arcno_cache;
    }
            
    /*!\brief
     * instance of internal Curve_kernel_2 instance
     *
     * \return 
     */
    Curve_kernel_2& kernel() const {
        return _m_kernel;
    }

    //!@}

protected:
    //!\name private members
    //!@{
    
    //! an instance of \c Curve_kernel_2
    mutable Curve_kernel_2 _m_kernel;
    
    //! an instance of \c Curve_interval_arcno_cache
    mutable Curve_interval_arcno_cache _m_interval_arcno_cache;
    
    //!@}

public:
    //!\name Static Member to provide CKvA instance
    //!@{

    /*!\brief
     * a default instance of \c Curved_kernel_via_analysis_2
     *
     * \return static instance of \c Curved_kernel_via_analysis_2
     */
    static NewCKvA& instance() {
        return set_instance(_set_instance());
    }
    
    /*!\brief
     * sets static instance of \c Curved_kernel_via_analysis_2 to \c ckva
     *
     * \param ckva The instance that should be stored
     * \return the stored instance
     */
    static NewCKvA& set_instance(
            const NewCKvA& ckva
    ) {
        static NewCKvA instance;
        static NewCKvA binstance;
        
        if (&ckva == &_reset_instance()) {
            instance = binstance; 
        } else if (&ckva != &_set_instance()) {
            binstance = instance;
            instance = ckva;
        }
        return instance;
        
    }
    
    /*!\brief
     * resets static instance to original one
     */ 
    static void reset_instance() {
        set_instance(_reset_instance());
    }
    
private:
    /*!\brief 
     * sets instance to default for internal purposes
     */
    static NewCKvA& _set_instance() {
        static NewCKvA instance;
        return instance;
        
    }
    
    /*!\brief 
     * sets instance to default for internal purposes
     */
    static NewCKvA& _reset_instance() {
        static NewCKvA instance;
        return instance;
    }
    
    //!@}
};    
    
} // namespace internal

template < class CurveKernel_2 >
class Curved_kernel_via_analysis_2 :
    public internal::Curved_kernel_via_analysis_2_base<
        Curved_kernel_via_analysis_2<CurveKernel_2>,
        void, CurveKernel_2 > 
{
public: 
    //! type of curve kernel
    typedef CurveKernel_2 Curve_kernel_2;

    //! this instance itself
    typedef Curved_kernel_via_analysis_2< Curve_kernel_2 > Self;
    
    //! redefine rebind to terminate unrolling nested templates 
    template < class X >
    struct rebind {

        typedef internal::Curved_kernel_via_analysis_2_functors< X >
             Functor_base;      
    };

public:
    //!\name Embedded types to fulfill \c ArrangementTraits_2 concept
    //!@{

    typedef internal::Point_2< Self > Point_2;
    
    typedef internal::Arc_2< Self > Arc_2;
    
    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    //!@}
protected:
    //!\name Protected base types
    //!@{
    
    //! class collecting basic types
    typedef internal::Curved_kernel_via_analysis_2_base < Self, void,
         CurveKernel_2 > Base_kernel;

    //!@}
public:
    //! \name Constructors
    //!@{

    /*!\brief
     * default constructor
     */
    Curved_kernel_via_analysis_2() :
      Base_kernel(Curve_kernel_2::get_static_instance()) {
    }
    
    /*!\brief
     * construct from \c kernel
     *
     * \param kernel Kernel to use internally
     */
    Curved_kernel_via_analysis_2(const Curve_kernel_2& kernel) :
        Base_kernel(kernel) {
    }
    
    //!@}
}; 

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_IMPL_H
// EOF
