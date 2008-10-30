// Copyright (c) 2007-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#ifndef CGAL_ARRANGEMENT_2l_SURFACE_3_ENVELOPE_TRAITS
#define CGAL_ARRANGEMENT_2l_SURFACE_3_ENVELOPE_TRAITS 1

/*!\file include/CGAL/Arrangement_2l/Surface_3_envelope_traits.h
 * \brief Model for CGAL's EnvelopeTraits_3 concept.
 */

#ifndef CGAL_ENVELOPE_3_USE_EDGE_HANDLES
#define CGAL_ENVELOPE_3_USE_EDGE_HANDLES 0
#endif

#include <CGAL/config.h>

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Envelope_3/Envelope_base.h>

#include <CGAL/Arrangement_2l/Restricted_cad_3_accessor.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

#if CGAL_ENVELOPE_3_USE_EDGE_HANDLES
template < class SurfacePair_3 >
struct Intersection_info : public CGAL::Handle_with_policy<
std::map< SurfacePair_3,  
          typename SurfacePair_3::Restricted_cad_3::Edge_const_handle,
          CGAL::Handle_id_less_than< SurfacePair_3 >
>
> {
    
    typedef SurfacePair_3 Surface_pair_3;
    
    typedef typename Surface_pair_3::Restricted_cad_3::Edge_const_handle 
    Edge_const_handle;

    typedef CGAL::Handle_id_less_than< Surface_pair_3 > Less;
    
    typedef std::map< SurfacePair_3, Edge_const_handle, Less > Rep;
    
    typedef CGAL::Handle_with_policy< Rep > Base;

    typedef Intersection_info< Surface_pair_3 > Self;

    typedef typename Rep::const_iterator const_iterator;

    typedef typename Rep::value_type value_type;

private:
    typedef typename Rep::iterator iterator;
    
public:

    Intersection_info() {
    }

    Intersection_info(const Surface_pair_3& surface, 
                      const Edge_const_handle& eh) {
        CGAL_precondition(this->ptr()->empty());
        this->ptr()->insert(std::make_pair(surface, eh));
    }

    const_iterator begin() const {
        return this->ptr()->begin();
    }

    const_iterator end() const {
        return this->ptr()->end();
    }

    const_iterator find(const Surface_pair_3& pair) const {
        return this->ptr()->find(pair);
    }
    
    Self operator() (Self d1, Self d2) {
        Self tmp = d2;
        tmp.copy_on_write();
        
        for (iterator it1 = d1.ptr()->begin();
             it1 != d1.ptr()->end(); it1++) {
            iterator it2 = tmp.ptr()->find(it1->first);
            if (it2 == tmp.ptr()->end()) {
                tmp.ptr()->insert(it2, *it1);
            }
        }
        return tmp;
    }
    
    bool operator== (Self self) const {
        return this->id() == self.id();
    }
};
#endif

} // namespace CGALi


// TODO implement Apollonius-Mode
// TODO implicitly requires Surface_pair_3::surface_pair_cache()
template < class SurfacePair_3 >
class Surface_3_envelope_traits : public
#if CGAL_ENVELOPE_3_USE_EDGE_HANDLES
CGAL::Arr_curve_data_traits_2<
typename SurfacePair_3::Surface_z_at_xy_isolator_traits::Arrangement_traits_2,
CGAL::CGALi::Intersection_info< SurfacePair_3 >,
CGAL::CGALi::Intersection_info< SurfacePair_3 >
>
#else
SurfacePair_3::Surface_z_at_xy_isolator_traits::Arrangement_traits_2
#endif
{
public:
    //! this instance template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the class itself
    typedef Surface_3_envelope_traits< Surface_pair_3 > Self;

    // types for EnvelopeTraits
    //! type of surfaces
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of Xy_monotone_surface_3
    typedef Surface_3 Xy_monotone_surface_3;
    
    //! type of multiplicity
    typedef unsigned int Multiplicity;

private:
    //! type of Restricted_cad
    typedef typename Surface_pair_3::Restricted_cad_3 Restricted_cad_3;

    //! type of surface traits
    typedef typename Surface_pair_3::Surface_z_at_xy_isolator_traits
    Surface_z_at_xy_isolator_traits;
    
    //! type of basic traits
    typedef typename Surface_z_at_xy_isolator_traits::Arrangement_traits_2
    Arrangement_traits_2;
    
    //! type of Edge_const_handle
    typedef typename Restricted_cad_3::Edge_const_handle Edge_const_handle;

#if CGAL_ENVELOPE_3_USE_EDGE_HANDLES
    typedef CGAL::CGALi::Intersection_info< Surface_pair_3 > Info;

    //! type of base class
    typedef 
    CGAL::Arr_curve_data_traits_2< Arrangement_traits_2, Info, Info > Base;
#else
    //! type of base class
    typedef Arrangement_traits_2 Base;
#endif
    
    //! type of Z_Stack
    typedef typename Restricted_cad_3::Z_stack Z_stack;


public:
    
    //! class of point
    typedef typename Base::Point_2 Point_2;
    
    //! class of x-monotone curve
    typedef typename Base::X_monotone_curve_2 X_monotone_curve_2;
    
public:

    /*!\brief
     * Subdivide the given surface into envelope relevant xy-monotone 
     * parts, and insert them into the output iterator.
     * 
     * The iterator value-type is Xy_monotone_surface_3
     */
    class Make_xy_monotone_3  
    {
    protected:
        const Self *parent; 
    public:
        Make_xy_monotone_3(const Self* p)
            : parent(p)
	{}


        template <class OutputIterator>
        OutputIterator operator()(
                const Surface_3& s,
                bool is_lower,
                OutputIterator oi) 
	{
            //parent->total_timer.start();
            
            // TODO surfaces must be coprime? 
            // TASK Ask RW, MM for triangles, and planes

            // we just apply the sophisticated stuff in sqff_3
            // that also deals with the case of vertical components
            typename 
                Surface_z_at_xy_isolator_traits::Square_free_factorization_3
                sqff_3; // TODO unique instance!
            
            std::list< int > mults;
            sqff_3(s.f(), oi, std::back_inserter(mults)); 
            
            //parent->total_timer.stop();
            
            return oi;
        }
    };

    /*! Get a Make_xy_monotone_3 functor object. */
    Make_xy_monotone_3 make_xy_monotone_3_object()
    {
        return Make_xy_monotone_3(this);
    }

    /*!\brief
     * Insert all 2D curves, which form of the boundary of the 
     * vertical projection of s onto the xy-plane, into the output iterator.
     * The iterator value-type is Curve_2.
     */
    class Construct_projected_boundary_2 {
    protected:
        const Self *parent; 
    public:
        Construct_projected_boundary_2(const Self* p)
	    : parent(p)
	{}

        template <class OutputIterator>
        OutputIterator operator()(
                const Xy_monotone_surface_3& s,
                OutputIterator oi) {
           
#if !NDEBUG
            std::cout << "Construct_projected_boundary ... " << std::flush;
#endif
            //parent->total_timer.start();
	    //parent->pboundary_timer.start();
            
            // create cached instance
            Restricted_cad_3 cad = Restricted_cad_3::cad_cache()(s);
            
            // we run over all segments/isolated points of cad and
            // check whether they belong to boundary of s
            
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     cad.vertices_begin();
                 vit != cad.vertices_end(); vit++) {
                if (vit->is_isolated() && cad.has_silhouette(vit)) {
                    Z_stack z_stack = cad.z_stack(vit);
                    // we are only interested in the lowest boundary!
                    if (!z_stack.is_empty()) {
                        if (z_stack.level_of_surface_in_z_cell(s,0) == 0) {
                            *oi++ = CGAL::make_object(vit->point());
                        }
                    }
                }
            }
            for (typename Restricted_cad_3::Edge_const_iterator eit =
                     cad.edges_begin();
                 eit != cad.edges_end(); eit++) {
                if (cad.has_silhouette(eit)) {
                    //std::cout << "eit->curve(): " 
                    //          << eit->curve() << std::endl;
                    Z_stack z_stack = cad.z_stack(eit);
                    if (!z_stack.is_empty()) {
                        // we are only interested in the lowest boundary
                        // TODO check whether it is cheaper to check z_stacks 
                        // of incident faces first, instead of 
                        // checking z_stacks for eit
                        if (z_stack.level_of_surface_in_z_cell(s,0) == 0) {
                            CGAL::Arr_halfedge_direction dir = 
                                eit->direction();
                            // check both incident faces of eit:
                            // 1) both stacks are empty: 
                            //    return ON_ORIENTED_BOUNDARY
                            // 2) both stacks are non-empty:
                            //    this part of the curve isn't 
                            //    a true boundary of s
                            // 3) one is empty, the other not
                            //    we indicate the non-empty side!
                            int k = cad.z_stack(
                                    eit->face()
                            ).number_of_z_cells();
                            int l = cad.z_stack(
                                    eit->twin()->face()
                            ).number_of_z_cells();
                            CGAL_assertion(k >= 0);
                            CGAL_assertion(l >= 0);
                            //std::cout << "k: " << k << std::endl;
                            //std::cout << "l: " << l << std::endl;
                            CGAL::Oriented_side side = 
                                CGAL::ON_ORIENTED_BOUNDARY;
                            if (k > 0) {
                                if (l > 0) {
                                    continue;
                                } else {
                                    CGAL_assertion(l == 0);
                                    side = CGAL::ON_NEGATIVE_SIDE;
                                }
                            } else {
                                CGAL_assertion(k == 0);
                                if (l > 0) {
                                    side = CGAL::ON_POSITIVE_SIDE;
                                }
                            }
                            if (dir == CGAL::ARR_LEFT_TO_RIGHT) {
                                side = -side;
                            }
                            // TODO deal with vertical-case
                            if (k + l > 0) {
                                *oi++ = CGAL::make_object(
                                        std::make_pair(
                                                X_monotone_curve_2(
                                                        eit->curve()
                                                ),
                                                side
                                        )
                                );
                            }
                        }
                    }
                }
            }
            
            //parent->total_timer.stop();
	    //parent->pboundary_timer.stop();

#if !NDEBUG
            std::cout << "done." << std::endl;
#endif
            return oi;
        }
    };

    /*! Get a Construct_projected_boundary_2 functor object. */
    Construct_projected_boundary_2
    construct_projected_boundary_2_object() const
    {
        return Construct_projected_boundary_2(this);
    }


    /*!\brief
     * Insert all the 2D projections (onto the xy-plane) of the 
     * intersection objects between s1 and s2 into the output iterator.
     * 
     * The iterator value-type is Object. An Object may be:
     * 1. A pair<Curve_2,Intersection_type>, where the intersection 
     * type is an enumeration that can take the values
     * {Transversal, Tangency, Unknown}.
     * 2. A Point_2 instance (in degenerate cases).
     */
    class Construct_projected_intersections_2 {
    protected:
        const Self *parent; 
    public:
        Construct_projected_intersections_2(const Self* p)
	    : parent(p)
        {}

        template <class OutputIterator>
        OutputIterator operator()(
                const Xy_monotone_surface_3& s1,
                const Xy_monotone_surface_3& s2,
                OutputIterator oi) {
#if !NDEBUG
            std::cout << "Construct_projected_intersections ... " 
                      << std::flush;
#endif
            //parent->total_timer.start();
	    //parent->intersection_timer.start();
            
            Surface_pair_3 pair = 
                Surface_pair_3::surface_pair_cache()(std::make_pair(s1, s2));

            Restricted_cad_3 cad = pair.silhouettes_cut();

            // run over all segments/isolated points of cad and
            // check whether they belong to intersection of s1 and s2
            for (typename Restricted_cad_3::Vertex_const_iterator vit =
                     cad.vertices_begin();
                 vit != cad.vertices_end(); vit++) {
                if (vit->is_isolated() && cad.has_cut(vit)) {
                    Z_stack z_stack = cad.z_stack(vit);
                    if (!z_stack.is_empty()) {
                        if (z_stack.level_of_surface_in_z_cell(s1,0) == 0 &&
                            z_stack.level_of_surface_in_z_cell(s2,0) == 0) {
                            *oi++ = CGAL::make_object(vit->point());
                        }
                    }
                }
            }
            
            for (typename Restricted_cad_3::Edge_const_iterator eit =
                     cad.edges_begin();
                 eit != cad.edges_end(); eit++) {
                if (cad.has_cut(eit)) {  
                    Z_stack z_stack = cad.z_stack(eit);
                    if (!z_stack.is_empty()) {
                        if (z_stack.level_of_surface_in_z_cell(s1,0) == 0 &&
                            z_stack.level_of_surface_in_z_cell(s2,0) == 0) {
                            // TODO compute multiplicity of intersection
                            // if mult = 1 -> give 1
                            // if mult = 2 
                            // -> check if unique intersection - give 2
                            // if mult > 2 
                            // -> complex intersections destroy picture
                            Multiplicity multiplicity = 0; 
#if CGAL_ENVELOPE_3_USE_EDGE_HANDLES
                            Info info(pair, eit);
#endif
                            *oi++ = CGAL::make_object(
                                    std::make_pair(
                                            // store eit in info
#if CGAL_ENVELOPE_3_USE_EDGE_HANDLES
                                            X_monotone_curve_2(eit->curve(),
                                                               info),
#else
                                            X_monotone_curve_2(eit->curve()),
#endif
                                            multiplicity
                                    )
                            );
                        }
                    }
                }
            }
            
            //parent->total_timer.stop();
	    //parent->intersection_timer.stop();
            
#if !NDEBUG
            std::cout << "done." << std::endl;
#endif

            return oi;
        }
    };

    /*! Get a Construct_projected_intersections_2 functor object. */
    Construct_projected_intersections_2
    construct_projected_intersections_2_object() const
    {
        return Construct_projected_intersections_2(this);
    }

private:

    //! compares the vertical alignment of \c s1 and \c s2 along \c z_stack
    static
    CGAL::Comparison_result _intern_compare(const Xy_monotone_surface_3& s1,
                                            const Xy_monotone_surface_3& s2,
                                            const Z_stack& z_stack,
                                            bool reverse) {
        CGAL_assertion(!z_stack.is_empty());
        int i1 = z_stack.level_of_surface_in_z_cell(s1,0);
        int i2 = z_stack.level_of_surface_in_z_cell(s2,0);
        
        CGAL_assertion(i1 >= -1 || i2 >= -1);
        CGAL_assertion(i1 < 1);
        CGAL_assertion(i2 < 1);
        if (i1 == i2) {
            CGAL_assertion(i1 == 0);
            CGAL_assertion(i2 == 0);
            return CGAL::EQUAL;
        } else if (i1 == 0) {
            if (reverse) {
                return CGAL::LARGER;
            }
            return CGAL::SMALLER;
        }
        // else 
        CGAL_assertion(i2 == 0);
        if (reverse) {
            return CGAL::SMALLER;
        }
        return CGAL::LARGER;
    }

public:
    
    /*!\brief
     * Check if the surface s1 is closer/equally distanced/farther 
     * from the envelope with respect to s2 at the xy-coordinates of p/c.
     */
    class Compare_z_at_xy_3 {
    protected:
        const Self *parent; 
    public:
	Compare_z_at_xy_3(const Self* p)
	  : parent(p)
	{}

        //! compare over point
        CGAL::Comparison_result operator()(
                const Point_2& p,
                const Xy_monotone_surface_3& s1,
                const Xy_monotone_surface_3& s2) const {
#if !NDEBUG
            std::cout << "Compare over point ... " << std::flush;
#endif
            //parent->total_timer.start();
	    //parent->compare_timer.start();      
            
            Surface_pair_3 pair = 
                Surface_pair_3::surface_pair_cache()(std::make_pair(s1, s2));
            
            bool reverse = (s1 == pair.surface2());
            
            Restricted_cad_3 cad = pair.silhouettes_cut();

            // locate point
            // TODO develop strategy to avoid point location!
            Z_stack z_stack = cad.z_stack_for(p);
            
            // and compare lowest z-cells of s1 and s2
            CGAL::Comparison_result result = 
                _intern_compare(s1, s2, z_stack, reverse);
            
            //parent->total_timer.stop();
	    //parent->compare_timer.stop();
#if !NDEBUG
            std::cout << "done." << std::endl;
#endif

	    return result;
        }
        
        //! compare over curve 
        CGAL::Comparison_result operator() (
                const X_monotone_curve_2& cv,
                const Xy_monotone_surface_3& s1,
                const Xy_monotone_surface_3& s2) const {
            
#if !NDEBUG
            std::cout << "Compare over curve ... " << std::flush;
#endif	    
            //parent->total_timer.start();
	    //parent->compare_on_cv_timer.start();	    
	    
            Surface_pair_3 pair =
                Surface_pair_3::surface_pair_cache()(std::make_pair(s1, s2));
            
            bool reverse = (s1 == pair.surface2());
            
            Restricted_cad_3 cad = pair.silhouettes_cut();
            
            typedef Restricted_cad_3_accessor< Restricted_cad_3 > Accessor;
            Accessor acc(cad);

            // construct point on c
            Point_2 p = acc.point_in_interior(cv);
            
            // locate point 
            // TODO develop strategy to avoid construction/point location!
            Z_stack z_stack = cad.z_stack_for(p);
            
            // and compare lowest z-cells of s1 and s2
            CGAL::Comparison_result result = 
                _intern_compare(s1, s2, z_stack, reverse);

            // the two surfaces are not allowed to intersect ove c
            CGAL_postcondition(result != CGAL::EQUAL);

	    //parent->total_timer.stop();
	    //parent->compare_on_cv_timer.stop();
            
#if !NDEBUG
            std::cout << "done." << std::endl;
#endif
            
	    return result;
        }
        
        //! compare over full surfaces
        CGAL::Comparison_result operator() (
                const Xy_monotone_surface_3& s1,
                const Xy_monotone_surface_3& s2) const {
#if !NDEBUG
            std::cout << "Compare over face ... " << std::flush;
#endif            
            //parent->total_timer.start();
	    //parent->compare_on_face_timer.stare();
            
            Surface_pair_3 pair =
                Surface_pair_3::surface_pair_cache()(std::make_pair(s1, s2));
            
            bool reverse = (s1 == pair.surface2());

            Restricted_cad_3 cad = pair.silhouettes_cut();
            
            // use z_stack of only existing face
            Z_stack z_stack = cad.z_stack(cad.faces_begin());
            
            // and compare lowest z-cells of s1 and s2
            CGAL::Comparison_result result = 
                _intern_compare(s1, s2, z_stack, reverse);
            
            // the two surfaces are not allowed to intersect at all
            CGAL_postcondition(result != CGAL::EQUAL);
            
            //parent->total_timer.stop();
	    //parent->compare_on_face_timer.stop();
            
#if !NDEBUG
            std::cout << "done." << std::endl;
#endif
            return result;
        }
    };
    
    /*! Get a Compare_z_at_xy_3 functor object. */
    Compare_z_at_xy_3 
    compare_z_at_xy_3_object() const
    {
        return Compare_z_at_xy_3(this);
    }

    /*!\brief 
     * Check if the surface s1 is closer/equally distanced/farther 
     * from the envelope with
     * respect to s2 immediately below the curve c. 
     */
    class Compare_z_at_xy_below_3 {
    protected:
        const Self *parent; 
    public:
        Compare_z_at_xy_below_3(const Self* p)
           : parent(p)
        {}
        CGAL::Comparison_result operator() (
                const X_monotone_curve_2& c,
                const Xy_monotone_surface_3& s1,
                const Xy_monotone_surface_3& s2) const {
#if !NDEBUG
            std::cout << "Compare below curve ... " << std::flush;
#endif
            //parent->total_timer.start();
	    //parent->compare_side_timer.start();
	    
            Surface_pair_3 pair =
                Surface_pair_3::surface_pair_cache()(std::make_pair(s1, s2));
            
            bool reverse = (s1 == pair.surface2());

            Restricted_cad_3 cad = pair.silhouettes_cut();
            
            // not invalid, i.e., formed by an intersection and not a boundary
#if CGAL_ENVELOPE_3_USE_EDGE_HANDLES
            typename Info::const_iterator dit = c.data().find(pair);
            CGAL_assertion(dit != c.data().end());
            Edge_const_handle eh = dit->second;

            CGAL_assertion(cad.has_cut(eh,s1,s2));
            
            // locate halfedge of c in cad using stored information
            // as we want to compare "below" we chose the face incident the
            // the halfegde of the pair that is directed from RIGHT_TO_LEFT
            // and compare the lowest z-cells over this face
            // Remark: This also works when c is vertical!
            typename Restricted_cad_3::Face_const_handle fh =
                (eh->direction() == CGAL::ARR_RIGHT_TO_LEFT ? 
                 eh->face() : eh->twin()->face());
            Z_stack z_stack = cad.z_stack(fh);
#else
            typedef Restricted_cad_3_accessor< Restricted_cad_3 > Accessor;
            Accessor acc(cad);
            
            // construct point on c
            Point_2 p = acc.point_in_interior(c);

            // TODO develop strategy to avoid point location!
            CGAL::Object obj = cad.locate(p);
            typename Accessor::Halfedge_const_handle heh;
            CGAL_assertion_code(bool check = )
                CGAL::assign(heh, obj);
            CGAL_assertion(check);
            
            // locate halfedge of c in cad using stored information
            // as we want to compare "below" we chose the face incident the
            // the halfegde of the pair that is directed from RIGHT_TO_LEFT
            // and compare the lowest z-cells over this face
            // Remark: This also works when c is vertical!
            typename Restricted_cad_3::Face_const_handle fh =
                (heh->direction() == CGAL::ARR_RIGHT_TO_LEFT ? 
                 heh->face() : heh->twin()->face());
            Z_stack z_stack = cad.z_stack(fh);
#endif
            
            // and compare lowest z-cells of s1 and s2
            CGAL::Comparison_result result = 
                _intern_compare(s1, s2, z_stack, reverse);
            
            // the two surfaces are not allowed to intersect in a 
            // two-dimensional patch
            CGAL_postcondition(result != CGAL::EQUAL);

	    //parent->total_timer.stop();
	    //parent->compare_side_timer.stop();
            
#if !NDEBUG
            std::cout << "done." << std::endl;
#endif
	    return result;
        }
        
    };

    /*! Get a Compare_z_at_xy_below_3 functor object. */
    Compare_z_at_xy_below_3
    compare_z_at_xy_below_3_object() const
    {
        return Compare_z_at_xy_below_3(this);
    }
    
    /*!\brief 
     * Check if the surface s1 is closer/equally distanced/farther 
     * from the envelope with
     * respect to s2 immediately above the curve c. 
     */
    class Compare_z_at_xy_above_3 {
    protected:
        const Self *parent; 
    public:
        Compare_z_at_xy_above_3(const Self* p)
	  : parent(p)
	{}

        CGAL::Comparison_result operator()(
                const X_monotone_curve_2& c,
                const Xy_monotone_surface_3& s1,
                const Xy_monotone_surface_3& s2) {
#if !NDEBUG
            std::cout << "Compare above curve ... " << std::flush;
#endif
            //parent->total_timer.start();
	    //parent->compare_side_timer.start();

            Surface_pair_3 pair =
                Surface_pair_3::surface_pair_cache()(std::make_pair(s1, s2));
            
            bool reverse = (s1 == pair.surface2());

            Restricted_cad_3 cad = pair.silhouettes_cut();
            
#if CGAL_ENVELOPE_3_USE_EDGE_HANDLES
            // not invalid, i.e., formed by an intersection and not a boundary
            typename Info::const_iterator dit = c.data().find(pair);
            CGAL_assertion(dit != c.data().end());
            Edge_const_handle eh = dit->second;
            
            CGAL_assertion(cad.has_cut(eh,s1,s2));

            // locate halfedge of c in cad using stored information
            // as we want to compare "above" we chose the face incident the
            // the halfegde of the pair that is directed from LEFT_TO_RIGHT
            // and compare the lowest z-cells over this face
            // Remark: This also works when c is vertical!
            typename Restricted_cad_3::Face_const_handle fh =
                (eh->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 eh->face() : eh->twin()->face());
            Z_stack z_stack = cad.z_stack(fh);
#else
            typedef Restricted_cad_3_accessor< Restricted_cad_3 > Accessor;
            Accessor acc(cad);
            
            // construct point on c
            Point_2 p = acc.point_in_interior(c);

            // TODO develop strategy to avoid point location!
            CGAL::Object obj = cad.locate(p);
            typename Accessor::Halfedge_const_handle heh;
            CGAL_assertion_code(bool check = )
                CGAL::assign(heh, obj);
            CGAL_assertion(check);
            
            // locate halfedge of c in cad using stored information
            // as we want to compare "below" we chose the face incident the
            // the halfegde of the pair that is directed from RIGHT_TO_LEFT
            // and compare the lowest z-cells over this face
            // Remark: This also works when c is vertical!
            typename Restricted_cad_3::Face_const_handle fh =
                (heh->direction() == CGAL::ARR_LEFT_TO_RIGHT ? 
                 heh->face() : heh->twin()->face());
            Z_stack z_stack = cad.z_stack(fh);
#endif
            // and compare lowest z-cells of s1 and s2
            CGAL::Comparison_result result = 
                _intern_compare(s1, s2, z_stack, reverse);
            
            // the two surfaces are not allowed to intersect in a 
            // two-dimensional patch
            CGAL_postcondition(result != CGAL::EQUAL);
            
	    //parent->total_timer.stop();
	    //parent->compare_side_timer.stop();
            
#if !NDEBUG
            std::cout << "done." << std::endl;
#endif
            return result;
        }
    };

    /*! Get a Compare_z_at_xy_above_3 functor object. */
    Compare_z_at_xy_above_3
    compare_z_at_xy_above_3_object() const
    {
        return Compare_z_at_xy_above_3(this);
    }

}; // Surface_3_envelope_traits

CGAL_END_NAMESPACE

#endif // CGAL_ARRANGEMENT_2l_SURFACE_3_ENVELOPE_TRAITS
// EOF
