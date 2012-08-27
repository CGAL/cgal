// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ACK_SHEAR_TRANSFORMATION
#define CGAL_ACK_SHEAR_TRANSFORMATION 1

#include <CGAL/basic.h>

#include <vector>

#include <CGAL/Algebraic_kernel_d/macros.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_coefficient_kernel_at_alpha.h>
#include <CGAL/Algebraic_kernel_d/shear.h>

namespace CGAL {

/*!
 * The class is a functor, getting an algebraic curve and some
 * shear factor, and returning the sheared curve.
 */
template<typename AlgebraicKernelWithAnalysis_2> class Shear_transformation {

public:
      
    typedef AlgebraicKernelWithAnalysis_2  Algebraic_kernel_with_analysis_2;

    typedef typename Algebraic_kernel_with_analysis_2::Curve_analysis_2 
        Curve_analysis_2;

    typedef typename  AlgebraicKernelWithAnalysis_2::Polynomial_traits_2 
        Polynomial_traits_2;

    CGAL_ACK_SNAP_ALGEBRAIC_CURVE_KERNEL_2_TYPEDEFS(Curve_analysis_2);

    typedef std::pair<Bound,Bound> Point;

    typedef std::vector< Algebraic_real_1 > Root_container;

    typedef typename Root_container::iterator Root_iterator;
    
    typedef typename Curve_analysis_2::Event_line_iterator 
    Status_line_1_iterator;

private:

    struct Y_structure_element;

    typedef std::list<Y_structure_element> Y_structure;

    // TODO replace by something that we already have
    enum Coor_type { MINUS_INFTY,FINITE,PLUS_INFTY};

public:

    Shear_transformation(Algebraic_kernel_with_analysis_2* kernel)
        : _m_kernel(kernel),
          x_extreme_index_counter(0),
          disc_roots_computed(false),
          sh_disc_roots_computed(false)
    {}

    template<typename InputIterator>      
    void report_sheared_disc_roots(InputIterator begin,
                                   InputIterator end) {
        std::copy(begin,end,std::back_inserter(sh_disc_roots));
        sh_disc_roots_computed=true;
    }

    Curve_analysis_2 operator() (const Curve_analysis_2& C, Integer s,
                      bool use_primitive_curve=true) {
        Curve_analysis_2 D;
        this->operator() (C,s,D,use_primitive_curve);
        return D;
    }

    void operator() (const Curve_analysis_2& C, Integer s, Curve_analysis_2& D, 
                     bool use_primitive_curve=true) {
        this->C=C;
        this->s=s;
        this->use_primitive_curve = use_primitive_curve;
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Curve_analysis_2: " 
                             << C.polynomial_2() << std::endl;
        CGAL_ACK_DEBUG_PRINT << "num events: " 
                             << C.number_of_status_lines_with_event() 
                             << std::endl;
        CGAL_ACK_DEBUG_PRINT << "s: " << s << std::endl;
#endif
*/
        x_structure.clear();
        /*
          sh_disc_roots.clear();
          x_structure_info.clear();

          ev_res_roots_mults.clear();
          sh_ev_indices.clear();
          stripe_values.clear();
          pre_vert_lines.clear();
          sh_intermediate_lines.clear();
        */
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Compute the polynomials.." << std::flush;
#endif

        if(this->use_primitive_curve) {
            pol = C.primitive_polynomial_2();
        } else {
            pol=C.polynomial_2();
        }
        sh_pol=CGAL::internal::shear(pol,Coefficient(s));
        if(CGAL::degree(typename Polynomial_traits_2
                  ::Univariate_content_up_to_constant_factor()( sh_pol ))>0) {
            throw CGAL::internal::Non_generic_position_exception();
        }
        if(! D.has_defining_polynomial()) {
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "set f.." << std::flush;
#endif
            D.set_f(sh_pol);
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "done.." << std::flush;
#endif
        }
        
        der_sh_pol = typename Polynomial_traits_2::Differentiate() (sh_pol,1);
        sh_der_sh_pol = CGAL::internal::shear(der_sh_pol,Coefficient(-s));
      
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif

        Solve_1 solve_1;
        
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Store the discriminant roots.." 
                             << std::flush;
#endif

        Root_container disc_roots;
        for(Status_line_1_iterator it=C.event_begin();
            it!=C.event_end();
            it++) {
            disc_roots.push_back(it->x());
        }
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif

        if(! sh_disc_roots_computed) {
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Compute the sheared discriminant.." 
                                 << std::flush;
#endif

            if(typename Polynomial_traits_2::Degree() (sh_pol) > 0) {
            

                Polynomial_1 sh_disc 
                    = CGAL::resultant(sh_pol,der_sh_pol);
                
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "root isolation.." << std::flush;
#endif
                solve_1(sh_disc,std::back_inserter(sh_disc_roots),false);
                
            }
             
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
            sh_disc_roots_computed=true;
        }

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Compute the event resultant.." << std::flush;
#endif
        Root_container ev_res_roots;
        if(typename Polynomial_traits_2::Degree() (sh_pol) > 0) {
            Polynomial_1 ev_res = CGAL::resultant(pol,sh_der_sh_pol);
            
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "root isolation.." << std::flush;
#endif
            solve_1(ev_res,std::back_inserter(ev_res_roots),false);
            
        }
        
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done, " << ev_res_roots.size() 
                             << " roots found" << std::endl;
#endif

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Merge both root sets..." << std::flush;
#endif
        typename CGAL::Real_embeddable_traits<Algebraic_real_1>::Compare 
            x_compare;

        CGAL::internal::set_union_with_source
            (disc_roots.begin(),
             disc_roots.end(),
             ev_res_roots.begin(),
             ev_res_roots.end(),
             std::back_inserter(x_structure),
             std::back_inserter(x_structure_info),
             x_compare);
 
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
        CGAL_ACK_DEBUG_PRINT << "Take the stripe values..." << std::flush;
#endif

        CGAL::internal::find_intermediate_values
	  (kernel(),
	   sh_disc_roots.begin(),
	   sh_disc_roots.end(),
	   std::back_inserter(stripe_values));
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
        CGAL_ACK_DEBUG_PRINT << "Search sheared event points..." << std::flush;
#endif
        find_sheared_event_points();

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
        CGAL_ACK_DEBUG_PRINT << "Find start- and endpoints for sweep.." 
                             << std::flush;
#endif
        find_far_points(D);

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
        CGAL_ACK_DEBUG_PRINT << "Start sweep..." << std::flush;
#endif
        
        sweep();
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
        CGAL_ACK_DEBUG_PRINT << "vert lines info:" << std::endl;

        for(int i=0;i<static_cast<int>(pre_vert_lines.size());i++) {
            
            CGAL_ACK_DEBUG_PRINT << "At: " 
                                 << CGAL::to_double(sh_disc_roots[i]) << ", "
                                 << pre_vert_lines[i].number_of_non_event_roots
                                 << " non-event-roots, and " 
                                 << pre_vert_lines[i].event_points.size()
                                 << std::endl;  
        }
        CGAL_ACK_DEBUG_PRINT << "Vert_lines.." << std::flush;
#endif

        CGAL_assertion(sh_disc_roots.size()==pre_vert_lines.size());
        std::vector<Status_line_1> sh_ev_lines;
        for(int i=0;i<static_cast<int>(pre_vert_lines.size());i++) {
            sh_ev_lines.push_back(create_event_line(D,i));
        }
        D.set_event_lines(sh_ev_lines.begin(),sh_ev_lines.end(),
                          sh_intermediate_lines.begin(),
                          sh_intermediate_lines.end());
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif

    }

private:
    
    // X-coordinate of the shear of p
    Bound x_sheared(Point p,Integer sh) {
        return p.first-sh*p.second;
    }
    Bound x_sheared(Bound x,Bound y,Integer sh) {
        return x-sh*y;
    }

    int compute_stripe(Status_line_1& ev, int index) {
        int left_index = -1, 
            right_index = static_cast<int>(stripe_values.size()-1);
        Algebraic_real_1 xv = ev.x();
        Bound lx = xv.low(), rx=xv.high(),
            x_iv_size = rx-lx;
        Bound ly = ev.lower_bound(index),
            ry = ev.upper_bound(index);;
        while(left_index < right_index) {
            if(x_iv_size > ry-ly) {
                xv.refine();
                lx = xv.low();
                rx=xv.high();
                x_iv_size=rx-lx;
                continue;
            }
            ev.refine(index);
            ly = ev.lower_bound(index);
            ry = ev.upper_bound(index);
            Bound right(0), left(0);
            left  = (s>0) ? x_sheared(lx,ry,s) : x_sheared(lx,ly,s);
            right = (s>0) ? x_sheared(rx,ly,s) : x_sheared(rx,ry,s);
            CGAL_assertion(left<right);
            while(left_index<right_index && stripe_values[left_index+1]<left) {
                ++left_index;
            }
            while(left_index<right_index && right<stripe_values[right_index]) {
                --right_index;
            }
        }
        CGAL_assertion(left_index==right_index);
        return left_index;

    }

    void find_far_points(Curve_analysis_2& D) {
        int n = static_cast<int>(stripe_values.size());
        Bound upper_bound,lower_bound;
        Bound left_bound = stripe_values[0],
            right_bound=stripe_values[n-1];
        lower_bound = upper_bound = Bound(0);
        for(int i=0;i<n;i++) {
            Algebraic_real_1 curr_bound(stripe_values[i]);
            Bitstream_traits traits(Bitstream_coefficient_kernel
                                    (kernel(),curr_bound));
            CGAL::internal::Square_free_descartes_tag tag;
            Bitstream_descartes descartes(tag,sh_pol,traits);
            int m = descartes.number_of_real_roots();
            if(m>0) {
                if(descartes.left_bound(0)<lower_bound) {
                    lower_bound = descartes.left_bound(0);
                }
                if(descartes.right_bound(m-1) 
                   > upper_bound) {
                    upper_bound = descartes.right_bound(m-1);
                }
            }
            // Create intermediate line for later use
            Algebraic_real_1 xval(curr_bound);
            Status_line_1 inter_line(xval,i,D,m);
            inter_line.set_isolator(descartes);
            sh_intermediate_lines.push_back(inter_line);
        }
        far_left=(s<0) ? x_sheared(left_bound,upper_bound,-s)
            : x_sheared(left_bound,lower_bound,-s)-1;
        far_right=(s<0) ? x_sheared(right_bound,lower_bound,-s)
            : x_sheared(right_bound,upper_bound,-s)+1;
        if(C.number_of_status_lines_with_event()>0) {
            if(far_left>C.status_line_at_event(0).x().low()) {
                far_left = C.status_line_at_event(0).x().low();
            }   
            if(far_right<C.status_line_at_event
               (C.number_of_status_lines_with_event()-1).x().high()) {
                far_right = C.status_line_at_event
                    (C.number_of_status_lines_with_event()-1).x().high();
            }
        }
        // just to be sure...
        far_left = far_left-1;
        far_right = far_right+1;
        y_in_box = (upper_bound + lower_bound)/2;
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "far_left: " << far_left << "\nfar_right: " 
                             << far_right << "\ny_in_box: " << y_in_box 
                             << std::endl;
#endif
*/
    }

    template<class OutputIterator>
    void find_sheared_event_points_at_x(const Status_line_1& ev,
                                        Status_line_1& left,
                                        Status_line_1& right,
                                        OutputIterator out) {

        int ev_id = 0, left_id=0, right_id=0;
        int ev_n = ev.number_of_events(),
            left_n = left.number_of_events(),
            right_n = right.number_of_events();
        (void)left_n;
        (void)right_n;
        // Simple if a vertical line exists
        if(ev.covers_line() && ! use_primitive_curve) {
            for(int i=0;i<ev_n;i++) {
                out++=i;
            }
            return;
        }

        // left side
        
        CGAL_assertion(left.x().is_rational());
        Bound left_x = left.x().rational();
        
	Polynomial_1 left_pol = kernel()->evaluate_utcf_2_object()
	  (typename Polynomial_traits_2::Swap() (pol, 0, 1),
	   left_x);
        Polynomial_1 left_sh_der_sh_pol = kernel()->evaluate_utcf_2_object()
	  (typename Polynomial_traits_2::Swap() (sh_der_sh_pol, 0, 1),
	   left_x);

        // right side

      
        CGAL_assertion(right.x().is_rational());
        Bound right_x = right.x().rational();

	Polynomial_1 right_pol = kernel()->evaluate_utcf_2_object()
	  (typename Polynomial_traits_2::Swap() (pol, 0, 1),
	   right_x);
        Polynomial_1 right_sh_der_sh_pol = kernel()->evaluate_utcf_2_object()
	  (typename Polynomial_traits_2::Swap() (sh_der_sh_pol, 0, 1),
	   right_x);

        int asym_left_minus,asym_left_plus,asym_right_minus,asym_right_plus;
        
        typedef typename Status_line_1::Arc_pair Arc_pair;
        
        Arc_pair apair1 = ev.number_of_branches_approaching_minus_infinity();
        Arc_pair apair2 = ev.number_of_branches_approaching_plus_infinity();

        asym_left_minus = apair1.first;
        asym_right_minus = apair1.second;

        asym_left_plus = apair2.first;
        asym_right_plus = apair2.second;

        left_id += asym_left_minus;
        right_id += asym_right_minus;
        while(ev_id != ev_n) {
        
            typename Status_line_1::Arc_pair arc_pair = 
                ev.number_of_incident_branches(ev_id);

            int arcs_left = arc_pair.first;
            int arcs_right = arc_pair.second;
            if(arcs_left+arcs_right!=2) {
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "Sheared event point found at " 
                                     << ev.x().to_double() << ", index " 
                                     << ev_id << std::endl;
#endif
*/
                out++=ev_id;
                left_id+=arcs_left;
                right_id+=arcs_right;
                ev_id++;
            }
            else {
                if(arcs_left==1 && arcs_right==1) {

                    Algebraic_real_1 left_y(left_pol,
                                          left.lower_bound(left_id),
                                          left.upper_bound(left_id));

                    CGAL::Sign left_sign 
		      = kernel()->sign_at_1_object()
		          (left_sh_der_sh_pol,left_y,true);

                    Algebraic_real_1 right_y(right_pol,
                                           right.lower_bound(right_id),
                                           right.upper_bound(right_id));

                    CGAL::Sign right_sign 
		      = kernel()->sign_at_1_object()
		          (right_sh_der_sh_pol,right_y,true);

                    if(left_sign!=right_sign) {
/*
#if CGAL_ACK_DEBUG_FLAG
                        CGAL_ACK_DEBUG_PRINT << "Sheared ev point found at " 
                                             << ev.x().to_double() 
                                             << ", index " << ev_id 
                                             << std::endl;
#endif
*/
                        out++=ev_id;
                    }
                    ev_id++;
                    left_id++;
                    right_id++;
                }
                else if(arcs_left==2 && arcs_right==0) {
                    Algebraic_real_1 left_y_1(left_pol,
                                          left.lower_bound(left_id),
                                          left.upper_bound(left_id));

                    CGAL::Sign left_sign_1 
		      = kernel()->sign_at_1_object()
		          (left_sh_der_sh_pol,left_y_1,true);

                    left_id++;
                    Algebraic_real_1 left_y_2(left_pol,
                                          left.lower_bound(left_id),
                                          left.upper_bound(left_id));
                    CGAL::Sign left_sign_2 
		      = kernel()->sign_at_1_object()
		          (left_sh_der_sh_pol,left_y_2,true);
                    if(left_sign_1!=left_sign_2) {
/*
#if CGAL_ACK_DEBUG_FLAG
                        CGAL_ACK_DEBUG_PRINT << "Sheared ev point found at " 
                                             << ev.x().to_double() 
                                             << ", index " << ev_id 
                                             << std::endl;
#endif
*/
                        out++=ev_id;
                    }
                    ev_id++;
                    left_id++;
                }
                else if(arcs_left==0 && arcs_right==2) {
                    Algebraic_real_1 right_y_1(right_pol,
                                               right.lower_bound(right_id),
                                               right.upper_bound(right_id));

                    CGAL::Sign right_sign_1 
		      = kernel()->sign_at_1_object()
		          (right_sh_der_sh_pol,right_y_1,true);
                    right_id++;
                    Algebraic_real_1 right_y_2(right_pol,
                                               right.lower_bound(right_id),
                                               right.upper_bound(right_id));
                    CGAL::Sign right_sign_2 
		      = kernel()->sign_at_1_object()
		          (right_sh_der_sh_pol,right_y_2,true);
                    if(right_sign_1!=right_sign_2) {
/*
#if CGAL_ACK_DEBUG_FLAG
                        CGAL_ACK_DEBUG_PRINT << "Sheared ev point found at " 
                                             << ev.x().to_double()  
                                             << ", index " << ev_id 
                                             << std::endl;
#endif
*/
                        out++=ev_id;
                    }
                    ev_id++;
                    right_id++;
                }
            }
        }
        left_id += asym_left_plus;
        right_id += asym_right_plus;
        CGAL_assertion(ev_id==ev_n);
        CGAL_assertion(left_id==left_n);
        CGAL_assertion(right_id==right_n);
    }

    void find_sheared_event_points() {

        sh_ev_indices.resize(x_structure.size());
        std::vector<Bound> intermediate_values;
        find_intermediate_values(kernel(),
				 x_structure.begin(),
                                 x_structure.end(),
                                 std::back_inserter(intermediate_values));
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "interline.." << std::flush;
#endif
*/
        std::vector<Status_line_1> intermediate_lines(intermediate_values.size());

        int i=0;
        for(typename std::vector<Bound>::iterator it 
                = intermediate_values.begin();
            it!=intermediate_values.end();it++) {
            intermediate_lines[i]=C.status_line_at_exact_x(*it);
            i++;
        }
        int event_count=0;
/*
#if CGAL_ACK_DEBUG_FLAG        
        CGAL_ACK_DEBUG_PRINT << "at some x.." << std::flush;
#endif
*/
        for(int i=0;i<static_cast<int>(x_structure.size());i++) {
            CGAL::internal::Three_valued info = x_structure_info[i];
/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT <<  i << "th of " << x_structure.size() 
                                 << std::endl;
            CGAL_ACK_DEBUG_PRINT << "To_double approx" << std::endl;
            CGAL_ACK_DEBUG_PRINT << "x_struct: " 
                                 << CGAL::to_double(x_structure[i]) 
                                 << std::endl;
            CGAL_ACK_DEBUG_PRINT << "Info: " << info << std::endl;
#endif
*/
            if(info==CGAL::internal::ROOT_OF_SECOND_SET || 
               info==CGAL::internal::ROOT_OF_BOTH_SETS) {
                const Status_line_1& event_line_at_x = 
                    (info==CGAL::internal::ROOT_OF_BOTH_SETS)
                    ? C.status_line_at_event(event_count) : C.status_line_at_exact_x(x_structure[i]);
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "now really at x.." << std::flush;
#endif
*/
                find_sheared_event_points_at_x(event_line_at_x,
                                               intermediate_lines[i],
                                               intermediate_lines[i+1],
                                               std::back_inserter
                                               (sh_ev_indices[i]));

/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
*/
            }
            if(info==CGAL::internal::ROOT_OF_FIRST_SET || 
               info==CGAL::internal::ROOT_OF_BOTH_SETS) {
                event_count++;
            }
        }

        CGAL_assertion(event_count==C.number_of_status_lines_with_event());
    }
    
    struct Sh_ev_point_info {

        Sh_ev_point_info(Status_line_1 ev,int index) 
            : ev(ev),index(index),
              incident_left(0),
              incident_right(0)
        {}

        Status_line_1 ev;
        int index;
        int incident_left;
        int incident_right;
    };

    struct Sh_ev_line_info {
        int asym_left_plus,asym_left_minus,asym_right_plus,asym_right_minus;
        int number_of_non_event_roots;
        std::vector<Sh_ev_point_info> event_points;
      
        Bound lower_bound(int i) {
            Sh_ev_point_info p= event_points[i];
            return p.ev.lower_bound(p.index);
        }
        Bound upper_bound(int i) {
            Sh_ev_point_info p= event_points[i];
            return p.ev.upper_bound(p.index);
        }
        void refine(int i) {
            Sh_ev_point_info p= event_points[i];
            p.ev.refine(p.index);
        }

        int num_arcs_left() {
            int sum=0;
            sum+=asym_left_plus+asym_left_minus;
            sum+=number_of_non_event_roots;
            for(int k=0;k<static_cast<int>(event_points.size());k++) {
                sum+=event_points[k].incident_left;
            }
            return sum;
        }

        int num_arcs_right() {
            int sum=0;
            sum+=asym_right_plus+asym_right_minus;
            sum+=number_of_non_event_roots ;
            for(int k=0;k<static_cast<int>(event_points.size());k++) {
                sum+=event_points[k].incident_right;
            }
            return sum;
        }

    };

    struct Y_structure_element {
        bool one_event_known;
        Coor_type x_type;
        int x_index;
        Coor_type y_type;
        int y_index;
    };

    Y_structure_element create_unbounded_element(Status_line_1& ev, int i) {
        int n = static_cast<int>(stripe_values.size());
        int stripe = compute_stripe(ev,i);
        Y_structure_element y_el;
        y_el.one_event_known=true;
        if(stripe==-1) {
            y_el.x_type=MINUS_INFTY;
        } else if(stripe==n-1) {
            y_el.x_type=PLUS_INFTY;
        } else {
            y_el.x_type=FINITE;
            y_el.x_index=stripe;
            while((ev.upper_bound(i)>y_in_box)  &&
                  (ev.lower_bound(i)<y_in_box)) {
                ev.refine(i);
            }
            if(ev.upper_bound(i)<y_in_box) {
                y_el.y_type=MINUS_INFTY;
            }
            else {
                y_el.y_type=PLUS_INFTY;
            }
        }
        return y_el;
    }

    Y_structure_element create_event(Status_line_1& ev, int i) {
        int n = static_cast<int>(stripe_values.size());
        (void)n;
        int stripe = compute_stripe(ev,i);
        CGAL_assertion(stripe>=0 && stripe<=n);
        Y_structure_element y_el;
        y_el.one_event_known=true;
        y_el.x_type=FINITE;
        y_el.x_index=stripe;
        Sh_ev_point_info ev_info(ev,i);
        pre_vert_lines[stripe].event_points.push_back(ev_info);
        y_el.y_type=FINITE;
        y_el.y_index=static_cast<int>
            (pre_vert_lines[stripe].event_points.size()-1);
        return y_el;
    }

    void start_sweep() {
        y_structure.clear();
        for(int i=0;i<static_cast<int>(sh_disc_roots.size());i++) {
            Sh_ev_line_info info;
            info.number_of_non_event_roots=0;
            info.asym_left_plus=info.asym_left_minus=
                info.asym_right_plus=info.asym_right_minus=0;
            pre_vert_lines.push_back(info);
        }
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "X-coordinate (far left) " 
                             << CGAL::to_double(far_left) << std::endl;
#endif
*/
        Status_line_1 far_left_line 
            = C.status_line_at_exact_x(Algebraic_real_1(far_left));
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "No. arcs " 
                             << far_left_line.number_of_events() << std::endl;
#endif
*/
        for(int i=0;i<far_left_line.number_of_events();i++) {
            y_structure.push_back(create_unbounded_element(far_left_line,i));
        }
    }

    
    void end_sweep() {
        Status_line_1 far_right_line 
            = C.status_line_at_exact_x(Algebraic_real_1(far_right));
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "X-coordinate (far right) " 
                             << CGAL::to_double(far_right) << std::endl;
#endif
*/
        CGAL_assertion(far_right_line.number_of_events()
                       ==static_cast<int>(y_structure.size()));
      
        typename Y_structure::iterator y_it=y_structure.begin();
        for(int i=0;i<far_right_line.number_of_events();i++) {
            Y_structure_element y_el 
                = create_unbounded_element(far_right_line,i);
            handle_edge(*y_it,y_el);
            y_it++;
        }
    }

    void sweep_at_x_coordinate(int index) {
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "X-coordinate " 
                             << CGAL::to_double(x_structure[index]) 
                             << std::endl;
#endif
*/
        std::vector<int>::iterator sh_ev_it = sh_ev_indices[index].begin();
        Status_line_1 ev=C.status_line_at_exact_x(x_structure[index]);
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "EV: " << std::endl << ev << std::endl;
#endif
*/
        int ev_id=0, ev_n=ev.number_of_events();
        typename Y_structure::iterator y_it=y_structure.begin();
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Y-structure has " << y_structure.size() 
                             << " elements" << std::endl;
#endif
*/
        // needed for vertical components
        std::vector<Y_structure_element> events_at_x; 
        bool vert=ev.covers_line() && ! this->use_primitive_curve;
        Y_structure_element below,above;
        Y_structure_element minus_x_inf,plus_x_inf;
        minus_x_inf.one_event_known = plus_x_inf.one_event_known=true;
        minus_x_inf.x_type=MINUS_INFTY;
        plus_x_inf.x_type = PLUS_INFTY;
        below = (s>0) ? plus_x_inf : minus_x_inf;
        above = (s>0) ? minus_x_inf: plus_x_inf;
        events_at_x.push_back(below);
        int minus_left,minus_right,plus_left,plus_right;
        
        typedef typename Status_line_1::Arc_pair Arc_pair;
        
        Arc_pair apair1 = ev.number_of_branches_approaching_minus_infinity();
        Arc_pair apair2 = ev.number_of_branches_approaching_plus_infinity();

        minus_left = apair1.first;
        minus_right = apair1.second;

        plus_left = apair2.first;
        plus_right = apair2.second;

        for(int i=0;i<minus_left;i++) {
            handle_edge(*y_it,below);
            y_it=y_structure.erase(y_it);
        }
        for(int i=0;i<minus_right;i++) {
            y_it=y_structure.insert(y_it,below);
            y_it++;
        }
        while(ev_id<ev_n) {
            typename Status_line_1::Arc_pair arc_pair 
                = ev.number_of_incident_branches(ev_id);
            int left_arcs=arc_pair.first;
            int right_arcs=arc_pair.second;
            if(sh_ev_it!=sh_ev_indices[index].end() && ev_id==*sh_ev_it) {
                Y_structure_element y_ev = create_event(ev,ev_id);
                events_at_x.push_back(y_ev);
                //y_struct_info(y_ev);
                for(int i=0;i<left_arcs;i++) {
                    handle_edge(*y_it,y_ev);
                    y_it=y_structure.erase(y_it);
                }
                for(int i=0;i<right_arcs;i++) {
                    y_it=y_structure.insert(y_it,y_ev);
                    y_it++;
                }
                sh_ev_it++;
            } else {
                CGAL_assertion(left_arcs+right_arcs==2);
                if(left_arcs==1 && right_arcs==1) {
                    y_it++;
                }
                else if(left_arcs==2 && right_arcs==0) {
                    Y_structure_element y1,y2;
                    y1=*y_it;
                    y_it=y_structure.erase(y_it);
                    y2=*y_it;
                    y_it=y_structure.erase(y_it);
                    handle_edge(y1,y2);
                }
                else {
                    Y_structure_element new_y;
                    new_y.one_event_known=false;
                    new_y.x_index=x_extreme_index_counter;
            
                    // Prevent compiler warnings:
                    new_y.y_index=-1;
                    new_y.y_type=FINITE;
                    new_y.x_type=FINITE;
            
                    y_it=y_structure.insert(y_it,new_y);
                    y_it++;
                    y_it=y_structure.insert(y_it,new_y);
                    y_it++;
                    x_extreme_index_counter++;
                }     
            }
            ev_id++;
        }
        for(int i=0;i<plus_left;i++) {
            handle_edge(*y_it,above);
            y_it=y_structure.erase(y_it);
        }
        for(int i=0;i<plus_right;i++) {
            y_it=y_structure.insert(y_it,above);
            y_it++;
        }
        CGAL_assertion(ev_id==ev_n);
        CGAL_assertion(y_it==y_structure.end());
        CGAL_assertion(sh_ev_it==sh_ev_indices[index].end());
        if(vert) {
            events_at_x.push_back(above);
            // edges corresponding to vertical segments...
            for(int i=1;i<static_cast<int>(events_at_x.size());i++) {
                handle_edge(events_at_x[i-1],events_at_x[i]);
            }
        }
    }

    void y_struct_info(Y_structure_element e1,
                       std::ostream& out) {
        if(!e1.one_event_known) {
            out << "dummy node with id " << e1.x_index << std::endl;
        } else {
            if(e1.x_type==MINUS_INFTY) {
                out << "Point at -infty" << std::endl;
            } else if(e1.x_type==PLUS_INFTY) {
                out << "point at +infty" << std::endl;
            } else {
                out << "point at index" << e1.x_index;
                if(e1.y_type==MINUS_INFTY) {
                    out << " y-coor: -infty" << std::endl;
                } else if(e1.y_type==PLUS_INFTY) {
                    out << " y-coor: +infty" << std::endl;
                } else {
                    out << " y-id: " << e1.y_index << std::endl;
                }

            }
        }
    }
    
    void handle_edge(Y_structure_element& e1,
                     Y_structure_element& e2) {
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Y-STRUCT: " << std::endl;
        for(typename Y_structure::iterator it=y_structure.begin();
            it!=y_structure.end();it++) {
            y_struct_info(*it,CGAL_ACK_DEBUG_PRINT);
        }
        CGAL_ACK_DEBUG_PRINT << "Y-STRUCT done" << std::endl;
      
        CGAL_ACK_DEBUG_PRINT << "handle edge..." << std::flush;
        CGAL_ACK_DEBUG_PRINT << "info for e1: ";
        y_struct_info(e1,CGAL_ACK_DEBUG_PRINT);
        CGAL_ACK_DEBUG_PRINT << "info for e2: ";
        y_struct_info(e2,CGAL_ACK_DEBUG_PRINT);
#endif
*/
        if(! e1.one_event_known) {
            int id = e1.x_index;
            typename Y_structure::iterator it=y_structure.begin();
            while(it!=y_structure.end()) {
                if((! it->one_event_known) && it->x_index==id) {
                    //it=y_structure.erase(it);
                    //it=y_structure.insert(it,e2);
                    *it=e2;
                }
                it++;
            }
        } else if(! e2.one_event_known) {
            int id = e2.x_index;
            typename Y_structure::iterator it=y_structure.begin();
            while(it!=y_structure.end()) {
                if((! it->one_event_known) && it->x_index==id) {
                    //it=y_structure.erase(it);
                    //it=y_structure.insert(it,e1);
                    *it=e1;
                }
                it++;
            }
        } else {
            CGAL_assertion(e1.x_type!=e2.x_type || 
                           e1.x_type==FINITE || 
                           e2.x_type==FINITE);
            if(e2.x_type==MINUS_INFTY || e1.x_type==PLUS_INFTY) {
                handle_edge(e2,e1);
                return;
            } else if(e1.x_type==FINITE && e2.x_type==FINITE) {
                CGAL_assertion(e1.x_index!=e2.x_index);
                if(e1.x_index>e2.x_index) {
                    handle_edge(e2,e1);
                    return;
                }
            }
            int left_stripe = (e1.x_type==MINUS_INFTY) ? -1 : e1.x_index;
            int right_stripe = (e2.x_type==PLUS_INFTY) 
                ? static_cast<int>(stripe_values.size()-1) : e2.x_index;
            for(int i=left_stripe+1;i<right_stripe;i++) {
                pre_vert_lines[i].number_of_non_event_roots++;
            }
            if(e1.x_type==FINITE) {
                switch(e1.y_type) {
                case(PLUS_INFTY) : {
                    pre_vert_lines[e1.x_index].asym_right_plus++;
                    break;
                }
                case(MINUS_INFTY) : {
                    pre_vert_lines[e1.x_index].asym_right_minus++;
                    break;
                }
                case(FINITE): {
                    pre_vert_lines[e1.x_index].
                        event_points[e1.y_index].incident_right++;
                    break;
                }
                }
            }
            if(e2.x_type==FINITE) {
                switch(e2.y_type) {
                case(PLUS_INFTY) : {
                    pre_vert_lines[e2.x_index].asym_left_plus++;
                    break;
                }
                case(MINUS_INFTY) : {
                    pre_vert_lines[e2.x_index].asym_left_minus++;
                    break;
                }
                case(FINITE): {
                    pre_vert_lines[e2.x_index].
                        event_points[e2.y_index].incident_left++;
                    break;
                }
                }
            }   
        }
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "handle edge done" << std::endl;
#endif
*/
      
    }

    void sweep() {
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Init.." << std::flush;
#endif
*/
        start_sweep();
      
/*
#if CGAL_ACK_DEBUG_FLAG
        for(typename Y_structure::iterator it=y_structure.begin();
            it!=y_structure.end();
            it++) {
            y_struct_info(*it,CGAL_ACK_DEBUG_PRINT);
          }

        CGAL_ACK_DEBUG_PRINT << "End of y-struct" << std::endl;
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
*/
        for(int i=0;i<static_cast<int>(x_structure.size());i++) {
/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Coordinate.." 
                                 << CGAL::to_double(x_structure[i]) 
                                 << ", index " << i << std::endl;
            CGAL_ACK_DEBUG_PRINT <<  << "Y-struct before x" << std::endl;
            for(typename Y_structure::iterator it=y_structure.begin();
                it!=y_structure.end();
                it++) {
                y_struct_info(*it,CGAL_ACK_DEBUG_PRINT);
            }
            
            CGAL_ACK_DEBUG_PRINT << "End of y-struct" << std::endl;
#endif
*/
            sweep_at_x_coordinate(i);
        }
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "Terminate.." << std::flush;
#endif
        end_sweep();
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
    }

    Status_line_1 create_event_line(Curve_analysis_2& D,int i) {
        Algebraic_real_1 xval = sh_disc_roots[i];
        Bitstream_traits traits(Bitstream_coefficient_kernel(kernel(),xval));
        int number_of_events 
            = static_cast<int>(pre_vert_lines[i].event_points.size());
        int number_of_roots 
            = pre_vert_lines[i].number_of_non_event_roots+number_of_events; 
        Polynomial_2 sh_pol_with_correct_degree 
            = CGAL::internal::poly_non_vanish_leading_term(kernel(),sh_pol,xval);
        Bitstream_descartes descartes(CGAL::internal::Backshear_descartes_tag(),
                                      sh_pol_with_correct_degree,
                                      number_of_roots,
                                      number_of_events,
                                      pre_vert_lines[i],
                                      traits);
        typename Status_line_1::Arc_container arc_container;
        for(int j=0;j<descartes.number_of_real_roots();j++) {
            if(descartes.is_certainly_multiple_root(j)) {
                int n = static_cast<int>(pre_vert_lines[i].event_points.size());
                int k=0;
                while(k<n) {
                    if((pre_vert_lines[i].lower_bound(k)
                        <=descartes.right_bound(j)) &&
                       (pre_vert_lines[i].upper_bound(k)
                        >=descartes.left_bound(j))) {
                        break;
                    }
                    else {
                        k++;
                    }
                }
                CGAL_assertion(k<n);
                int left_arcs=pre_vert_lines[i].event_points[k].incident_left;
                int right_arcs=pre_vert_lines[i].event_points[k].incident_right;
                arc_container.push_back(std::make_pair(left_arcs,right_arcs));
            } else {
                arc_container.push_back(std::make_pair(1,1));
            }
        }
        Status_line_1 ev(xval,i,D,
                         pre_vert_lines[i].num_arcs_left(),
                         pre_vert_lines[i].num_arcs_right(),
                         arc_container);
        ev.set_isolator(descartes);
        ev._set_number_of_branches_approaching_infinity
            (std::make_pair(pre_vert_lines[i].asym_left_minus,
                            pre_vert_lines[i].asym_right_minus),
             std::make_pair(pre_vert_lines[i].asym_left_plus,
                            pre_vert_lines[i].asym_right_plus));
        return ev;
    }

    Algebraic_kernel_with_analysis_2* kernel() const {
        return _m_kernel;
    }

    Algebraic_kernel_with_analysis_2* _m_kernel;
  
    Curve_analysis_2 C;

    Integer s;

    Polynomial_2 pol, sh_pol, der_sh_pol,sh_der_sh_pol;

    Root_container sh_disc_roots,x_structure;
    std::vector<CGAL::internal::Three_valued> x_structure_info;
    
    std::vector<std::vector<int> > sh_ev_indices;

    std::vector<Bound> stripe_values;                   

    Bound far_left, far_right, y_in_box;

    Y_structure y_structure;

    std::vector<Sh_ev_line_info> pre_vert_lines;

    int x_extreme_index_counter;

    std::vector<Status_line_1> sh_intermediate_lines;

    bool use_primitive_curve;

    bool disc_roots_computed,sh_disc_roots_computed;

};
 
} //namespace CGAL

#endif
