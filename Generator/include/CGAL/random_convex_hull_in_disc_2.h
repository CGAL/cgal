// Copyright (c) 2014
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
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
// Author(s)     : Remy Thomasse  <remy.thomasse@inria.fr>


#ifndef  CGAL_RANDOM_CONVEX_HULL_DISC_H
#define  CGAL_RANDOM_CONVEX_HULL_DISC_H
#include <iostream>
#include <list>
#include <algorithm>

#ifndef Q_MOC_RUN
#include <boost/random.hpp>
#endif

#include <cmath>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/convex_hull_traits_2.h>
#include <CGAL/function_objects.h>
namespace CGAL{
    namespace internal{
    template<class P >
        struct compare_points_angle{
            bool operator()(const P&  p, const P&  q ){
                typedef typename Kernel_traits<P>::Kernel Traits;
                Traits ch_traits;
                typedef  typename Traits::Left_turn_2  Left_turn;
                Left_turn left_turn = ch_traits.left_turn_2_object();

                if ( to_double(p.y())>0 )
                {
                    if (to_double(q.y())>0) return left_turn(ORIGIN,p,q);
                    else return false;

                }
                else {
                    if (to_double(q.y())>0) return true;
                    else return left_turn(ORIGIN,p,q);

                }
            }
        };


    //////////////////////////////////////
    template <class P, class GEN>
    void generate_points_annulus(long n,double a, double b,double small_radius, double big_radius,std::list<P> & l,GEN & gen){ //generate n points between a and b
        if (n>1){
            boost::random::binomial_distribution<long>di(n,.5);
            boost::random::variate_generator<GEN &,boost::binomial_distribution<long> >d(gen,di);
            long nb=d();
            generate_points_annulus(nb,a,(a+b)/2.0,small_radius,big_radius,l,gen);
            generate_points_annulus(n-nb,(a+b)/2.0,b,small_radius,big_radius, l,gen);
        }
        if (n==1)//generation of a point
        {
            boost::random::uniform_real_distribution<double> gd(small_radius*small_radius/(big_radius*big_radius),1);
            
            boost::random::uniform_real_distribution<double> hd(a,b);
            boost::random::variate_generator<GEN &,boost::random::uniform_real_distribution<double> >h(gen,hd);
            boost::random::variate_generator<GEN &,boost::random::uniform_real_distribution<double> >g(gen,gd);
            double alpha=h();
            double r=big_radius*std::sqrt(g());
            //typedef Creator_uniform_2<typename Kernel_traits<P>::Kernel::FT,P> Creator;
            typedef Creator_uniform_2<double, P> Creator;
            Creator creator;
            typedef typename Creator::argument_type T;
            l.push_back(creator(T(r*cos(alpha)) ,T(r*std::sin(alpha))));
            
        }
        
    }
    
    
    template <class P>
    void Cyclic_increment_iterator(typename std::list<P>::iterator & it,std::list<P> & l){
        ++it;
        if (it==l.end()) 
        {
            it=l.begin();
        }
        
    }
    //////////////////////////////////////////////////////////////////////////////
    template<class P, class Traits >
    void Graham_without_sort_2( std::list<P> &l, const Traits& ch_traits){
        if (l.size()>3){
            typedef  typename Traits::Left_turn_2  Left_turn;
            Left_turn left_turn = ch_traits.left_turn_2_object();
            typename std::list<P>::iterator  pmin=l.begin();
            for (typename std::list<P>::iterator it=l.begin(); it!=l.end(); ++it) {
                if ((*pmin).x()>(*it).x())
                {
                    pmin=it;
                }
        }//*pmin is the extremal point on the left
        typename std::list<P>::iterator u=pmin;
        typename std::list<P>::iterator u_next=u;
        Cyclic_increment_iterator(u_next,l);
        
        
        typename std::list<P>::iterator u_next_next=u_next;
        Cyclic_increment_iterator(u_next_next,l);
        
        while (u_next !=pmin){

            if (left_turn(*u,*u_next,*u_next_next)){
                Cyclic_increment_iterator(u,l);
                Cyclic_increment_iterator(u_next,l);
                Cyclic_increment_iterator(u_next_next,l);
            }
            else{
                u_next=l.erase(u_next);
                if (u_next==l.end()) u_next=l.begin();
                if(u !=pmin){
                    u_next=u;
                    if(u==l.begin()){u=l.end();}
                    --u;
                }
                else{
                    u_next_next=u_next;
                    Cyclic_increment_iterator(u_next_next,l);
                }
            }
        }
    }

}



    //////////////////////////////////////////////////////////////////////////////
    template<class P,class GEN>
void random_convex_hull_in_disc_2(std::size_t n,  double radius, std::list<P> & l,GEN & gen, bool fast=true ){
    CGAL_precondition( n >= 3);
    typedef typename Kernel_traits<P>::Kernel K;
    std::size_t simulated_points=0;
    std::size_t generated_points=0;
    do
        { //Initialisation
            //std::size_t init=std::min( (std::size_t)100,n-simulated_points );
            std::size_t init=std::min( static_cast<std::size_t>(100), n-simulated_points );
            generate_points_annulus(init,-CGAL_PI, CGAL_PI,0,to_double(radius),l,gen);
            
            simulated_points+=init;
            generated_points+=init;
            Graham_without_sort_2(l,K());
        } while ((bounded_side_2(l.begin(),l.end(),P  (0,0),K())!=ON_BOUNDED_SIDE)&&(simulated_points<n)); //initialisation such that 0 in P_n
        std::size_t T=n;
        //if (!fast)  T=(size_t)std::floor(n/std::pow(log(n),2));
        if (!fast) T=static_cast<std::size_t>( std::floor( n/std::pow(log(n),2) ) );
        while (simulated_points<n)
        {
            //l is a list coming from a convex hull operation. we are moving the points s.t the angles are from -pi to pi.
            {
                typename std::list<P>::iterator it=l.begin();
                while(to_double((*it).y())>0){
                    l.push_back(*it);
                    l.pop_front();
                    it=l.begin();
                }
                it=l.end();
                --it;//last element
                while(to_double((*it).y())<0){
                    l.push_front(*it);
                    l.pop_back();
                    it=l.end();
                    --it;//last element
                }
                
            }
            double squared_radius=radius*radius;
            double squared_small_radius=squared_radius;
            
            {

                P zero(0,0);
                typename std::list<P>::iterator it=l.begin();
                typename std::list<P>::iterator it2=++it;
                for(;it!=l.end();++it,Cyclic_increment_iterator(it2,l)){ //computation of  annulus
                    typename K::Segment_2 s(*it,*it2);
                    double temp=to_double(squared_distance(s,zero));
                    if (squared_small_radius>temp)  squared_small_radius=temp;
                }
            }//squared_small_radius=squared small radius of the annulus
            
            
            double p_disc=squared_small_radius/squared_radius;
            std::size_t nb;
            if (simulated_points< T){nb=std::min(simulated_points,n-simulated_points);}
            else {nb=std::min(T,n-simulated_points); }
            boost::random::binomial_distribution<long> dbin(nb,p_disc);
            boost::random::variate_generator<GEN&, boost::random::binomial_distribution<long> >bin(gen,dbin);
            
              //How many points are falling in the small disc and wont be generated:
            long k_disc=bin();
            simulated_points+=k_disc;
            
            std::list<P> m;
            internal::generate_points_annulus(nb-k_disc,-CGAL_PI, CGAL_PI,std::sqrt(squared_small_radius),radius,m,gen);
            l.merge(m,internal::compare_points_angle<P>());
            generated_points+=nb-k_disc;
            simulated_points+=nb-k_disc;
            m.clear();
            Graham_without_sort_2(l,K());
        }
    }
} //namespace CGAL::internal

    ///

    template <class OutputIterator, class Traits, class Generator>
void random_convex_hull_in_disc_2(std::size_t n, double radius, Generator & gen, OutputIterator it, const Traits & traits, bool fast=true)
{
    typedef Point_2<Traits> Points;
    std::list<Points> l;
    internal::random_convex_hull_in_disc_2(n,radius, l,gen,fast);
    // for(typename std::list<Points>::iterator i=l.begin();i!=l.end();++i)
    // {
    //     *it=*i;
    //     ++*it;
    // }
    std::copy(l.begin(),l.end(), it);

}

}//namespace CGAL
#endif