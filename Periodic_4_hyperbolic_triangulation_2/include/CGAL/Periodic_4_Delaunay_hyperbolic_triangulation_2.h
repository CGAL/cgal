// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/candidate-packages/Triangulation_2/include/CGAL/Delaunay_triangulation_2.h $
// $Id: Delaunay_triangulation_2.h 57509 2010-07-15 09:14:09Z sloriot $
//
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_PERIODIC_4_DELAUNAY_HYPERBOLIC_TRIANGULATION_2_H
#define CGAL_PERIODIC_4_DELAUNAY_HYPERBOLIC_TRIANGULATION_2_H

#include <vector>

// I needed to use a non-reserved word for this thing, there were problems when I (stupidly) used TYPE
#define CLEARLY_MY_TYPE 2 // 1 = Cross, 2 = Diametric
#define USE_TEST_THINGS 1 // 1 = True, 0 = False

#if USE_TEST_THINGS == 1
#include <CGAL/Periodic_4_hyperbolic_triangulation_2.h>
#else
#include <CGAL/Delaunay_triangulation_2.h>
#endif

/*
#if CLEARLY_MY_TYPE == 1
#include <CGAL/Cross_translations.h>
#else
#include <CGAL/Diametric_translations.h>
#endif
*/
//#include <CGAL/Periodic_2_hyperbolic_triangulation_dummy.h>

namespace CGAL {
    
    template < class Gt,
    class Tds =  Triangulation_data_structure_2 <
    Triangulation_vertex_base_2<Gt> > >
#if USE_TEST_THINGS == 1
    class Periodic_4_Delaunay_hyperbolic_triangulation_2 : public Periodic_4_hyperbolic_triangulation_2<Gt,Tds>
#else
    class Periodic_4_Delaunay_hyperbolic_triangulation_2 : public Delaunay_triangulation_2<Gt,Tds>
#endif
    {
    public:
        typedef Periodic_4_Delaunay_hyperbolic_triangulation_2<Gt, Tds>     Self;
#if USE_TEST_THINGS == 1
        typedef Periodic_4_hyperbolic_triangulation_2<Gt,Tds>               Base;
#else
        typedef Delaunay_triangulation_2<Gt,Tds>                            Base;
#endif
        typedef typename Base::Vertex_handle                                Vertex_handle;
        typedef typename Base::Face_handle                                  Face_handle;
        typedef typename Base::Edge                                         Edge;
        typedef typename Base::Locate_type                                  Locate_type;
        
        typedef typename Base::Finite_edges_iterator                        Finite_edges_iterator;
        typedef typename Base::Face_circulator                              Face_circulator;
        
        typedef Gt                                                          Geom_traits;
        typedef typename Geom_traits::FT                                    FT;
        typedef typename Geom_traits::Point_2                               Point_2;
        typedef typename Geom_traits::Segment_2                             Segment;
        
        //void insert_dummy_points(std::vector<typename GT::Point_2>& );
        
        void Set_recursion_depth(int new_depth) {
            this->recursion_depth = new_depth;
        }
        
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
        using Base::tds;
#endif
  
  /*      
        Periodic_4_Delaunay_hyperbolic_triangulation_2(const Gt& gt = Gt())
#if USE_TEST_THINGS == 1
        : Periodic_4_hyperbolic_triangulation_2<Gt,Tds>(gt)
#else
        : Delaunay_triangulation_2<Gt,Tds>(gt)
#endif
        {
            recursion_depth = 0;
        }
        
        Periodic_4_Delaunay_hyperbolic_triangulation_2(const Periodic_4_Delaunay_hyperbolic_triangulation_2<Gt,Tds> &tr)
#if USE_TEST_THINGS == 1
        : Periodic_4_hyperbolic_triangulation_2<Gt,Tds>(tr)
#else
        : Delaunay_triangulation_2<Gt,Tds>(tr)
#endif
        { CGAL_triangulation_postcondition( this->is_valid() ); }
        
        
        */  

        /************************ INSERT OVERLOADS **************************/
        
      /*  
        
        Vertex_handle insert(const Point_2  &p,
                             Face_handle start = Face_handle() )
        {
            
            std::vector<Point_2> copies;

            Vertex_handle vh = Base::insert(p, start);
            
            recursive_translate(g, copies, p, recursion_depth);

            for(int i = 0; i < copies.size(); i++) {
                vh = Base::insert(copies[i], start);
            }
            
            return vh;
        }
        
        
        Vertex_handle insert(const Point_2& p,
                             typename Base::Locate_type lt,
                             Face_handle loc, int li )
        {
            return Base::insert(p, lt, loc, li);
        }
        
        
        template < class InputIterator >
        std::ptrdiff_t
        insert(InputIterator first, InputIterator last)
        {
            return Base::insert(first, last);
        }
        
    
    */

        Object
        dual(const Finite_edges_iterator& ei) const {
            return this->dual(*ei);
        }
        
        Object
        dual(const Edge &e) const {
            // default implementation
            //assert(false);
            return make_object(e);
        }
        
        // Clears the triangulation and initializes it again.
        void clear() {
            tds().clear();
            init_tds();
        }
        
    private:
        
        /*

        void recursive_translate(
#if CLEARLY_MY_TYPE == 1 
                                 Cross_translations<Gt> g,
#else
                                 Diametric_translations<Gt> g,
#endif
                                 std::vector<Point_2>& points,
                                 Point_2 o,
                                 int depth,
                                 int start = -1,
                                 int end = -1)
        {
            if (depth > 0) {

                if (start == -1 && end == -1) {
                    int start = points.size();
                    int end   = start + 8;
                    
                    Point_2 subject = o;

                    // Add points in the order indicated by the group -- not necessary, but seems logical
                    points.push_back( g.a().DoAction(subject) );
                    points.push_back( g.b().inverse().DoAction(subject) );
                    points.push_back( g.c().DoAction(subject) );
                    points.push_back( g.d().inverse().DoAction(subject) );
                    points.push_back( g.a().inverse().DoAction(subject) );
                    points.push_back( g.b().DoAction(subject) );
                    points.push_back( g.c().inverse().DoAction(subject) );
                    points.push_back( g.d().DoAction(subject) );
                
                    recursive_translate(g, points, o, depth - 1, start, end);
                }
                else
                {
                    int my_start = points.size();
                    int my_end   = start;
                
                    for (int i = start; i <= end; i++){
                        Point_2 subject = points[i];
                        points.push_back( g.a().DoAction(subject) );
                        points.push_back( g.b().inverse().DoAction(subject) );
                        points.push_back( g.c().DoAction(subject) );
                        points.push_back( g.d().inverse().DoAction(subject) );
                        points.push_back( g.a().inverse().DoAction(subject) );
                        points.push_back( g.b().DoAction(subject) );
                        points.push_back( g.c().inverse().DoAction(subject) );
                        points.push_back( g.d().DoAction(subject) );
                        my_end += 8;
                    }
                
                    recursive_translate(g, points, o, depth - 1, my_start, my_end);
                }
            }
            
        }
        
        */



        // Initializes the triangulation data structure
        void init_tds() {
            this->_infinite_vertex = tds().insert_first();
        }
        
        void paste_together_opposite_sides(const std::vector<Vertex_handle>& on_vertex, const std::vector<Vertex_handle>& on_boundary);
        
        // The following variable defines how many copies to insert on the Poincaré disk (for visualisation purposes only!)
        int recursion_depth;

/*

#if CLEARLY_MY_TYPE == 1 // This means that we must do Cross-translations (alternative identification)
        Cross_translations<Gt> g;
#else
        Diametric_translations<Gt> g;
#endif

*/


    };
    
} // namespace CGAL

//#include <Periodic_2_hyperbolic_triangulation_dummy.h>

#endif // CGAL_PERIODIC_4_DELAUNAY_HYPERBOLIC_TRIANGULATION_2_H
