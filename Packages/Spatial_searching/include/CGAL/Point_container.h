// Copyright (c) 2002 Utrecht University (The Netherlands).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)


// custom point container
#ifndef CGAL_POINT_CONTAINER_H
#define CGAL_POINT_CONTAINER_H
#include <list>
#include <vector>
#include <functional>
#include <algorithm>
#include <CGAL/Kd_tree_rectangle.h>
namespace CGAL {

  template <class SearchTraits> class Point_container {

  private:
    typedef typename SearchTraits::Point_d Point_d;
    typedef std::vector<Point_d*> Point_vector;
    typedef Point_container<Point_d> Self;

  public:
   typedef typename SearchTraits::FT FT;
   typedef std::list<Point_d*> Point_list; 
    typedef typename Point_list::iterator iterator;

  private:
    Point_list p_list; // list of pointers to points
    int built_coord;    // a coordinate for which the pointer list is built
    unsigned int the_size;
    Kd_tree_rectangle<SearchTraits> bbox;       // bounding box, i.e. rectangle of node
    Kd_tree_rectangle<SearchTraits> tbox;       // tight bounding box, 
				      // i.e. minimal enclosing bounding
	                	      // box of points
	                	    
  public:

    inline const Kd_tree_rectangle<SearchTraits>& bounding_box() const { return bbox; }

    inline const Kd_tree_rectangle<SearchTraits>& tight_bounding_box() const 
    { return tbox; }

    inline int dimension() const { return bbox.dimension(); } 

    inline int built_coordinate() const { return built_coord; } 

    // coordinate of the maximal span
    inline int max_span_coord() const { return bbox.max_span_coord(); }

    // coordinate of the maximal tight span
    inline int max_tight_span_coord() const { return tbox.max_span_coord(); }

    inline FT  max_span_lower() const 
	{ return bbox.min_coord(max_span_coord());}

    inline FT  max_tight_span_lower() const {
      return tbox.min_coord(max_tight_span_coord());}

    inline FT  max_span_upper() const 
	{ return bbox.max_coord(max_span_coord());}

    inline FT  max_tight_span_upper() const {
      return tbox.max_coord(max_tight_span_coord());}

    inline FT max_spread() const 
	{ return  max_span_upper() -  max_span_lower(); }

    inline FT max_tight_spread() const {
      return  max_tight_span_upper() -  max_tight_span_lower(); }


	int max_tight_span_coord_balanced(FT Aspect_ratio) const {
		int cut_dim(-1);
		FT max_spread_points(FT(-1));
		FT max_length=max_spread();  // length of longest side of box
		int dim=dimension();
		for (int d=0; d<dim; d++) {
			FT length=bbox.max_coord(d)-bbox.min_coord(d);
                     
		        if (FT(2)*max_length/length <= Aspect_ratio) {
			        FT spread=tbox.max_coord(d)-tbox.min_coord(d);
                                
			        if (spread > max_spread_points) {
				        max_spread_points = spread;
				        cut_dim = d;
			        }
                        }
		}
		// assert(cut_dim >= 0);
		return cut_dim;
	}

	FT max_span_upper_without_dim(int d) const {
		FT max_span(FT(0));
        	int dim=dimension();
		for (int i=0; i<dim; i++) {
			FT span = bbox.max_coord(i)-bbox.min_coord(i);
			if (d != i && span > max_span) max_span=span;
		}
		return max_span;
	}

	FT balanced_fair(int d, FT Aspect_ratio) {
	  	FT small_piece = 
		max_span_upper_without_dim(d) / Aspect_ratio;
	  	FT low_cut = 
		bbox.min_coord(d) + small_piece; // lowest legal cut;
	  	FT high_cut = 
		bbox.max_coord(d) - small_piece; //highest legal cut;
	  	// assert (high_cut >= low_cut);
        	FT split_value = median(d);
		if (split_value < low_cut) split_value=low_cut;
		if (split_value > high_cut) split_value=high_cut;
		return split_value;
	}

	FT balanced_sliding_fair(int d, FT Aspect_ratio) {
		FT small_piece = 
		max_span_upper_without_dim(d) / Aspect_ratio;
		FT low_cut = 
		bbox.min_coord(d) + small_piece; // lowest legal cut;
		FT high_cut = 
		bbox.max_coord(d) - small_piece; //highest legal cut;
		// assert (high_cut >= low_cut);
                FT split_value = median(d);
		FT max_span_lower = tbox.min_coord(d);
		FT max_span_upper = tbox.max_coord(d);
		if (split_value < low_cut) split_value= max_span_lower; 
		if (split_value > high_cut) split_value = max_span_upper; 
		return split_value;
	}

    //  points
    inline unsigned int size() const {
    	return the_size;
    }
    
    inline void set_size() {the_size=p_list.size(); }
    
    inline typename Point_list::const_iterator begin() const {
      return p_list.begin();
    }

    inline typename Point_list::const_iterator end() const {
      return p_list.end();
    }
    
    // building the container from a sequence of points
    template <class InputIterator>
    Point_container(const int d, InputIterator begin, InputIterator end) :
       bbox(d), tbox(d)  {

        

      bbox = Kd_tree_rectangle<SearchTraits>(d, begin, end);
      tbox = bbox;

      // build list 
      InputIterator it;
      for (it=begin; it != end; ++it) p_list.push_back(&(*it));

      built_coord = max_span_coord();
      set_size();
    }

	// building an empty container 
	Point_container(const int d) :
	the_size(0), bbox(d), tbox(d)  {}

	void swap(Point_container<SearchTraits>& c) {

		swap(p_list,c.p_list);

        // Borland generates strange compile errors
		// swap(built_coord,c.built_coord);
		// swap(bbox,c.bbox);
		// swap(tbox,c.tbox);


                // work-around
                int h=built_coord;
                built_coord = c.built_coord;
                c.built_coord = h;

                // work-around
                Kd_tree_rectangle<SearchTraits> h_bbox(bbox);
                bbox = c.bbox;
                c.bbox = h_bbox;

                // work-around
                Kd_tree_rectangle<SearchTraits> h_tbox(tbox);
                tbox = c.tbox;
                c.tbox = h_tbox;
                
                //work-around
                h=the_size;
                the_size = c.the_size;
                c.the_size = h;
                
	}

       

    void recompute_tight_bounding_box() {
		tbox.update_from_point_pointers(p_list.begin(),
		     p_list.end(),p_list.empty());
	}


      // note that splitting is restricted to the built coordinate
      template <class Separator>
      void split(Point_container<SearchTraits>& c, Separator& sep,  
	bool sliding=false) {

	assert(dimension()==c.dimension());
		
        Point_list l_lower, l_upper;

        c.bbox=bbox;
        // bool test_validity=false;

        const int split_coord = sep.cutting_dimension();
        const FT cutting_value = sep.cutting_value();

        built_coord=split_coord;
	c.built_coord=split_coord;
		
	
	iterator pt=p_list.begin();
			
	typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
	typename SearchTraits::Cartesian_const_iterator_d ptit;

	for (; (pt != p_list.end()); ++pt) {
                        
	  ptit = construct_it(*(*pt));
	if ( *(ptit+split_coord) < cutting_value) 
			l_lower.push_back (*pt);
		else
			l_upper.push_back (*pt);
	};
	
	if (sliding) { // avoid empty lists 
		if (l_lower.empty()) {
		  iterator pt_min=l_upper.begin();
		  FT min_value=bbox.max_coord(built_coord);
		  for (pt=l_upper.begin(); (pt != l_upper.end()); ++pt) {
		    ptit = construct_it((*(*pt)));	
				if ( *(ptit+split_coord) < min_value) {
				  min_value= *(ptit+split_coord);
					pt_min=pt;
				}
			}
			l_lower.splice(l_lower.end(), l_upper, pt_min);
		}
		if (l_upper.empty()) {
		  iterator pt_max=l_lower.begin();
		  FT max_value=bbox.min_coord(built_coord);
		  for (pt=l_lower.begin(); (pt != l_lower.end()); ++pt) {
		    ptit = construct_it((*(*pt)));	
				if ( *(ptit+split_coord) > max_value) {
					max_value= *(ptit+split_coord) ;
					pt_max=pt;
				}
			}
			l_upper.splice(l_upper.end(), l_lower, pt_max);
		}
        }
	

	p_list.clear();
        c.p_list.clear();
        p_list.splice(p_list.end(),l_upper);
        c.p_list.splice(c.p_list.end(),l_lower);
        
		
	// adjusting boxes
	bbox.set_lower_bound(split_coord, cutting_value);
	tbox.update_from_point_pointers(p_list.begin(),
	p_list.end(),p_list.empty());
	c.bbox.set_upper_bound(split_coord, cutting_value);
	c.tbox.update_from_point_pointers(
	c.p_list.begin(),
	c.p_list.end(),c.p_list.empty());
        
        c.set_size();
        set_size();
        
       
    }



    template <class SearchTraits2, class Value>
    struct comp_coord_val {
      
    private:
      Value coord;   
      
      typedef typename SearchTraits2::Point_d Point_d;
    public:
      comp_coord_val (const Value& coordinate) : coord(coordinate) {}
      
      bool operator() (const Point_d *a, const Point_d *b) {
	typename SearchTraits2::Construct_cartesian_const_iterator_d construct_it;
	typename SearchTraits2::Cartesian_const_iterator_d ait = construct_it(*a),
	  bit = construct_it(*b);
	return *(ait+coord) < *(bit+coord);
    }
  };
  

  FT median(const int split_coord) {
      
    #ifdef CGAL_CFG_RWSTD_NO_MEMBER_TEMPLATES
        Point_vector p_vector;
    	std::copy(p_list.begin(), p_list.end(), std::back_inserter(p_vector));
    	std::sort(p_vector.begin(),p_vector.end(),comp_coord_val<SearchTraits,int>(split_coord));
    	p_list.clear();
    	std::copy(p_vector.begin(), p_vector.end(), std::back_inserter(p_list));
    #else
        p_list.sort(comp_coord_val<SearchTraits,int>(split_coord));
    #endif 
      
      iterator median_point_ptr = p_list.begin();
      for (unsigned int i = 0; i < the_size/2-1; i++, 
		   median_point_ptr++) {}
      
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
      typename SearchTraits::Cartesian_const_iterator_d mpit = construct_it((*(*median_point_ptr)));
      FT val1= *(mpit+split_coord);
      median_point_ptr++;
      mpit = construct_it((*(*median_point_ptr)));
      FT val2= *(mpit+split_coord);
      
      
      
      return (val1+val2)/FT(2); 
    };


    inline bool empty() const { return the_size == 0;}
     
     

  private:
    explicit Point_container() {} // disable default constructor
  
  };

template <class Point>
    std::ostream& operator<< (std::ostream& s, Point_container<Point>& c) {
    s << "Points container of size " << the_size << "\n cell:";
    s << bbox; // bbox.print(s);
    s << "\n minimal box enclosing points:"; s << tbox; 
    return s;
  }

} // namespace CGAL

#endif // CGAL_POINT_CONTAINER_H


