// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Point_container.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


// custom point container
#ifndef CGAL_POINT_CONTAINER_H
#define CGAL_POINT_CONTAINER_H
#include <list>
#include <set>
#include <functional>
#include <algorithm>
#include <CGAL/Kd_tree_rectangle.h>
namespace CGAL {

  template <class Item> class Point_container;

  /* problem with platform sparc_SunOS-5.6_CC-5.30
  template <class Item> struct comp_coord_val {
    int coord; 
    comp_coord_val(const int i) : coord(i) {}
    bool operator() (const Item *a, const Item *b) {
      return (*a)[coord] < (*b)[coord];
    }
  }; */

  template <class Item>
    std::ostream& operator<< (std::ostream& s, Point_container<Item>& c) {
    return c.print(s);
  }

  template <class Item> class Point_container {
  public:
    typedef std::list<Item*> Point_list; 
    typedef typename Item::R::FT NT;
    typedef Point_container<Item> Self;

    

  private:
    Point_list p_list; // list of pointers to points
    int built_coord;    // a coordinate for which the pointer list is built
    Kd_tree_rectangle<NT> bbox;       // bounding box, i.e. rectangle of node
    Kd_tree_rectangle<NT> tbox;       // tight bounding box, 
				      // i.e. minimal enclosing bounding
	                	      // box of points
  public:
    std::ostream& print(std::ostream& s) {
      s << "Points container of size " << size() << "\n cell:";
      s << bbox; // bbox.print(s);
      s << "\n minimal box enclosing points:"; s << tbox; // tbox.print(s);
      return s;
    };

    inline const Kd_tree_rectangle<NT>& bounding_box() const { return bbox; }

    inline const Kd_tree_rectangle<NT>& tight_bounding_box() const 
    { return tbox; }

    inline int dimension() const { return bbox.dimension(); } 

    inline int built_coordinate() const { return built_coord; } 

    // coordinate of the maximal span
    inline int max_span_coord() const { return bbox.max_span_coord(); }

    // coordinate of the maximal tight span
    inline int max_tight_span_coord() const { return tbox.max_span_coord(); }

    inline NT  max_span_lower() const { return bbox.lower(max_span_coord());}

    inline NT  max_tight_span_lower() const {
      return tbox.lower(max_tight_span_coord());}

    inline NT  max_span_upper() const { return bbox.upper(max_span_coord());}

    inline NT  max_tight_span_upper() const {
      return tbox.upper(max_tight_span_coord());}

    inline NT max_spread() const 
	{ return  max_span_upper() -  max_span_lower(); }

    inline NT max_tight_spread() const {
      return  max_tight_span_upper() -  max_tight_span_lower(); }


	int max_tight_span_coord_balanced(NT Aspect_ratio) const {
		int cut_dim(-1);
		NT max_spread_points(-1.0);
		NT max_length=max_spread();  // length of longest side of box
		int dim=dimension();
		for (int d=0; d<dim; d++) {
			NT length=bbox.upper(d)-bbox.lower(d);
		        if (2.0*max_length/length <= Aspect_ratio) {
			        NT spread=tbox.upper(d)-tbox.lower(d);
			        if (spread > max_spread_points) {
				        max_spread_points = spread;
				        cut_dim = d;
			        }
                        }
		}
		// assert(cut_dim >= 0);
		return cut_dim;
	}

	NT max_span_upper_without_dim(int d) const {
		NT max_span(0.0);
        	int dim=dimension();
		for (int i=0; i<dim; i++) {
			NT span = bbox.upper(i)-bbox.lower(i);
			if (d != i && span > max_span) max_span=span;
		}
		return max_span;
	}

	NT balanced_fair(int d, NT Aspect_ratio) {
		NT small_piece = max_span_upper_without_dim(d) / Aspect_ratio;
		NT low_cut = bbox.lower(d) + small_piece; // lowest legal cut;
		NT high_cut = bbox.upper(d) - small_piece; //highest legal cut;
		// assert (high_cut >= low_cut);
        	NT split_value = median(d);
		if (split_value < low_cut) split_value=low_cut;
		if (split_value > high_cut) split_value=high_cut;
		return split_value;
	}

	NT balanced_sliding_fair(int d, NT Aspect_ratio) {
		NT small_piece = max_span_upper_without_dim(d) / Aspect_ratio;
		NT low_cut = bbox.lower(d) + small_piece; // lowest legal cut;
		NT high_cut = bbox.upper(d) - small_piece; //highest legal cut;
		// assert (high_cut >= low_cut);
                NT split_value = median(d);
		NT max_span_lower = tbox.lower(d);
		NT max_span_upper = tbox.upper(d);
		if (split_value < low_cut) split_value= max_span_lower; 
		if (split_value > high_cut) split_value = max_span_upper; 
		return split_value;
	}

    //  points
    inline unsigned int size() const { return p_list.size(); }

    inline typename Point_list::const_iterator begin() const {
      return p_list.begin();
    }

    inline typename Point_list::const_iterator end() const {
      return p_list.end();
    }
    
    // building the container from a sequence of points
    template <class Iter>
    Point_container(const int d, Iter begin, Iter end) :
      bbox(d), tbox(d) {

      

      bbox = Kd_tree_rectangle<NT>(d, begin, end);
      tbox = bbox;

      // build list 
      p_list.clear();
      Iter it;
      for (it=begin; it != end; ++it) p_list.push_back(&(*it));

      // p_list[max_span_coord()].sort(comp_coord_val<Item>(max_span_coord()));
      built_coord = max_span_coord();
    }

	// building an empty container 
	Point_container(int d) :
	bbox(d), tbox(d) {}

	void swap(Point_container<Item>& c) {

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
                Kd_tree_rectangle<NT> h_bbox(bbox);
                bbox = c.bbox;
                c.bbox = h_bbox;

                // work-around
                Kd_tree_rectangle<NT> h_tbox(tbox);
                tbox = c.tbox;
                c.tbox = h_tbox;
	}

	void add_points_from_container(Point_container<Item>& c) {
	  // assert(built_coord==c.built_coord);
	  merge(p_list, c.p_list); 
		// Less_lexicographically_d());
	}

    void recompute_tight_bounding_box() {
		tbox.update_from_point_pointers(p_list.begin(),
		     p_list.end(),p_list.empty());
	}


      // note that splitting is restricted to the built coordinate
      template <class Separator>
      void split_container(Point_container<Item>& c, Separator* sep,  
	bool sliding=false) {

	//assert(dimension()==c.dimension());
		
        Point_list l_lower, l_upper;

	l_lower.clear();
	l_upper.clear();

        c.bbox=bbox;
        // bool test_validity=false;

        const int split_coord = sep->cutting_dimension();
        const NT cutting_value = sep->cutting_value();

        built_coord=split_coord;
	c.built_coord=split_coord;
		
	
	typename Point_list::iterator pt=p_list.begin();
				
	for (; (pt != p_list.end()); ++pt) {
                        
	if ( (*(*pt))[split_coord] < cutting_value) 
			l_lower.push_back (*pt);
		else
			l_upper.push_back (*pt);
	};
	
	if (sliding) { // avoid empty lists 
		if (l_lower.empty()) {
			  typename Point_list::iterator pt_min=l_upper.begin();
			  NT min_value=bbox.upper(built_coord);
			  for (pt=l_upper.begin(); (pt != l_upper.end()); ++pt) {
				if ( (*(*pt))[split_coord] < min_value) {
					min_value=(*(*pt))[split_coord];
					pt_min=pt;
				}
			}
			l_lower.splice(l_lower.end(), l_upper, pt_min);
		}
		if (l_upper.empty()) {
			typename Point_list::iterator pt_max=l_lower.begin();
			NT max_value=bbox.lower(built_coord);
			for (pt=l_lower.begin(); (pt != l_lower.end()); ++pt) {
				if ( (*(*pt))[split_coord] > max_value) {
					max_value=(*(*pt))[split_coord];
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
        
        // assert(is_valid()); 
        // assert(c.is_valid());
    }

/*
template <class Item_>
	struct comp_coord_val {
    		int coord; 
    		comp_coord_val(const int i) : coord(i) {}
    		bool operator() (const Item_ *a, const Item_ *b) {
      		return (*a)[coord] < (*b)[coord];
    		}
  	};
*/ 
// replaced by

template <class Item_, class Value>
	struct comp_coord_val {
        
	private:
                Value coord;   

	public:
		comp_coord_val (const Value& coordinate) : coord(coordinate) {}

    		bool operator() (const Item_ *a, const Item_ *b) {
  	    		return (*a)[coord] < (*b)[coord];
    		}
  	};


      NT median(const int split_coord) {
        
      p_list.sort(comp_coord_val<Item,int>(split_coord));
      
      typename Point_list::iterator 
      median_point_ptr=p_list.begin();
      for (unsigned int i = 0; i < p_list.size()/2-1; i++, 
		   median_point_ptr++) {}
      NT val1=(*(*median_point_ptr))[split_coord];
      median_point_ptr++;
      NT val2=(*(*median_point_ptr))[split_coord];
      return (val1+val2)/2;
    };


    ~Point_container() { 
			p_list.clear(); 
    }

    inline bool empty() const { return p_list.size() == 0;}
     
     

  private:
    explicit Point_container() {} // disable default constructor
  
  };

} // namespace CGAL

#endif // CGAL_POINT_CONTAINER_H


