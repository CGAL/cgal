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
// file          : include/CGAL/Points_container.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


// custom point container
#ifndef CGAL_POINTS_CONTAINER_H
#define CGAL_POINTS_CONTAINER_H
#include <list>
#include <set>
#include <functional>
#include <algorithm>
#include <CGAL/Box.h>
namespace CGAL {

  template <class Item> class Points_container;

  template <class Item> struct comp_coord_val {
    int coord; 
    comp_coord_val(const int i) : coord(i) {}
    bool operator() (const Item *a, const Item *b) {
      return (*a)[coord] < (*b)[coord];
    }
  };

  template <class Item>
    std::ostream& operator<< (std::ostream& s, Points_container<Item>& c) {
    return c.print(s);
  }

  template <class Item> class Points_container {
  public:
    typedef std::list<Item*> Points_list; 

    typedef typename Kernel_traits<Item>::Kernel K;
    typedef typename K::FT NT;
    typedef Points_container<Item> Self;

  private:
    Points_list *p_list;// array of sorted lists of pointers to points
    int built_coord;    // a coordinate for which the pointer list is built
    //    Points points;// points container
    Box<NT> bbox;       // bounding box, i.e. cell of node
    Box<NT> tbox;       // tight bounding box, i.e. minimal enclosing bounding
	                // box of points

    struct build_max_span_list : public std::unary_function<Item, void> {
      Points_list *x;
      build_max_span_list(Points_list *p) :
	x(p) {}
      void operator() (Item& p) {
	x->push_back(&p);
      }
    };

  public:
    std::ostream& print(std::ostream& s) {
      s << "Points container of size " << size() << "\n cell:";
      s << bbox; // bbox.print(s);
      s << "\n minimal box enclosing points:"; s << tbox; // tbox.print(s);
      return s;
    };

    inline const Box<NT>& bounding_box() const { return bbox; }

    inline const Box<NT>& tight_bounding_box() const { return tbox; }

    inline const int dimension() const { return bbox.dimension(); } 

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

	NT max_span_upper_without_dim(int d) {
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
    inline unsigned int size() const { return p_list[built_coord].size(); }

    inline typename Points_list::const_iterator begin() const {
      return p_list[built_coord].begin();
    }

    inline typename Points_list::const_iterator end() const {
      return p_list[built_coord].end();
    }

    // building the container from a sequence of points
    template <class Iter>
    Points_container(const int d, Iter begin, Iter end) :
      p_list(new Points_list[d]), bbox(d), tbox(d) {

      

      bbox = Box<NT>(d, begin, end);
      tbox = bbox;
      std::for_each(begin, end,
		    build_max_span_list(p_list + max_span_coord()));


      p_list[max_span_coord()].sort(comp_coord_val<Item>(max_span_coord()));
      built_coord = max_span_coord();
    }

	// building an empty container 
	Points_container(int d) :
	p_list(new Points_list[d]), bbox(d), tbox(d) {}

	void swap(Points_container<Item>& c) {

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
                Box<NT> h_bbox(bbox);
                bbox = c.bbox;
                c.bbox = h_bbox;

                // work-around
                Box<NT> h_tbox(tbox);
                tbox = c.tbox;
                c.tbox = h_tbox;
	}

	void add_points_from_container(Points_container<Item>& c) {
	  // assert(built_coord==c.built_coord);
	  merge(p_list[built_coord], c.p_list[built_coord], 
		Less_lexicographically_d());
	}

    void recompute_tight_bounding_box() {
		tbox.update_from_point_pointers(p_list[built_coord].begin(),
		     p_list[built_coord].end(),p_list[built_coord].empty());
	}


      // note that splitting is restricted to the built coordinate
      template <class Separator>
      void split_container(Points_container<Item>& c, Separator* sep, 
	  bool sliding=false) {

	//	assert(dimension()==c.dimension());
		

        c.bbox=bbox;
        // bool test_validity=false;

        const int split_coord = sep->cutting_dimension();
        const NT cutting_value = sep->cutting_value();

		// if necessary prepare the coordinate and 
		// clear old built_coord.
		if (p_list[split_coord].empty()) {
			// copy p_list[built_coord] to p_list[split_coord]
			p_list[split_coord]=p_list[built_coord];
			// sort p_list[split_coord]
		p_list[split_coord].sort(comp_coord_val<Item>(split_coord));
	        // clear old built coord
		}
		built_coord=split_coord;
		c.built_coord=split_coord;
		

		for (int i = 0; i < dimension(); ++i) {
			// if (! p_list[i].empty()) {
			if (i==split_coord) {    
                // find iterator to split the list
				// avoid empty list by moving first 
				c.p_list[i].clear();
		typename Points_list::iterator pt;
                for (pt = p_list[i].begin();
				( (sep->side(*(*pt)) == ON_NEGATIVE_SIDE) 
				&&
				(pt != p_list[i].end())
					); ++pt) {}
                // avoid empty lists 
				if (sliding) {
			    	if (pt==p_list[i].begin()) pt++;
					if (pt==p_list[i].end())  pt--;
                };
				// move points on negative side to c.p_list
			    c.p_list[i].splice(c.p_list[i].end(), 
				p_list[i], p_list[i].begin(), pt); 
				
		    } // end of if
				else { 
				  p_list[i].clear(); c.p_list[i].clear(); 
				}
		} // end of for
		
		// adjusting boxes
		bbox.set_lower(split_coord, cutting_value);
		tbox.update_from_point_pointers(p_list[built_coord].begin(),
		p_list[built_coord].end(),p_list[built_coord].empty());
		c.bbox.set_upper(split_coord, cutting_value);
		c.tbox.update_from_point_pointers(
		c.p_list[c.built_coord].begin(),
		c.p_list[c.built_coord].end(),c.p_list[c.built_coord].empty());
        
        // assert(is_valid()); 
        // assert(c.is_valid());
	}

	NT median(const int split_coord) {
      if (p_list[split_coord].empty()) {
         // copy p_list[built_coord] to p_list[split_coord]
         p_list[split_coord]=p_list[built_coord];
         // sort p_list[split_coord]
         p_list[split_coord].sort(comp_coord_val<Item>(split_coord));
      }
      typename Points_list::iterator 
      median_point_ptr=p_list[split_coord].begin();
      for (unsigned int i = 0; i < p_list[split_coord].size()/2-1; i++, 
		   median_point_ptr++) {}
      NT val1=(*(*median_point_ptr))[split_coord];
      median_point_ptr++;
      NT val2=(*(*median_point_ptr))[split_coord];
      return (val1+val2)/2;
    }


    ~Points_container() { delete [] p_list; }

    inline bool empty() const { return size() == 0;}
     /*
     bool is_valid() {

	 assert(! p_list[built_coord].empty());

      // checking that all the lists are of the same size
     for (int i = 0; i < dimension(); ++i)
	 if (! p_list[i].empty()) assert(p_list[i].size() == size());
     
     // checking that all the points lie in the box specified

      for (int i = 0; i < dimension(); ++i) {
      if (! p_list[i].empty())   {
        typename Points_list::iterator p1=p_list[i].begin();
		if (! belongs(*(*p1), bbox)) {
			std::cout << "error: point " << (*(*p1)) 
						     << std::endl;
			std::cout << "error: does not belong to box " 
			<< bbox << std::endl;
		}
        assert(belongs(*(*p1), bbox));
        if  (p_list[i].size()>1) {
           typename Points_list::iterator p2 = p1;
           p2++;
	   for (; p2 != p_list[i].end(); p1++, p2++ ) {
                assert ((*(*p1))[i] <= (*(*p2))[i]);
                assert(belongs(*(*p2), bbox));
           }
        }
      // else std:: cout << "p_list is empty" << std:: endl;
      // end added by Hans
     } 
	 }
     return true;
  } */

     

  private:
    explicit Points_container() {} // disable default constructor
  
  };

} // namespace CGAL

#endif // CGAL_POINTS_CONTAINER_H


