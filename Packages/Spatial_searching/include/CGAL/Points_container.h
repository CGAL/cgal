// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// package       : APSPAS
// revision      : 1.0 
// revision_date : 2001/06/12 
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
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
  };

  template <class Item> class Points_container {
  public:
    typedef std::list<Item*> Points_list; 

    typedef typename double NT; // Item::FT NT;
    typedef Points_container<Item> Self;
  private:
    Points_list *p_list; // array of sorted lists of pointers to points
    int built_coord;     // a coordinate for which the pointer list is built
    //    Points points;// points container
    Box<NT> bbox;        // bounding box, i.e. cell of node
    Box<NT> tbox;        // tight bounding box, i.e. minimal enclosing bounding box of points

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

	inline int built_coordinate() { return built_coordinate(); } 

    // coordinate of the maximal span
    inline int max_span_coord() { return bbox.max_span_coord(); }

    // coordinate of the maximal tight span
    inline int max_tight_span_coord() { return tbox.max_span_coord(); }

    inline NT  max_span_lower() { return bbox.lower(max_span_coord());}

    inline NT  max_tight_span_lower() {
      return tbox.lower(max_tight_span_coord());}

    inline NT  max_span_upper() { return bbox.upper(max_span_coord());}

    inline NT  max_tight_span_upper() {
      return tbox.upper(max_tight_span_coord());}

    inline NT max_spread() { return  max_span_upper() -  max_span_lower(); }

    inline NT max_tight_spread() {
      return  max_tight_span_upper() -  max_tight_span_lower(); }


	int max_tight_span_coord_balanced(NT Aspect_ratio) {
		int cut_dim = -1;
		NT max_spread_points = -1.0;
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
		assert(cut_dim >= 0);
		return cut_dim;
	}

	NT max_span_upper_without_dim(int d) {
		NT max_span=0.0;
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
		NT high_cut = bbox.upper(d) - small_piece; // hihgest legal cut;
		assert (high_cut >= low_cut);
        NT split_value = median(d);
		if (split_value < low_cut) split_value=low_cut;
		if (split_value > high_cut) split_value=high_cut;
		return split_value;
	}

	NT balanced_sliding_fair(int d, NT Aspect_ratio) {
		NT small_piece = max_span_upper_without_dim(d) / Aspect_ratio;
		NT low_cut = bbox.lower(d) + small_piece; // lowest legal cut;
		NT high_cut = bbox.upper(d) - small_piece; // hihgest legal cut;
		assert (high_cut >= low_cut);
                NT split_value = median(d);
		NT max_span_lower = tbox.lower(d);
		NT max_span_upper = tbox.upper(d);
		if (split_value < low_cut) split_value= max_span_lower; 
		if (split_value > high_cut) split_value = max_span_upper; 
		return split_value;
	}

    //  points
    inline unsigned int size() { return p_list[built_coord].size(); }

    inline typename Points_list::iterator begin() {
      return p_list[built_coord].begin();
    }

    inline typename Points_list::iterator end() {
      return p_list[built_coord].end();
    }

    // building the container from a sequence of points
    template <class Iter>
    Points_container(const int d, Iter begin, Iter end) :
      p_list(new Points_list[d]), bbox(d), tbox(d) {

      //      std::for_each(begin, end, points.push_back);
      //      for (; begin != end; ++begin) points.push_back(*begin);

      //      bbox = Box<NT>(d, points.begin(), points.end());
      //      tbox = bbox;
      //      std::for_each(points.begin(), points.end(),
      //		    build_max_span_list(p_list + max_span_coord()));

      bbox = Box<NT>(d, begin, end);
      tbox = bbox;
      std::for_each(begin, end,
		    build_max_span_list(p_list + max_span_coord()));


      p_list[max_span_coord()].sort(comp_coord_val<Item>(max_span_coord()));
      built_coord = max_span_coord();
    }

	// building an empty container 
//    template <class Iter>
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
	  assert(built_coord=c.built_coord);
	  merge(p_list[built_coord], c.p_list[built_coord], Less_lexicographically_d());
	  /* alternative implemenetation
	  bool something_done=true;
	  for (int i = 0; i < dimension(); ++i) {
	    // add  c.p_list[i] to p_list[i]
		  if (!(p_list[i].empty()) && (!(c.p_list[i].empty))) {
			p_list[i].splice(p_list[i].end(), c.p_list(i));
			something_done=true;
		  }
	  }
	  assert(something_done);*/

	}

    void recompute_tight_bounding_box() {
		tbox.update_from_point_pointers(p_list[built_coord].begin(),
		     p_list[built_coord].end(),p_list[built_coord].empty());
	}


    
      template <class Separator>
      void split_container(Points_container<Item>& c, Separator* sep, bool sliding=false) {

		assert(dimension()==c.dimension());
		

        c.bbox=bbox;
        // bool test_validity=false;

        const int split_coord = sep->cutting_dimension();
        const NT cutting_value = sep->cutting_value();

		/*
		std::cout << "container size =" << size() << std::endl;
		std::cout << "split_coord=" << split_coord << std::endl;
		std::cout << "cutting_value=" << cutting_value << std::endl;
		std::cout << "dimension=" << dimension() << std::endl;
        */

		// prepare the coordinate, if necessary
		if (p_list[split_coord].empty()) {
			// copy p_list[built_coord] to p_list[split_coord]
			p_list[split_coord]=p_list[built_coord];
			// sort p_list[split_coord]
			p_list[split_coord].sort(comp_coord_val<Item>(split_coord));
		}
        
		//splitting builds list along the split_coord
		built_coord=split_coord;
		c.built_coord=split_coord;

		// splitting the lists in two; can be done better for
		// i == split_coord...

		Points_list tmp_list(0);
		for (int i = 0; i < dimension(); ++i) {
			if (! p_list[i].empty()) {
                NT min_val=bbox.upper(split_coord); //init with upperbound
				NT max_val=bbox.lower(split_coord);  //init with lowerbound;
                Points_list::iterator pt_min=p_list[i].begin();
                Points_list::iterator pt_max=p_list[i].begin();
				tmp_list.clear();
				c.p_list[i].clear();
				
				for (typename Points_list::iterator pt = p_list[i].begin();
				pt != p_list[i].end(); pt= p_list[i].begin()) {
					if ((*(*pt))[split_coord] < min_val) {
						pt_min=pt; min_val= (*(*pt))[split_coord];
                    }
					if ((*(*pt))[split_coord] > max_val) {
						pt_max=pt; max_val= (*(*pt))[split_coord];
                    }
					if (sep->side(*(*pt)) == ON_NEGATIVE_SIDE) {
						c.p_list[i].splice(c.p_list[i].end(), p_list[i], pt); 
					}
					else {
						tmp_list.splice(tmp_list.end(), p_list[i], pt); 
					}
				}
				// in-place copy tmp_list to p_list[i]
				p_list[i].splice(p_list[i].end(), tmp_list);
				
				if (sliding) { // avoid empty list
					if (p_list[i].empty()) {// move maximal value to p_list
						p_list[i].splice(p_list[i].end(), c.p_list[i], pt_max);
					} 
					if (c.p_list[i].empty()) {// move minimal value to c.p_list
						c.p_list[i].splice(c.p_list[i].end(), p_list[i], pt_min);
					} 
				}
		}
		}
		
		


        /*
		if (sliding) { // then each list should contain at least one element
			if (p_list[split_coord].empty()) { // move last element of c.p_list to p_list
				// std::cout << "moved last element of c.p_list to p_list" << std::endl;
                                // test_validity=true;
				Points_list::const_reference Back_c_p_list=c.p_list[split_coord].back();
				for (int i = 0; i < dimension(); ++i) {
					if (! c.p_list[i].empty()) {
						p_list[i].push_front(Back_c_p_list);
						c.p_list[i].remove(Back_c_p_list);
					}
				}
			};

			{ for (int i = 0; i < dimension(); ++i)
	    if (! p_list[i].empty()) assert(p_list[i].size() == size());
        }

		{ for (int i = 0; i < dimension(); ++i)
	    if (! c.p_list[i].empty()) assert(c.p_list[i].size() == c.size());
        }

			if (c.p_list[split_coord].empty()) { // move first element of p_list to c.p_list
				// std::cout << "moved first element of p_list to c.p_list" << std::endl;
                                // test_validity=true;
				Points_list::const_reference Front_p_list=p_list[split_coord].front();
				for (int i = 0; i < dimension(); ++i) {
					if (! p_list[i].empty()) {
						c.p_list[i].push_back(Front_p_list);
						std::cout << "i in loop = " << i << std::endl;
						std::cout << "p_list size before remove=" << p_list[i].size() << std::endl;
						p_list[i].remove(Front_p_list);
						std::cout << "p_list size after remove=" << p_list[i].size() << std::endl;
					}
				}
			}

			{ for (int i = 0; i < dimension(); ++i)
	    if (! p_list[i].empty()) 
			if (! (p_list[i].size() == size())) {
				
				
				std::cout << "split_coord=" << split_coord << std::endl;
				std::cout << "size=" << size() << std::endl;
				std::cout << "i=" << i << std::endl;
				std::cout << "p_list size=" << p_list[i].size() << std::endl;
				std::cout << "points:" << std::endl;
				for (typename Points_list::iterator p = p_list[i].begin();
					p != p_list[i].end(); ++p) std::cout << **p << " "; std::cout << std::endl;

			}
        } */

		{ for (int i = 0; i < dimension(); ++i)
	    if (! c.p_list[i].empty()) assert(c.p_list[i].size() == c.size());
        }
		// }

		// Alternatively split list only in two for i == split_coord
		// adjusting boxes
		bbox.set_lower(split_coord, cutting_value);
		tbox.update_from_point_pointers(p_list[built_coord].begin(),
		p_list[built_coord].end(),p_list[built_coord].empty());
		c.bbox.set_upper(split_coord, cutting_value);
		c.tbox.update_from_point_pointers(c.p_list[c.built_coord].begin(),
		c.p_list[c.built_coord].end(),c.p_list[c.built_coord].empty());
        if (true) {is_valid(); c.is_valid();};
	}

	NT median(const int split_coord) {
      if (p_list[split_coord].empty()) {
         // copy p_list[built_coord] to p_list[split_coord]
         p_list[split_coord]=p_list[built_coord];
         // sort p_list[split_coord]
         p_list[split_coord].sort(comp_coord_val<Item>(split_coord));
      }
      Points_list::iterator median_point_ptr=p_list[split_coord].begin();
      for (unsigned int i = 0; i < p_list[split_coord].size()/2-1; i++, median_point_ptr++) {}
      NT val1=(*(*median_point_ptr))[split_coord];
      median_point_ptr++;
      NT val2=(*(*median_point_ptr))[split_coord];
      return (val1+val2)/2;
    }


    ~Points_container() { delete [] p_list; }

    inline bool empty() { return size() == 0;}

    bool is_valid() {return true;}

    /*
    bool is_valid() {
      // checking that all the lists are of the same size
     { for (int i = 0; i < dimension(); ++i)
	if (! p_list[i].empty()) assert(p_list[i].size() == size());
     }
      // checking that all the points lie in the box specified
     {

	  // std::cout << "warning: this call of Poins_cointaner::is_valid() takes much computing time" << std::endl;
      // std::cout << bbox; // bbox.print(std::cout);


      // extra test added by Hans to see if all pointer lists are sorted

      for (int i = 0; i < dimension(); ++i) {
      if (! p_list[i].empty())   {
        typename Points_list::iterator p1=p_list[i].begin();
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
     }
     // std::cout << std::endl;
      // checking that all the lists point to the same Points (debug-only...)
     {
      //  define RT for lexicographically_smaller
      // typename RT::Less_lexicographically_d lt;  error
      // lexicographically_smaller<RT> Compare;
      // std::set<Item, Less_lexicographically_d>  t,t1;


      // test did not run using lexicographically smaller order

      // std::set<Item> t,t1;
      /*
      typedef typename Item::R  RT;
      std::set <Item, lexicographically_smaller<RT>() >;
      for (typename Points_list::iterator p = p_list[built_coord].begin();
	   p != p_list[built_coord].end(); ++p)
	   t.insert(*(*p));
      assert(t.size() == size());

      for (int i = 0; i < dimension(); ++i)
	if (! p_list[i].empty()) {
	  t1.clear();
	  for (typename Points_list::iterator p = p_list[i].begin();
	       p != p_list[i].end(); ++p) t1.insert(*(*p));
	  assert(t == t1);
	}
     
     }
      return true;
    }*/

      // friend std::ostream& operator<< CGAL_NULL_TMPL_ARGS
      // (std::ostream&, Points_container<P>&);
    // Self&);


  private:
    explicit Points_container() {} // disable default constructor
  
  };

} // namespace CGAL

#endif // CGAL_POINTS_CONTAINER_H


