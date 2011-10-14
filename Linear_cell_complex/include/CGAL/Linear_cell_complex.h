// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_COMBINATORIAL_MAP_WITH_POINTS_H
#define CGAL_COMBINATORIAL_MAP_WITH_POINTS_H 1

#include <CGAL/Combinatorial_map.h>
#include <CGAL/Linear_cell_complex_min_items.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/Direction_2.h>
#include <CGAL/predicates_d.h>

namespace CGAL {

  template <unsigned int d>
  struct Default_template_argument
  { typedef Linear_cell_complex_traits<d,Cartesian_d<double> > type; };
  template <>
  struct Default_template_argument<2>
  { typedef Linear_cell_complex_traits<2,Exact_predicates_inexact_constructions_kernel> type; };
  template <>
  struct Default_template_argument<3>
  { typedef Linear_cell_complex_traits<3,Exact_predicates_inexact_constructions_kernel> type; };
	
  /** @file Combinatorial_map_with_points.h
   * Definition of combinatorial map with points associated to all vertices.
   */

  /** Combinatorial map with point associated to each vertex.
   * The Combinatorial_map_with_points class describes a
   * combinatorial map with point associated to each vertex.
   */
  template < unsigned int d_, unsigned int ambiant_dim = d_,
						 class Traits_ = typename Default_template_argument<d_>::type,
						 class Items_ = Combinatorial_map_with_points_min_items<d_>,
						 class Alloc_ = CGAL_ALLOCATOR(int),
						 template<unsigned int,class,class,class> class CMap =  Combinatorial_map_base >
  class Combinatorial_map_with_points:
    /*    public Combinatorial_map_base<d_, 
				  Combinatorial_map_with_points<d_, ambiant_dim,
								Items, Alloc>,
								Items, Alloc>,*/
    public CMap<d_,Combinatorial_map_with_points<d_, ambiant_dim, Traits_,
						 Items_, Alloc_>,Items_,Alloc_>,
    public Linear_cell_complex_traits<ambiant_dim, Traits_>
  {
  public:
    typedef Combinatorial_map_with_points<d_, ambiant_dim, Traits_, Items_, Alloc_>  Self;
    typedef Combinatorial_map_base<d_, Self, Items_, Alloc_>                Base;

		typedef Traits_ Traits;
		typedef Items_ Items;
		typedef Alloc_ Alloc;
		
		
    static const unsigned int ambiant_dimension = ambiant_dim;

    typedef typename Base::Dart_handle                    Dart_handle;
    typedef typename Base::Dart_const_handle              Dart_const_handle;

    typedef typename Base::Helper             Helper;

		//    typedef typename Items::Traits Traits;
		//    typedef basic_types<ambiant_dimension, Traits> Traits_types;
    
    // typedef typename Traits::Kernel Kernel;
		//    typedef typename basic_types<ambiant_dimension, Traits>::Point  Point;
		//    typedef typename basic_types<ambiant_dimension, Traits>::Vector Vector;
		typedef typename Traits::Point  Point;
		typedef typename Traits::Vector Vector;
    typedef typename Traits::FT     FT;

    typedef typename Base::Dart_range             Dart_range;

    typedef typename Helper::template Attribute_type<0>::type 
                Vertex_attribute;
    typedef typename Helper::template Attribute_handle<0>::type 
                Vertex_attribute_handle;
    typedef typename Helper::template Attribute_const_handle<0>::type 
                Vertex_attribute_const_handle;

    typedef typename Helper::template Attribute_range<0>::type 
          Vertex_attribute_range;

    typedef typename Helper::template Attribute_const_range<0>::type 
          Vertex_attribute_const_range;

    typedef typename Base::size_type size_type;

    /// To use previous definition of create_dart methods.
    using Base::create_dart;

    /// Reserve vertices and darts in compact containers
    void reserve( size_type vertices, size_type darts )
    {
      /*      this->mdarts.reserve(darts);
      CGAL::cpp0x::get<Helper::template Dimension_index<0>::value>
      (this->mattribute_containers).reserve(vertices);*/
    }

    /**
     * Create a new dart associated with an handle through an attribute.
     * @param ahandle the point handle to associated with the dart.
     * @return a Dart_handle on the new dart.
     */
    Dart_handle create_dart(Vertex_attribute_handle ahandle)
    {
      Dart_handle res = create_dart();
      this->template set_attribute_of_dart<0>(res,ahandle);
      return res;
    }

    /** Create a new dart associated with a point.
     * @param apoint the point to associated with the dart.
     * @return a Dart_handle on the new dart.
     */
    Dart_handle create_dart(const Point& apoint)
    { return create_dart(create_vertex_attribute(apoint)); }

    /** Create a vertex attribute.
     * @return an handle on the new attribute.
     */
    Vertex_attribute_handle create_vertex_attribute()
    { return Base::template create_attribute<0>(); }

    /** Create a vertex attribute associated with a point.
     * @param point the point to associated with the dart.
     * @return an handle on the new attribute.
     */
    Vertex_attribute_handle create_vertex_attribute(const Point& apoint)
    { return Base::template create_attribute<0>(apoint); }

    /** Create a vertex attribute associated with a point.
     * @param point the point to associated with the dart.
     * @return an handle on the new attribute.
     */
    void erase_vertex_attribute(Vertex_attribute_handle ahandle)
    { Base::template erase_attribute<0>(ahandle); }

    /// @return the Vertex_attribute_range for all vertex_attributes.
    Vertex_attribute_range& vertex_attributes()
    { return this->template attributes<0>(); }

    /// @return the Vertex_attribute_const_range for all vertex_attributes.
    Vertex_attribute_const_range& vertex_attributes() const
    { return this->template attributes<0>(); }

    /// @return the size of the vertex_attribute container.
    typename Base::size_type number_of_vertex_attributes() const
    { return Base::template number_of_attributes<0>(); }

    /** Set the vertex attribute of the given dart.
     * @param adart a dart.
     * @param ah the attribute to set.
     */
    void set_vertex_attribute_of_dart(Dart_handle adart, 
				      Vertex_attribute_handle ah)
    { return Base::template set_attribute_of_dart<0>(adart,ah); }

    /** Set the vertex attribute of all the darts of the vertex.
     * @param adart a dart of the vertex.
     * @param ah the attribute to set.
     */
    void set_vertex_attribute(Dart_handle adart, 
			      Vertex_attribute_handle ah)
    { return Base::template set_attribute<0>(adart,ah); }

    /// Get the vertex_attribute associated with a dart.
    /// @param a dart
    /// @return the vertex_attribute.
    static Vertex_attribute_handle vertex_attribute(Dart_handle adart)
    { 
      CGAL_assertion(adart!=NULL);
      return adart->template attribute<0>();
    }

    /// Get the vertex_attribute associated with a const dart.
    /// @param a dart
    /// @return the vertex_const_attribute.
    static Vertex_attribute_const_handle vertex_attribute(Dart_const_handle adart)
    { 
      CGAL_assertion(adart!=NULL);
      return adart->template attribute<0>();
    }

    /// Get the point associated with a dart.
    /// @param a dart
    /// @return the point.
    static Point& point(Dart_handle adart)
    { 
      CGAL_assertion(adart!=NULL && adart->template attribute<0>()!=NULL );
      return *adart->template attribute<0>();
    }

    /// Get the point associated with a const dart.
    /// @param a dart
    /// @return the point.
    static const Point& point(Dart_const_handle adart)
    { 
      CGAL_assertion(adart!=NULL && adart->template attribute<0>()!=NULL );
      return *adart->template attribute<0>();
    }

    /** Test if the map is valid.
     * A Combinatorial_map_with_points is valid if it is a valid 
     * Combinatorial_map with an attribute associated to each dart.
     * @return true iff the map is valid.
     */
    bool is_valid() const
    {
      bool valid = Base::is_valid();
      for (typename Dart_range::const_iterator it(this->darts().begin()), 
	     itend(this->darts().end()); valid && it != itend; ++it)
      {
        if ( vertex_attribute(it) == NULL )
        {
          std::cerr << "Map not valid: dart "<<&(*it)
		    <<" does not have a vertex."<< std::endl;
          valid = false;
        }
      }
      return valid;
    }

    /** test if the two given facets have the same geometry
     * @return true iff the two facets have the same geometry.
     */
    bool are_facets_same_geometry(Dart_const_handle d1, 
				 Dart_const_handle d2) const
    {
      typename Base::template Dart_of_orbit_range<1>::const_iterator 
	it1(*this,d1);
      typename Base::template Dart_of_orbit_range<0>::const_iterator 
	it2(*this,d2);
      bool samegeometry = true;
      for ( ; samegeometry && it1.cont() && it2.cont(); ++it1, ++it2)
	{
	  if ( it2->other_extremity()!=NULL && 
	       point(it1)!=point(it2->other_extremity()) )
	    samegeometry = false;
	}
      if ( it1.cont() != it2.cont() ) samegeometry = false;
      return samegeometry;
    }

    /// Sew3 the facets having same geometry basic version in O(n^2)
    /*    unsigned int sew3_same_facets_basic()
    {
      unsigned int res = 0;
      for (typename Dart_range::iterator it(this->darts().begin()), 
	     itend(this->darts().end()); it!=itend; ++it )
	for (typename Dart_range::iterator it2(this->darts().begin());
	       it2!=itend; ++it2 )
	  {
	    if ( it!=it2 && it->is_free(3) && it2->is_free(3) && 
		 are_facets_same_geometry(it,it2) )
	      {
		++res;
		this->template sew<3>(it,it2);
	      }
	  }
        return res;
	};*/

    /// Sew3 the facets having same geometry - improved version with std::map
    unsigned int sew3_same_facets()
    {
      unsigned int res = 0;

      std::map<Point, std::vector<Dart_handle> > one_dart_per_facet;
      int mymark = this->get_new_mark();
      CGAL_assertion( mymark!=-1 );

      // First we fill the std::map by one dart per facet, and by using
      // the minimal point as index.
      for (typename Dart_range::iterator it(this->darts().begin()), 
	     itend(this->darts().end()); it!=itend; ++it )
	{
	  if ( !this->is_marked(it, mymark) )
	    {
	      Point min_point=*(it->template attribute<0>());
	      Dart_handle min_dart = it;
	      this->mark(it, mymark);
	      typename Base::template Dart_of_orbit_range<1>::iterator
		it2(*this,it);
	      ++it2;	      
	      for ( ; it2.cont(); ++it2 )
		{
		  Point cur_point=*(it2->template attribute<0>());
		  this->mark(it2, mymark);
		  if ( cur_point < min_point )
		    {
		      min_point = cur_point;
		      min_dart = it2;
		    }
		}
	      one_dart_per_facet[min_point].push_back(min_dart);
	    }
	}

      // Second we run through the map: candidates for sew3 have necessary the
      // same minimal point.
      typename std::map<Point, std::vector<Dart_handle> >::iterator
	itmap=one_dart_per_facet.begin(),
	itmapend=one_dart_per_facet.end();
      
      for ( ; itmap!=itmapend; ++itmap )
	{
	  for ( typename std::vector<Dart_handle>::iterator
		  it1=(itmap->second).begin(),
		  it1end=(itmap->second).end(); it1!=it1end; ++it1 )
	    {
	      typename std::vector<Dart_handle>::iterator it2=it1;
	      for ( ++it2; it2!= it1end; ++it2 )
		{
		  if ( *it1!=*it2 && (*it1)->is_free(3) &&
		       (*it2)->is_free(3) && 
		       are_facets_same_geometry(*it1,(*it2)->beta(0)) )
		    {
		      ++res;
		      this->template sew<3>(*it1,(*it2)->beta(0));
		    }
		}
	    }
	}

      CGAL_assertion( this->is_whole_map_marked(mymark) );
      this->free_mark(mymark);
      return res;
    }

		
  };

} // namespace CGAL

#endif // CGAL_COMBINATORIAL_MAP_WITH_POINTS_H //
// EOF //
