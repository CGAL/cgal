#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_3_H

#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/spatial_sort.h>

// Needed by remove to fill the hole.
#include <CGAL/Periodic_3_triangulation_remove_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>


namespace CGAL
{
template < class Gt,
		   class Tds = Triangulation_data_structure_3 < Triangulation_vertex_base_3<Gt, Periodic_3_triangulation_ds_vertex_base_3<> >,
		   	   	   	   Triangulation_cell_base_3<Gt, Periodic_3_triangulation_ds_cell_base_3<> > >
		 >
class Periodic_3_Regular_triangulation_3 : public Periodic_3_triangulation_3<Gt,Tds>
{
	typedef Periodic_3_Regular_triangulation_3<Gt,Tds>  Self;

public:
	typedef Periodic_3_triangulation_3<Gt,Tds>  Base;

	typedef Gt                                    Geometric_traits;
	typedef Tds                                   Triangulation_data_structure;

	typedef Geometric_traits                      Geom_traits;
	typedef typename Gt::FT                       FT;

	typedef typename Gt::Weighted_point_3         Weighted_point;
	typedef typename Gt::Bare_point               Bare_point;
	typedef typename Gt::Segment_3                Segment;
	typedef typename Gt::Triangle_3               Triangle;
	typedef typename Gt::Tetrahedron_3            Tetrahedron;

	typedef typename Base::Periodic_point         Periodic_point;
	typedef typename Base::Periodic_segment       Periodic_segment;
	typedef typename Base::Periodic_triangle      Periodic_triangle;
	typedef typename Base::Periodic_tetrahedron   Periodic_tetrahedron;

	typedef typename Base::Cell_handle            Cell_handle;
	typedef typename Base::Vertex_handle          Vertex_handle;

	typedef typename Base::Cell                   Cell;
	typedef typename Base::Vertex                 Vertex;
	typedef typename Base::Facet                  Facet;
	typedef typename Base::Edge                   Edge;

	typedef typename Base::Cell_circulator        Cell_circulator;
	typedef typename Base::Facet_circulator       Facet_circulator;
	typedef typename Base::Cell_iterator          Cell_iterator;
	typedef typename Base::Facet_iterator         Facet_iterator;
	typedef typename Base::Edge_iterator          Edge_iterator;
	typedef typename Base::Vertex_iterator        Vertex_iterator;

	typedef typename Base::All_cells_iterator     All_cells_iterator;
	typedef typename Base::All_facets_iterator    All_facets_iterator;
	typedef typename Base::All_edges_iterator     All_edges_iterator;
	typedef typename Base::All_vertices_iterator  All_vertices_iterator;

	typedef typename Base::size_type              size_type;
	typedef typename Base::difference_type        difference_type;

	typedef typename Base::Locate_type            Locate_type;
	typedef typename Base::Iterator_type          Iterator_type;

	typedef typename Base::Offset                 Offset;
	typedef typename Base::Iso_cuboid             Iso_cuboid;
	typedef typename Base::Covering_sheets        Covering_sheets;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
using Base::cw;
using Base::ccw;
using Base::domain;
using Base::geom_traits;
using Base::int_to_off;
using Base::number_of_sheets;
using Base::number_of_vertices;
using Base::number_of_edges;
using Base::number_of_facets;
using Base::number_of_cells;
using Base::cells_begin;
using Base::cells_end;
using Base::vertices_begin;
using Base::vertices_end;
using Base::facets_begin;
using Base::facets_end;
using Base::tds;
using Base::next_around_edge;
using Base::vertex_triple_index;
using Base::mirror_vertex;
using Base::orientation;
using Base::insert_dummy_points;
using Base::swap;
using Base::is_1_cover;
using Base::is_virtual;
using Base::point;
#endif

// For strict-ansi compliance
using Base::adjacent_vertices;
using Base::combine_offsets;
using Base::get_offset;
using Base::get_original_vertex;
using Base::get_location_offset;
using Base::get_neighbor_offset;
using Base::incident_edges;
using Base::incident_facets;
using Base::incident_cells;
using Base::is_valid_conflict;
using Base::locate;
using Base::periodic_point;
using Base::segment;

public:
	/** @name Creation */ //@{
	Periodic_3_Regular_triangulation_3 (const Iso_cuboid& domain = Iso_cuboid(0, 0, 0, 1, 1, 1),
			                            const Geometric_traits& gt = Geometric_traits())
	: Base(domain, gt)
	{
	}

protected:
	bool less_power_distance (const Bare_point &p, const Weighted_point &q, const Weighted_point &r)  const
	{
		return geom_traits().compare_power_distance_3_object()(p, q, r) == SMALLER;
	}

	bool less_power_distance (const Bare_point &p, const Weighted_point &q, const Weighted_point &r,
			                  const Offset &o1, const Offset &o2, const Offset &o3)  const
	{
		return geom_traits().compare_power_distance_3_object()(p, q, r, o1, o2, o3) == SMALLER;
	}

	Bare_point construct_weighted_circumcenter (const Weighted_point &p, const Weighted_point &q, const Weighted_point &r, const Weighted_point &s) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s);
	}

	Bare_point construct_weighted_circumcenter (const Weighted_point &p, const Weighted_point &q, const Weighted_point &r, const Weighted_point &s,
			                                    const Offset& o1, const Offset& o2, const Offset& o3, const Offset& o4) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s, o1,o2,o3,o4);
	}

	Bare_point construct_weighted_circumcenter(const Weighted_point &p, const Weighted_point &q, const Weighted_point &r) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r);
	}

	Bare_point construct_weighted_circumcenter(const Weighted_point &p, const Weighted_point &q, const Weighted_point &r,
			                                   const Offset& o1, const Offset& o2, const Offset& o3) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r);
	}
};
}// namespace CGAL

#endif
