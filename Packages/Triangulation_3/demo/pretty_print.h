#include <CGAL/basic.h>
#include <stdio.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <strstream.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Point_3.h>

#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_ds_iterators_3.h>

#include <CGAL/Triangulation_iterators_3.h>
#include <CGAL/Triangulation_circulators_3.h>
#include <CGAL/Triangulation_3.h>

#include <CGAL/IO/Geomview_stream.h>

typedef double coord_type;
typedef CGAL_Cartesian<coord_type>  Rep;

typedef CGAL_Point_3<Rep>  Point;

typedef CGAL_Triangulation_geom_traits_3<Rep> Gt;
typedef CGAL_Triangulation_vertex_base_3<Gt> Vb;
typedef CGAL_Triangulation_cell_base_3<Gt>  Cb;

typedef CGAL_Triangulation_data_structure_3<Vb,Cb> TDS;

typedef CGAL_Triangulation_ds_vertex_iterator_3<TDS> TDSVertex_iterator;
typedef CGAL_Triangulation_ds_edge_iterator_3<TDS> TDSEdge_iterator;
typedef CGAL_Triangulation_ds_facet_iterator_3<TDS> TDSFacet_iterator;
typedef CGAL_Triangulation_ds_cell_iterator_3<TDS> TDSCell_iterator;

typedef CGAL_Triangulation_vertex_iterator_3<Gt,TDS> Vertex_iterator;
typedef CGAL_Triangulation_edge_iterator_3<Gt,TDS> Edge_iterator;
typedef CGAL_Triangulation_facet_iterator_3<Gt,TDS> Facet_iterator;
typedef CGAL_Triangulation_cell_iterator_3<Gt,TDS> Cell_iterator;

typedef typename TDS::Cell TDSCell;
typedef typename TDS::Facet TDSFacet;
typedef typename TDS::Edge TDSEdge;
typedef typename TDS::Vertex TDSVertex;

typedef CGAL_Triangulation_data_structure_3<Vb,Cb> TDS;

typedef CGAL_Triangulation_3<Gt,TDS> Triangulation_3;

typedef typename Triangulation_3::Cell Cell;
typedef typename Triangulation_3::Facet Facet;
typedef typename Triangulation_3::Edge Edge;
typedef typename Triangulation_3::Vertex Vertex;

typedef typename Triangulation_3::Vertex_handle Vertex_handle;
typedef typename Triangulation_3::Cell_handle Cell_handle;

void pp_tds_vertex(const TDSVertex*);
void pp_tds_vertex(const TDSCell* c, int i);
void pp_tds_edge(const TDSEdge);
void pp_tds_edge(const TDSCell* c, int i, int j);
void pp_tds_facet(const TDSFacet);
void pp_tds_facet(const TDSCell* c, int i);
void pp_tds_cell(const TDSCell*);

void pp_point(const Point &);

void pp_vertex(const Vertex_handle);
void pp_edge(const Edge);
void pp_edge(const Cell_handle, int, int);
void pp_facet(const Facet);
void pp_facet(const Cell_handle, int);
void pp_cell(const Cell_handle);

void visu_face(CGAL_Geomview_stream & os, Triangulation_3 & T, Cell_handle c, int i);
