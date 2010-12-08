//build our internal TDS using the TDS of an already built triangulation
//returns the number of vertices
template <class TDS_src,class TDS_tgt>
void
copy_tds(const TDS_src& src,TDS_tgt& tgt,typename TDS_src::Vertex_handle s_infinite,typename TDS_tgt::Vertex_handle t_infinite)
{
  int n = src.number_of_vertices();
  if (n == 0)  return; 
  tgt.cells().clear();
  tgt.set_dimension(src.dimension());

  // Number of pointers to cell/vertex to copy per cell.
  int dim = (std::max)(1, tgt.dimension() + 1);

  // Create the vertices.
  std::vector<typename TDS_src::Vertex_handle> TV(n);
  int i = 0;

  for (typename TDS_src::Vertex_iterator vit = src.vertices_begin();
       vit != src.vertices_end(); ++vit)
    TV[i++] = vit; 
  
  CGAL_triangulation_assertion( i == n ); 

  std::map< typename TDS_src::Vertex_handle, typename TDS_tgt::Vertex_handle > V;
  std::map< typename TDS_src::Cell_handle, typename TDS_tgt::Cell_handle > F;
  
  assert(TV[0]==s_infinite);
  assert(tgt.vertices_begin()==t_infinite);
  V[ TV[0] ] = t_infinite;
  for (i=1; i <= n-1; ++i){
    typename TDS_tgt::Vertex_handle vh=
      tgt.create_vertex( typename TDS_tgt::Vertex(TV[i]->point()) );
    V[ TV[i] ] = vh;
  }

  // Create the cells.
  for (typename TDS_src::Cell_iterator cit = src.cells_begin();
	  cit != src.cells_end(); ++cit) {
      //WE ARE LOOSING ALL INFO INSIDE CELL (HIDDEN POINT ETC...)
      F[cit] = tgt.create_cell();
      for (int j = 0; j < dim; j++)
        F[cit]->set_vertex(j, V[cit->vertex(j)] );
  }

  // Link the vertices to a cell.
  for (typename TDS_src::Vertex_iterator vit2 = src.vertices_begin();
       vit2 != src.vertices_end(); ++vit2)
    V[vit2]->set_cell( F[vit2->cell()] );

  // Hook neighbor pointers of the cells.
  for (typename TDS_src::Cell_iterator cit2 = src.cells_begin();
	  cit2 != src.cells_end(); ++cit2) {
    for (int j = 0; j < dim; j++)
      F[cit2]->set_neighbor(j, F[cit2->neighbor(j)] );
  }

  CGAL_triangulation_postcondition( tgt.is_valid() );
}
