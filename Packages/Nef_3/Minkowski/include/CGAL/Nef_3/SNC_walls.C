

class SNC_walls {
  
  Vector_3 dir;

  Unique_hash_map<Vertex_handle, Vertex_handle> opposite_up;
  Unique_hash_map<Vertex_handle, Vertex_handle> opposite_down;
  Unique_hash_map<Vertex_handle, Volume_handle> volume_up;
  Unique_hash_map<Vertex_handle, Volume_handle> volume_down;

  Vertex_handle create_opposite_vertex(Vertex_const_handle v, const Vector_3 vec) {
    Object_handle o = pl()->shoot(vec);
    Vertex_handle vh;
    Halfedge_handle eh;
    Halffacet_handle fh;
    if(assign(fh,o))
      else if(assign(eh,o))
	else if(assign(vh,o))
	  else CGAL_assertion_msg(false, "wrong handle");
    
  }

  void create_walls() {

    Vertex_iterator vi;
    CGAL_forall_vertices(vi, *this) {

    Vertex_iterator vi;
    CGAL_forall_vertices(vi, *this) {
      opposite_up[vi] = create_opposite_vertex(vi, Vector_3(0,0,1));
      opposite_down[vi] = create_opposite_vertex(vi, Vector_3(0,0,-1));
    }

    Halfedge_iterator ei;
    CGAL_forall_vertices(ei, *this) {
      create_opposite_halfedge_up(opposite_up[ei->source()], 
				  opposite_up[ei->twin()->source()]);
      create_opposite_halfedge3_down(opposite_down[ei->source()], 
				     opposite_down[ei->twin()->source()]);
    }
  }

};
