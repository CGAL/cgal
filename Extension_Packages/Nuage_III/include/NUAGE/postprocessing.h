//=====================================================================
#ifndef NUAGE_POSTPROCESSING_H
#define NUAGE_POSTPROCESSING_H
//=====================================================================

//---------------------------------------------------------------------
// En principe, si l'allocateur de cellules etait bien fait on aurait pas besoin 
// de mettre a jour les valeurs rajoutees pour les cellules a  la main...

void re_init_for_free_cells_cache(const Triangulation_3& A, 
				  const Vertex_handle& vh)
{
  std::list<Cell_handle> ch_set;
  A.incident_cells(vh, std::back_inserter(ch_set));
  for (std::list<Cell_handle>::iterator c_it = ch_set.begin();
       c_it != ch_set.end(); 
       c_it++)
    (*c_it)->clear();
}

//---------------------------------------------------------------------

void swap_selected_facets_on_conflict_boundary(const Triangulation_3& A, 
					       const Vertex_handle& vh)
{
  std::list<Cell_handle> ch_set;
  A.incident_cells(vh, std::back_inserter(ch_set));
  for (std::list<Cell_handle>::iterator c_it = ch_set.begin();
       c_it != ch_set.end(); c_it++)
    {
      Cell_handle c = *c_it;
      int index = c->index(vh);
      Cell_handle neigh = c->neighbor(index);
      int n_ind = neigh->index(c);
      neigh->set_smallest_radius(n_ind, -1); // pour obliger le recalcul
      // si c est selectionnee c'est qu'elle est aussi le mem_IFacet renvoye par 
      // compute_value... donc a swapper aussi
      if (c->is_selected_facet(index))
      {
	c->unselect_facet(index);
	neigh->select_facet(n_ind);
	int i1 = (n_ind+1) & 3;
	int i2 = (n_ind+2) & 3;
	int i3 = (n_ind+3) & 3;
	Edge_like key(neigh->vertex(i1), neigh->vertex(i2));
	Border_elt result;
	if (is_border_elt(key, result))
	  {
	    Edge_IFacet ei_facet(void_Edge((void*) &*neigh, i1, i2), 
				 n_ind);
	    *get_border_IO_elt(key.first, key.second) =
	      IO_edge_type(ei_facet, ei_facet);
	  }
	key = Edge_like(neigh->vertex(i1), neigh->vertex(i3));
	if (is_border_elt(key, result))
	  {
	    Edge_IFacet ei_facet(void_Edge((void*) &*neigh, i1, i3), 
				 n_ind);
	    *get_border_IO_elt(key.first, key.second) =
	      IO_edge_type(ei_facet, ei_facet);
	  }
	key = Edge_like(neigh->vertex(i3), neigh->vertex(i2));
	if (is_border_elt(key, result))
	  {
	    Edge_IFacet ei_facet(void_Edge((void*) &*neigh, i3, i2), 
				 n_ind);
	    *get_border_IO_elt(key.first, key.second) =
	      IO_edge_type(ei_facet, ei_facet);
	  }
      }
    }
}

//---------------------------------------------------------------------

Facet get_next_selected_facet_around_edge(const Edge_IFacet& start)
{
  Edge_IFacet circ = inc_facet_circ(start);
  Cell_handle c = (Cell*) start.first.first;
  do
    {
      Cell_handle ch = (Cell*) circ.first.first;
      int ind = circ.second;
      Cell_handle neigh = ch->neighbor(ind);
      int n_ind = neigh->index(ch);
      if (ch->is_selected_facet(ind))
	return Facet(ch, ind);
      if (neigh->is_selected_facet(n_ind))
	return Facet(neigh, n_ind);
      circ = inc_facet_circ(circ);
    }
  while(Cell_handle ((Cell*) circ.first.first) != c);
  // si on passe par la, alors y a eu un probleme....
  std::cerr << "+++probleme dans la MAJ avant remove..." << std::endl;
  
  return Facet(c, start.second);
}

//---------------------------------------------------------------------

void retract_border_for_incident_facets(const Vertex_handle& vh)
{
  Next_border_elt border_elt =  *(vh->first_incident());
  int border_index = border_elt.second.second;
  Vertex_handle vh_succ = (Vertex*) border_elt.first;
  IO_edge_type io_edge = border_elt.second.first.second;
  Edge_IFacet i_facet = io_edge.first;
  Cell_handle c = (Cell*) i_facet.first.first;
  int i1 = c->index(vh);
  int i2 = c->index(vh_succ);
  int index = i_facet.second;
  int i3 = 6 - index - i1 - i2;
  Vertex_handle vh_int = c->vertex(i3);
  _ordered_map_erase(border_elt.second.first.first, 
		     get_border_IO_elt(vh, vh_succ));
  vh->remove_border_edge((void*) &*vh_succ);
  // 1- a virer au cas ou car vh va etre detruit
  vh_succ->remove_interior_edge((void*) &*vh);
  bool while_cond(true);
  do
    {
      _facet_number--;

      assert(c->is_selected_facet(index));
      c->unselect_facet(index);
 
//        if (!vh_succ->is_on_border())
// 	{
// 	  vh_succ->re_init();
// 	}
      Facet f32 = 
	get_next_selected_facet_around_edge(Edge_IFacet(void_Edge((void*) &*c, i3, i2), 
							index));

     if (!vh_int->is_on_border())
	{
	  vh_int->re_init(); 
	  vh_int->inc_mark();
// 	  std::list<Vertex_handle> vh_set;
// 	  A.incident_vertices(vh_int, std::back_inserter(vh_set));
// 	  for (std::list<Vertex_handle>::iterator v_it = vh_set.begin();
// 	       v_it != vh_set.end(); v_it++)
// 	    if((*v_it)->is_on_border())
// 	      {
// 		// pour retrouver cette info, on a besoin de savoir si l'arete
// 		// [vh_hint, *v_it] est une arete de la reconstruction...
// 		vh_int->set_interior_edge(*v_it);
// 	      }
	}

      Edge_IFacet e32(void_Edge((void*) &*f32.first, 
				f32.first->index(vh_int),
				f32.first->index(vh_succ)), f32.second);
      Radius_edge_type rad_elt_32(STANDBY_CANDIDATE, IO_edge_type(e32, e32)); 
      Border_elt result;
      if (is_ordered_border_elt(Edge_like(vh_int, vh), result))
	{
	  _ordered_map_erase(result.first.first, get_border_IO_elt(vh_int, vh));
	  vh_int->remove_border_edge((void*) &*vh);
	  // 1- a virer au cas ou car vh va etre detruit
	  vh_int->remove_interior_edge((void*) &*vh);
	  while_cond = false;
	}
      // a titre  preventif... on essaye de s'assurer de marquer les aretes
      // interieures au sens large...

      // 2- a virer a tout pris pour que maintenir le sens de interior edge
      vh_int->remove_interior_edge((void*) &*vh_succ);
      vh_succ->remove_interior_edge((void*) &*vh_int);
      
      IO_edge_type* p32 = set_border_elt(vh_int, vh_succ, 
					 Border_elt(rad_elt_32, border_index));
      _ordered_border->insert(Radius_ptr_type (STANDBY_CANDIDATE, p32));

      // incrementation...
      if (while_cond)
	{
	  Facet f31 = 
	    get_next_selected_facet_around_edge(Edge_IFacet(void_Edge((void*) &*c, i3, i1), 
							    index));

	  c = f31.first;
	  index = f31.second;
	  i1 = c->index(vh);
	  vh_succ = vh_int;
	  i2 = c->index(vh_int);
	  i3 = 6 - index - i1 - i2;
	  vh_int = c->vertex(i3);      
	}
    }
  while(while_cond);
}

//---------------------------------------------------------------------

inline bool create_singularity(const Triangulation_3& A, 
			       const Vertex_handle& vh)
{
  // Pour reperer le cas de triangle isole 
  if (vh->is_on_border())
    {
      // vh sommet 0
      Next_border_elt border_elt =  *(vh->first_incident());
      Vertex_handle vh_1 = (Vertex*) border_elt.first;// sommet 1
      border_elt =  *(vh_1->first_incident());
      Vertex_handle vh_2 = (Vertex*) border_elt.first;// sommet 2
      border_elt =  *(vh_2->first_incident());
      Vertex_handle vh_3 = (Vertex*) border_elt.first;// sommet 0 ???
      Cell_handle c;
      int i, j, k;
      if ((vh_3 == vh)&&(A.is_facet(vh, vh_1, vh_2, c, i ,j ,k)))
	{
	  int l = 6-i-j-k;
	  Cell_handle neigh = c->neighbor(l);
	  
	  if
	    (c->is_selected_facet(l)||neigh->is_selected_facet(neigh->index(c)))
	    return true;
	}
    }
  

  // Reperer le cas d'aretes interieures...
  std::list<Vertex_handle> vh_list;
  A.incident_vertices(vh, back_inserter(vh_list));

  for (std::list<Vertex_handle>::iterator v_it = vh_list.begin();
       v_it != vh_list.end(); v_it++)
    if ((*v_it)->is_on_border() && is_interior_edge(Edge_like(vh, *v_it)))
      return true;
  return false;
}


//---------------------------------------------------------------------
 
struct Remove : public std::unary_function<Vertex_handle, bool>
{
  Triangulation_3& T;

  Remove(Triangulation_3& T_) : T(T_) {}

  bool operator()(Vertex_handle vh) {
    if (vh->is_exterior())
      { 
	swap_selected_facets_on_conflict_boundary(T, vh);
	re_init_for_free_cells_cache(T, vh);
	if (!T.remove(vh))
	  std::cerr << "+++Delaunay_triangulation_3.remove(Vertex_handle) failed."  <<
	    convert()(vh->point()) << std::endl;
	return true;
      }
    else if (vh->is_on_border()&&(!create_singularity(T, vh)))
      {      
	swap_selected_facets_on_conflict_boundary(T, vh);
	retract_border_for_incident_facets(vh);
	re_init_for_free_cells_cache(T, vh);
	_vh_number--;
	if (!T.remove(vh))
	  std::cerr << "+++Delaunay_triangulation_3.remove(Vertex_handle) failed." <<
	    convert()(vh->point()) << std::endl;
	return true;
      }
    else
      { }
    return false;
  }
};


//---------------------------------------------------------------------

bool postprocessing(Triangulation_3& A, const int& NB_BORDER_MAX)
{  
  _postprocessing_cont++;
  t1.start();
  std::list<Vertex_handle> L_v;

  // Pour prendre en compte tous sommets exterieurs ou sur le bord
//   for(Finite_vertices_iterator v_it = A.finite_vertices_begin();
//       v_it != A.finite_vertices_end(); v_it++)
//     {
//       if (v_it->number_of_incident_border() != 0)
// 	{
// 	  L_v.push_back(v_it->handle());
// 	  v_it->erase_incidence_request();
// 	}
//     }

  //  Pour controler les sommets choisis sur le bord...
  
  // nombre d'aretes a partir duquel on considere que c'est irrecupperable NB_BORDER_MAX

  int vh_on_border_inserted(0);
  for(Finite_vertices_iterator v_it = A.finite_vertices_begin();
      v_it != A.finite_vertices_end(); 
      v_it++)
    {
      v_it->erase_incidence_request();
      if ((v_it->number_of_incident_border() > 0)&&
	  (!v_it->is_post_marked(_postprocessing_cont)))
	{
	  std::list<Vertex_handle> L_v_tmp;
	  Vertex_handle vprev_it, vh_it;
// 	  Vertex_handle vsucc_it;
	  int v_count(0);
	  vprev_it = v_it->handle();
	  do
	    {		      
	      vh_it = (Vertex*) vprev_it->first_incident()->first;
// 	      vsucc_it = (Vertex*) vh_it->first_incident()->first;
// 	      D_Point p1 = convert()(vprev_it->point());
// 	      D_Point p = convert()(vh_it->point());
// 	      D_Point p2 = convert()(vsucc_it->point());
	      // pour imposer une condition sur l'angle d'aretes...
// 	      if ((p1-p)*(p2-p) > 0)
		L_v_tmp.push_back(vh_it);
	      vh_it->set_post_mark(_postprocessing_cont);
	      vprev_it = vh_it;
	      v_count++;
	    }
	  while((vprev_it != v_it->handle())&&(v_count < NB_BORDER_MAX));

	  if (v_count < NB_BORDER_MAX)
	    {
	      L_v.insert(L_v.begin(), L_v_tmp.begin(), L_v_tmp.end());
	      vh_on_border_inserted += v_count;
	    }

	} 
      if (v_it->number_of_incident_border() < 0)
	L_v.push_back(v_it->handle());
    }

  unsigned int itmp, L_v_size_mem;
  L_v_size_mem = L_v.size();
  if ((vh_on_border_inserted != 0)&& // pour ne post-traiter que les bords
      (L_v.size() < .1*_A_size_before_postprocessing))
    {
      {
	//CGAL::Protect_FPU_rounding<true> P;
	do
	  {
	    itmp = L_v.size();
	    std::list<Vertex_handle>::iterator new_end =
	      std::remove_if(L_v.begin(), L_v.end(), Remove(A));
	    L_v.erase(new_end, L_v.end());
	  }
	while (!L_v.empty() && (L_v.size() < itmp));
      }
      std::cout << "   " << L_v.size() << " points non reguliers." << std::endl;
      re_compute_values(A);
    }
  else
    return false;
  // +10% de points retires ou trop d'etapes (>20) -> stop
  if ((L_v_size_mem == L_v.size())||
      ((_A_size_before_postprocessing-A.number_of_vertices()) >
       .1*_A_size_before_postprocessing)||
      (_postprocessing_cont > 20))
    return false;

  min_K = HUGE_VAL;
  t1.stop();
  // fin--
//   if (_postprocessing_cont < 5)
//     return true;
  return true;
}

//=====================================================================
#endif //NUAGE_POSTPROCESSING_H
//=====================================================================
