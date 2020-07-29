{ int i=0;
  typedef Triangulation::Full_cell_iterator Full_cell_iterator;
  typedef Triangulation::Facet Facet;

  for( Full_cell_iterator cit = t.full_cells_begin();
       cit != t.full_cells_end(); ++cit )
    {
      if( ! t.is_infinite(cit) )
        continue;
      Facet ft(cit, cit->index(t.infinite_vertex()));
      ++i;// |ft| is a facet of the convex hull
    }
  std::cout << "There are " << i << " facets on the convex hull."<< std::endl;
}
