{ int i=0;
  typedef Triangulation::Full_cell_handle Full_cell_handle;
  typedef Triangulation::Facet Facet;
  typedef std::vector<Full_cell_handle> Full_cells;

  Full_cells infinite_full_cells;
  std::back_insert_iterator<Full_cells> out(infinite_full_cells);

  t.incident_full_cells(t.infinite_vertex(), out);

  for( Full_cells::iterator sit = infinite_full_cells.begin();
       sit != infinite_full_cells.end(); ++sit )
    {
      Facet ft(*sit, (*sit)->index(t.infinite_vertex()));
      ++i; // |ft| is a facet of the convex hull
    }
  std::cout << "There are " << i << " facets on the convex hull."<< std::endl;
}
