//=====================================================================
#ifndef IOFILE_MANIPULATOR_H
#define IOFILE_MANIPULATOR_H
//=====================================================================

#include <iomanip>
#include <CGAL/algorithm.h>
//#define AF_CGAL_CLIB_STD std
#define AF_CGAL_CLIB_STD
//=====================================================================
//======================= INPUT DATA ==================================
//=====================================================================

void
file_output(char* foutput, std::vector<Point>& L) // to debug...
{
  std::ofstream os(foutput, std::ios::out);

  if(os.fail())
    std::cerr << "+++unable to open file for output" << std::endl;
  else
    std::cerr << "+++file for output : " << foutput << std::endl;

  os.clear();

  CGAL::set_ascii_mode(os);

  for(std::vector<Point>::iterator L_it=L.begin();
      L_it!=L.end(); L_it++)
    os << *L_it;
}

//---------------------------------------------------------------------

bool
file_input(char* finput, const int& number_of_points, std::vector<Point>& L)
{
  std::ifstream is(finput, std::ios::in);

  if(is.fail())
    {
      std::cerr << "+++unable to open file for input" << std::endl;
      std::exit(0);
      return false;
    }
  else
    std::cout << ">> input from file : " << finput << std::endl;

// pour selectionner le mode de lecture souhaite...
//   is.setf(std::ifstream::scientific);
//   is.setf(std::ifstream::showpos);
//   is.setf(std::ifstream::uppercase);

  CGAL::set_ascii_mode(is);

  int n;
  is >> n;
  std::cout << "   reading " << n << " points" << std::endl;
  Point p;
  L.reserve(n);
  CGAL::copy_n(std::istream_iterator<Point>(is), n, std::back_inserter(L));

  std::cout << "   random shuffling" << std::endl;
  std::random_shuffle(L.begin(), L.end());

  if ( (number_of_points > 0 ) && (number_of_points < n ))
    {
      L.erase(L.begin()+number_of_points, L.begin()+n);

      std::cout << std::endl 
		<< "   and randomize a sub-sample of " << number_of_points 
		<< " points." <<
	std::endl << std::endl;
    }
  //file_output("random_sample.data", L);
  return true;
}

//---------------------------------------------------------------------

bool
section_file_input(char* finput, const int& number_of_points, std::vector<Point>& L)
{
  std::ifstream is(finput, std::ios::in);

  if(is.fail())
    {
      std::cerr << "+++unable to open file for input" << std::endl;
      std::exit(0);
      return false;
    }
  else
    std::cout << ">> input from section-file : " << finput << std::endl;

// pour selectionner le mode de lecture souhaite...
//   is.setf(std::ifstream::scientific);
//   is.setf(std::ifstream::showpos);
//   is.setf(std::ifstream::uppercase);

  CGAL::set_ascii_mode(is);
      
  int i(0), N, points_num(0);
  char c;
  is >> c; // c == "S"
  is >> N;
  coord_type h;
  int n;

  do{
    is >> c; // c == "v"
    is >> n;
    is >> c; // c == "z"
    is >> h;
    is >> c; // c == "{"
    
    Rep::Point_2 p;

    for(; n > 0; n--)
      {
	is >> c;
	if (c == '}') 
	  {
	    is >> c; // c == "{"
	    n++;
	  }
	else
	  {
	    is.putback(c);
	    is >> p;
	    points_num++;
	    L.push_back(Point (p.x(), p.y(), h));
	  }
      }

    is >> c; // c == "}"
    i++;
  } while (i < N); 

  std::cout << "   reading " << points_num << " points";

  std::random_shuffle(L.begin(), L.end());

  if ( (number_of_points > 0 ) && (number_of_points < points_num ))
    {
      L.erase(L.begin()+number_of_points, L.begin()+points_num);

      std::cout << std::endl 
		<< "   and randomize a sub-sample of " << number_of_points 
		<< " points." <<
	std::endl << std::endl;
    }
  //file_output("core-dump.data", L);
  return true;
}

//=====================================================================
//======================= DUMP RESULT =================================
//=====================================================================

void
dump_in_file_medit_selected_facets(char* foutput, const Triangulation_3& A)
{ 
  char foutput_points[100];
  char foutput_faces[100];
  AF_CGAL_CLIB_STD::strcpy(foutput_points, foutput);
  AF_CGAL_CLIB_STD::strcpy(foutput_faces, foutput);
  std::strcat(foutput_points, ".points");
  std::strcat(foutput_faces, ".faces");
  std::ofstream os_points(foutput_points, std::ios::out);
  std::ofstream os_faces(foutput_faces, std::ios::out);
  if((os_points.fail())||(os_faces.fail()))
    std::cerr << "+++unable to open file for output" << std::endl;
  else
    std::cout << ">> files for output : " << foutput_points 
	      << ", " << foutput_faces << std::endl;

  os_points.clear();
  os_faces.clear();

  CGAL::set_ascii_mode(os_points);
  CGAL::set_ascii_mode(os_faces);
  

  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if ((!v_it->is_exterior()) && (v_it->get_visu_index()==-1))
      {
	v_it->set_visu_index(_vh_number);
	_vh_number++;
      }

  int _vh_size(_vh_number);

  os_points << _vh_size << std::endl;

  typedef const Point*  Const_point_star;
  Const_point_star * points_tab = new  Const_point_star[_vh_size];
  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if (v_it->get_visu_index() > -1)
      points_tab[v_it->get_visu_index()] = &v_it->point();

  for(int vh_i=0; vh_i<_vh_size; vh_i++)
    os_points << convert()(*points_tab[vh_i]) << " 0" << std::endl;

  os_faces << _facet_number << std::endl;
  for(Finite_facets_iterator f_it = A.finite_facets_begin(); 
      f_it != A.finite_facets_end(); 
      f_it++)
    {
      Cell_handle n, c = (*f_it).first;
      int ni, ci = (*f_it).second;
      n = c->neighbor(ci);
      ni = n->index(c);
      int i1, i2 ,i3;

      if (c->is_selected_facet(ci))
	{
	  i1 = (ci+1) & 3;
	  i2 = (ci+2) & 3;
	  i3 = (ci+3) & 3;
	  os_faces << 3 << " ";
	  os_faces << c->vertex(i1)->get_visu_index()+1 << " ";
	  os_faces << c->vertex(i2)->get_visu_index()+1 << " ";
	  os_faces << c->vertex(i3)->get_visu_index()+1 << " ";
	  os_faces << " 0 0 0 0" << std::endl; 
	}
      
       if (n->is_selected_facet(ni))
	{
	  i1 = (ni+1) & 3;
	  i2 = (ni+2) & 3;
	  i3 = (ni+3) & 3;
	  os_faces << 3 << " ";
	  os_faces << n->vertex(i1)->get_visu_index()+1 << " ";
	  os_faces << n->vertex(i2)->get_visu_index()+1 << " ";
	  os_faces << n->vertex(i3)->get_visu_index()+1 << " ";
	  os_faces << " 0 0 0 0" << std::endl; 
	}
    }
  // iterer sur _additional_facets_list pour rajouter les facttes manquantes
  for(Additional_facets_list_iterator add_f_it =
	_additional_facets_list->begin();
      add_f_it != _additional_facets_list->end(); add_f_it++)
    {
      os_faces << 3 << " ";
      os_faces << (*add_f_it).first->get_visu_index()+1 << " ";
      os_faces << (*add_f_it).second->get_visu_index()+1 << " ";
      os_faces << (*add_f_it).third->get_visu_index()+1 << " ";
      os_faces << " 0 0 0 0" << std::endl; 
    }
  std::cout << "-- medit result dumped." << std::endl;
  delete[] points_tab;
}

//---------------------------------------------------------------------

void
dump_in_file_gv_selected_facets(char* foutput, const Triangulation_3& A)
{ 
  char foutput_tmp[100];
  AF_CGAL_CLIB_STD::strcpy(foutput_tmp, foutput);

  std::strcat(foutput_tmp, ".oogl");
  std::ofstream os(foutput_tmp, std::ios::out);

  if(os.fail())
    std::cerr << "+++unable to open file for output" << std::endl;
  else
    std::cout << ">> file for output : " << foutput_tmp << std::endl;

  os.clear();

  CGAL::set_ascii_mode(os);

  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if ((!v_it->is_exterior()) && (v_it->get_visu_index()==-1))
      {
	v_it->set_visu_index(_vh_number);
	_vh_number++;
      }

  int _vh_size(_vh_number);

  // Header.
  os << "OFF" << std::endl
     << _vh_size << " " << _facet_number << " " << 0 << std::endl;

  typedef const Point*  Const_point_star;
  Const_point_star * points_tab = new  Const_point_star[_vh_size];
  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if (v_it->get_visu_index() > -1)
      points_tab[v_it->get_visu_index()] = &v_it->point();

  for(int vh_i=0; vh_i<_vh_size; vh_i++)
    os << convert()(*points_tab[vh_i]) << std::endl;
     
  for(Finite_facets_iterator f_it = A.finite_facets_begin(); 
      f_it != A.finite_facets_end(); 
      f_it++)
    {
      Cell_handle n, c = (*f_it).first;
      int ni, ci = (*f_it).second;
      n = c->neighbor(ci);
      ni = n->index(c);
      int i1, i2 ,i3;

      if (c->is_selected_facet(ci))
	{
	  i1 = (ci+1) & 3;
	  i2 = (ci+2) & 3;
	  i3 = (ci+3) & 3;
	  os << 3 << " ";
	  os << c->vertex(i1)->get_visu_index() << " ";
	  os << c->vertex(i2)->get_visu_index() << " ";
	  os << c->vertex(i3)->get_visu_index() << " ";
	  os << 0 << std::endl; // without color.
	  // os << 4 << drand48() << drand48() << drand48() << 1.0; // random
	  // color
	}

       if (n->is_selected_facet(ni))
	{
	  i1 = (ni+1) & 3;
	  i2 = (ni+2) & 3;
	  i3 = (ni+3) & 3;
	  os << 3 << " ";
	  os << n->vertex(i1)->get_visu_index() << " ";
	  os << n->vertex(i2)->get_visu_index() << " ";
	  os << n->vertex(i3)->get_visu_index() << " ";
	  os << 0 << std::endl; // without color.
	  // os << 4 << drand48() << drand48() << drand48() << 1.0; // random
	  // color 
	}
    }
  // iterer sur _additional_facets_list pour rajouter les facttes manquantes
  for(Additional_facets_list_iterator add_f_it =
	_additional_facets_list->begin();
      add_f_it != _additional_facets_list->end(); add_f_it++)
    {
      os << 3 << " ";
      os << (*add_f_it).first->get_visu_index() << " ";
      os << (*add_f_it).second->get_visu_index() << " ";
      os << (*add_f_it).third->get_visu_index() << " ";
      os << 0 << std::endl;  // without color.
      // os << 4 << drand48() << drand48() << drand48() << 1.0; // random
      // color 
    }
  std::cout << "-- oogl result dumped." << std::endl;
  delete[] points_tab;
}

//---------------------------------------------------------------------

void
dump_in_file_ply_selected_facets(char* foutput, const Triangulation_3& A)
{ 
  char foutput_tmp[100];
  AF_CGAL_CLIB_STD::strcpy(foutput_tmp, foutput);

  std::strcat(foutput_tmp, ".ply");
  std::ofstream os(foutput_tmp, std::ios::out);

  if(os.fail())
    std::cerr << "+++unable to open file for output" << std::endl;
  else
    std::cout << ">> file for output : " << foutput_tmp << std::endl;

  os.clear();

  CGAL::set_ascii_mode(os);

  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if ((!v_it->is_exterior()) && (v_it->get_visu_index()==-1))
      {
	v_it->set_visu_index(_vh_number);
	_vh_number++;
      }

  int _vh_size(_vh_number);

  // Header.
  os << "ply" << std::endl
     << "format ascii 1.0" << std::endl
     << "comment generated by ply_writer" << std::endl
     << "element vertex " << _vh_size << std::endl
     << "property float x" << std::endl
     << "property float y" << std::endl
     << "property float z" << std::endl
     << "element face " << _facet_number << std::endl
     << "property list uchar int vertex_indices" << std::endl
     << "end_header" << std::endl;

  typedef const Point*  Const_point_star;
  Const_point_star * points_tab = new  Const_point_star[_vh_size];
  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if (v_it->get_visu_index() > -1)
      points_tab[v_it->get_visu_index()] = &v_it->point();

  for(int vh_i=0; vh_i<_vh_size; vh_i++)
    os << convert()(*points_tab[vh_i]) << std::endl;
     
  for(Finite_facets_iterator f_it = A.finite_facets_begin(); 
      f_it != A.finite_facets_end(); 
      f_it++)
    {
      Cell_handle n, c = (*f_it).first;
      int ni, ci = (*f_it).second;
      n = c->neighbor(ci);
      ni = n->index(c);
      int i1, i2 ,i3;

      if (c->is_selected_facet(ci))
	{
	  i1 = (ci+1) & 3;
	  i2 = (ci+2) & 3;
	  i3 = (ci+3) & 3;
	  os << 3 << " ";
	  os << c->vertex(i1)->get_visu_index() << " ";
	  os << c->vertex(i2)->get_visu_index() << " ";
	  os << c->vertex(i3)->get_visu_index();
	  os << std::endl; // without color.
	  // os << 4 << drand48() << drand48() << drand48() << 1.0; // random
	  // color
	}

       if (n->is_selected_facet(ni))
	{
	  i1 = (ni+1) & 3;
	  i2 = (ni+2) & 3;
	  i3 = (ni+3) & 3;
	  os << 3 << " ";
	  os << n->vertex(i1)->get_visu_index() << " ";
	  os << n->vertex(i2)->get_visu_index() << " ";
	  os << n->vertex(i3)->get_visu_index();
	  os << std::endl; // without color.
	  // os << 4 << drand48() << drand48() << drand48() << 1.0; // random
	  // color 
	}
    }
  // iterer sur _additional_facets_list pour rajouter les facttes manquantes
  for(Additional_facets_list_iterator add_f_it =
	_additional_facets_list->begin();
      add_f_it != _additional_facets_list->end(); add_f_it++)
    {
      os << 3 << " ";
      os << (*add_f_it).first->get_visu_index() << " ";
      os << (*add_f_it).second->get_visu_index() << " ";
      os << (*add_f_it).third->get_visu_index();
      os << std::endl;  // without color.
      // os << 4 << drand48() << drand48() << drand48() << 1.0; // random
      // color 
    }
  std::cout << "-- ply result dumped." << std::endl;
  delete[] points_tab;
}

//---------------------------------------------------------------------

void 
dump_in_file_vrml_border_edges(const Triangulation_3& A, std::ofstream& os)
{
  typedef std::pair<Vertex_handle, int>  indiced_vh;
  std::map<Vertex_handle, int> _vh_vect;
  int _vh_bord_count(0);

  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if (v_it->number_of_incident_border() > 0)
      {
	_vh_vect.insert(indiced_vh (v_it->handle(), _vh_bord_count));
	_vh_bord_count++;
      }

  typedef const Point*  Const_point_star;
  Const_point_star * points_tab = new  Const_point_star[_vh_bord_count];
  for (std::map<Vertex_handle, int>::iterator vh_it = _vh_vect.begin();
       vh_it != _vh_vect.end(); vh_it++)
    points_tab[vh_it->second] = &vh_it->first->point();
  
  os << "  Separator {" << std::endl <<
"        Switch {" << std::endl  <<
"          whichChild 0" << std::endl <<
"          Separator {" << std::endl <<
"	     BaseColor {" << std::endl <<
"	                rgb 1 0 0" << std::endl <<
"        	       }" << std::endl << 
"              Coordinate3 {" << std::endl << 
"		 point [ ";
  bool first(true);	
  for(int vh_i=0; vh_i<_vh_bord_count; vh_i++)
    {
      if (!first) os << "," << std::endl <<
"		         ";
      else
	first=false;
      os << convert()(*points_tab[vh_i]);
    }
  os << " ]" << std::endl <<
"		}" << std::endl <<
"		IndexedLineSet {" << std::endl <<
"		  coordIndex [ ";

  first=true;
  for(Finite_edges_iterator e_it=A.finite_edges_begin();
      e_it!=A.finite_edges_end();
      e_it++)
    {
      Cell_handle c = (*e_it).first;
      int i1 = (*e_it).second, i2 = (*e_it).third;
      Edge_like key(c->vertex(i1), c->vertex(i2));
      Border_elt result;

      if (is_border_elt(key, result))
	{
	  if (!first) 
	    os << "," << std::endl << "	 	               ";
	  else
	    first=false;
	  os << _vh_vect.find(c->vertex(i1))->second << ", ";
	  os << _vh_vect.find(c->vertex(i2))->second << ", ";
	  os << -1; 
	}
    }
  os << " ]" << std::endl <<
"      }}" << std::endl <<
"    }}" << std::endl;
  delete[] points_tab;
}

//---------------------------------------------------------------------

void 
dump_in_file_vrml_remaining_points(const Triangulation_3& A, std::ofstream& os)
{
  typedef std::pair<Vertex_handle, int>  indiced_vh;
  std::map<Vertex_handle, int> _vh_vect;
  int _vh_bord_count(0);

  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if (v_it->is_exterior())
      {
	_vh_vect.insert(indiced_vh (v_it->handle(), _vh_bord_count));
	_vh_bord_count++;
      }

  typedef const Point*  Const_point_star;
  Const_point_star * points_tab = new  Const_point_star[_vh_bord_count];
  for (std::map<Vertex_handle, int>::iterator vh_it = _vh_vect.begin();
       vh_it != _vh_vect.end(); vh_it++)
    points_tab[vh_it->second] = &vh_it->first->point();
  
  os << "  Separator {" << std::endl <<
"        Switch {" << std::endl  <<
"          whichChild 0" << std::endl <<
"          Separator {" << std::endl <<
"	     BaseColor {" << std::endl <<
"	                rgb 0 0 1" << std::endl <<
"        	       }" << std::endl << 
"              Coordinate3 {" << std::endl << 
"		 point [ ";
  bool first(true);	
  for(int vh_i=0; vh_i<_vh_bord_count; vh_i++)
    {
      if (!first) os << "," << std::endl <<
"		         ";
      else
	first=false;
      os << convert()(*points_tab[vh_i]);
    }
  os << " ]" << std::endl <<
"		}" << std::endl <<
"		PointSet {" << std::endl <<
"		  startIndex  0" << std::endl <<
"		  numPoints  -1" << std::endl <<
"		}";

  os << " }" << std::endl <<
"      }}" << std::endl;
  delete[] points_tab;
  
}

//---------------------------------------------------------------------
// attention cette procedure produit un fichier tres sale... trop de sommets...

// !!!!  bizarre : ca a l'air de buggue pour hand.xyz (seg fault...)

void 
dump_in_file_vrml_border_facets(const Triangulation_3& A, std::ofstream& os)
{
  typedef std::pair<Vertex_handle, int>  indiced_vh;
  std::map<Vertex_handle, int> _vh_vect;
  int _vh_bord_count(0);

  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
//     if (v_it->number_of_incident_border() > 0)
      {
	_vh_vect.insert(indiced_vh (v_it->handle(), _vh_bord_count));
	_vh_bord_count++;
      }

  typedef const Point*  Const_point_star;
  Const_point_star * points_tab = new  Const_point_star[_vh_bord_count];
  for (std::map<Vertex_handle, int>::iterator vh_it = _vh_vect.begin();
       vh_it != _vh_vect.end(); vh_it++)
    points_tab[vh_it->second] = &vh_it->first->point();
  
  os << "  Separator {" << std::endl <<
"        Switch {" << std::endl  <<
"          whichChild 0" << std::endl <<
"          Separator {" << std::endl <<
"            ShapeHints {" << std::endl <<
"             vertexOrdering  CLOCKWISE" << std::endl <<
"             shapeType       UNKNOWN_SHAPE_TYPE" << std::endl <<
"             faceType        CONVEX" << std::endl <<
"             creaseAngle     1.0" << std::endl <<
"                        }" << std::endl <<
"	     BaseColor {" << std::endl <<
"	                rgb 0 0 1" << std::endl <<
"        	       }" << std::endl << 
"              Coordinate3 {" << std::endl << 
"		 point [ ";
  bool first(true);	
  for(int vh_i=0; vh_i<_vh_bord_count; vh_i++)
    {
      if (!first) os << "," << std::endl <<
"		         ";
      else
	first=false;
      os << convert()(*points_tab[vh_i]);
    }
  os << " ]" << std::endl <<
"		}" << std::endl <<
"		IndexedFaceSet {" << std::endl <<
"		  coordIndex [ ";

  first=true;
  for(Finite_facets_iterator f_it=A.finite_facets_begin();
      f_it!=A.finite_facets_end();
      f_it++)
    {
      Cell_handle c = (*f_it).first;
      int index = (*f_it).second;
      int i1 = (index+1) & 3;
      int i2 = (index+2) & 3;
      int i3 = (index+3) & 3;
      Edge_like key12(c->vertex(i1), c->vertex(i2));
      Edge_like key13(c->vertex(i1), c->vertex(i3));
      Edge_like key32(c->vertex(i3), c->vertex(i2));
      Border_elt result;

      // les trois aretes sur le bord...
//       if (is_border_elt(key12, result)&&
// 	  is_border_elt(key13, result)&&
// 	  is_border_elt(key32, result))

      // au moins 2 aretes sur le bord...
//       if (((is_border_elt(key12, result)&&
// 	    is_border_elt(key13, result)))||
// 	  ((is_border_elt(key32, result)&&
// 	    is_border_elt(key13, result)))||
// 	  ((is_border_elt(key12, result)&&
// 	    is_border_elt(key32, result))))

      // une arete sur le bord...
      if ((is_border_elt(key12, result)||
	   is_border_elt(key13, result)||
	   is_border_elt(key32, result))&&
	  (c->is_selected_facet(index)||
	   c->neighbor(index)->is_selected_facet(c->neighbor(index)->index(c))))

      // au moins 2 aretes sur le bord...
//       if (((is_border_elt(key12, result)&&
// 	    is_border_elt(key13, result)))||
// 	  ((is_border_elt(key32, result)&&
// 	    is_border_elt(key13, result)))||
// 	  ((is_border_elt(key12, result)&&
// 	    is_border_elt(key32, result))))
	{
	  if (!first) 
	    os << "," << std::endl << "	 	               ";
	  else
	    first=false;
	  os << _vh_vect.find(c->vertex(i1))->second << ", ";
	  os << _vh_vect.find(c->vertex(i2))->second << ", ";
	  os << _vh_vect.find(c->vertex(i3))->second << ", ";

	  os << -1; 
	}
    }
  os << " ]" << std::endl <<
"      }}" << std::endl <<
"    }}" << std::endl;	
  delete[] points_tab;
}

//---------------------------------------------------------------------

void
dump_in_file_vrml_selected_facets(char* foutput, const Triangulation_3& A,
				  const bool& contour)
{ 
  char foutput_tmp[100];
  AF_CGAL_CLIB_STD::strcpy(foutput_tmp, foutput);

  std::strcat(foutput_tmp, ".vrml");
  std::ofstream os(foutput_tmp, std::ios::out);

  if(os.fail())
    std::cerr << "+++unable to open file for output" << std::endl;
  else
    std::cout << ">> file for output : " << foutput_tmp << std::endl;

  os.clear();

  CGAL::set_ascii_mode(os);

  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if ((!v_it->is_exterior()) && (v_it->get_visu_index()==-1))
      {
	v_it->set_visu_index(_vh_number);
	_vh_number++;
      }

  int _vh_size(_vh_number);

  // Header.
  os << 
"#Inventor V2.1 ascii" << std::endl <<
"Separator {" << std::endl  <<
"  PerspectiveCamera {" << std::endl <<
"    position 0 0 2.41421" << std::endl <<
"    nearDistance 1.41421" << std::endl <<
"    farDistance 3.41421" << std::endl <<
"    focalDistance 2.41421" << std::endl <<
"  }" << std::endl <<
"  Group {" << std::endl <<
"    Rotation {" << std::endl <<
"             }" << std::endl <<
"    DirectionalLight {" << std::endl <<
"      direction 0.2 -0.2 -0.979796" << std::endl <<
"	 }" << std::endl <<
"    ResetTransform {" << std::endl <<
"	 }  }" << std::endl <<
"  Separator {" << std::endl <<
"    Switch {" << std::endl <<
"      whichChild 0" << std::endl <<
"        Separator {" << std::endl <<
"          ShapeHints {" << std::endl <<
"             vertexOrdering  CLOCKWISE" << std::endl <<
"             shapeType       UNKNOWN_SHAPE_TYPE" << std::endl <<
"             faceType        CONVEX" << std::endl <<
"             creaseAngle     1.0" << std::endl <<
"                      }" << std::endl <<
"	    BaseColor {" << std::endl <<
"	      rgb 0.6 0.6 0.48" << std::endl <<
"        	      }" << std::endl <<
"	    Coordinate3 {" << std::endl <<
"	       point [ ";

  typedef const Point*  Const_point_star;
  Const_point_star * points_tab = new  Const_point_star[_vh_size];
  for (Finite_vertices_iterator v_it = A.finite_vertices_begin();
       v_it != A.finite_vertices_end();
       v_it++)
    if (v_it->get_visu_index() > -1)
      points_tab[v_it->get_visu_index()] = &v_it->point();

  bool first(true);

  for(int vh_i=0; vh_i<_vh_size; vh_i++)
    {
      if (!first) os << "," << std::endl
      << "	               ";
      else
	first=false;
      os << convert()(*points_tab[vh_i]);
    }
  os << " ]" << std::endl <<
"		         }" << std::endl <<
"	     IndexedFaceSet {" << std::endl <<
"	       coordIndex [ ";

  first=true;
  for(Finite_facets_iterator f_it = A.finite_facets_begin(); 
      f_it != A.finite_facets_end(); 
      f_it++)
    {
      Cell_handle n, c = (*f_it).first;
      int ni, ci = (*f_it).second;
      n = c->neighbor(ci);
      ni = n->index(c);
      int i1, i2 ,i3;

      if (c->is_selected_facet(ci))
	{
	  if (!first) os << "," << std::endl <<
"		            ";
	  else
	    first=false;
	  i1 = (ci+1) & 3;
	  i2 = (ci+2) & 3;
	  i3 = (ci+3) & 3;
	  os << c->vertex(i1)->get_visu_index() << ", ";
	  os << c->vertex(i2)->get_visu_index() << ", ";
	  os << c->vertex(i3)->get_visu_index() << ", ";
	  os << -1;
	}

       if (n->is_selected_facet(ni))
	{
	  if (!first) os << "," << std::endl <<
"		            ";
	  else
	    first=false;
	  i1 = (ni+1) & 3;
	  i2 = (ni+2) & 3;
	  i3 = (ni+3) & 3;
	  os << n->vertex(i1)->get_visu_index() << ", ";
	  os << n->vertex(i2)->get_visu_index() << ", ";
	  os << n->vertex(i3)->get_visu_index() << ", ";
	  os << -1; 
	}
    }
  // iterer sur _additional_facets_list pour rajouter les facttes manquantes
  for(Additional_facets_list_iterator add_f_it =
	_additional_facets_list->begin();
      add_f_it != _additional_facets_list->end(); add_f_it++)
    {
      os << "," << std::endl <<
"		            ";
       os << (*add_f_it).first->get_visu_index() << ", ";
       os << (*add_f_it).second->get_visu_index() << ", ";
       os << (*add_f_it).third->get_visu_index() << ", ";
       os << -1;
    }

  os << " ]" << std::endl << 
"      }" << std::endl <<
"    }}" << std::endl << 
"  }" << std::endl;

  if (contour)
    {
      // pour visualiser les contours restant a la fin...
      dump_in_file_vrml_border_edges(A, os);
      
      // pour visualiser les facettes eventuellement candidates...
      //       dump_in_file_vrml_border_facets(A, os);

      // pour afficher les points non selectionnes, ~bruit???
      //      dump_in_file_vrml_remaining_points(A, os);
    }

  os << "}" << std::endl;  

  std::cout << "-- vrml result dumped." << std::endl;
  delete[] points_tab;
}

//---------------------------------------------------------------------

void
dump_in_file_selected_facets(char* foutput, const Triangulation_3& A,
			     const bool& contour, const int& out_format)
{
  for(Finite_vertices_iterator v_it = A.finite_vertices_begin();
      v_it != A.finite_vertices_end(); v_it++)
    v_it->set_visu_index(-1);
  _vh_number = 0;
  switch(out_format)
    {
    case -2:
      // no output file...
      return;
    case -1:      
      dump_in_file_vrml_selected_facets(foutput, A, contour);
      dump_in_file_gv_selected_facets(foutput, A);
      dump_in_file_medit_selected_facets(foutput, A);
      //dump_in_file_ply_selected_facets(foutput, A);
      return;
    case 0:
      return dump_in_file_vrml_selected_facets(foutput, A, contour);
    case 1:
      return dump_in_file_gv_selected_facets(foutput, A);
    case 2:
      return dump_in_file_medit_selected_facets(foutput, A);
    case 3:
      return dump_in_file_ply_selected_facets(foutput, A);
    }
}


//=====================================================================
#endif // IOFILE_MANIPULATOR_H
//=====================================================================
