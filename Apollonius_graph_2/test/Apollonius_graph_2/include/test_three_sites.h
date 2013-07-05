#ifndef TEST_AG_THREE_SITES_H
#define TEST_AG_THREE_SITES_H 1

typedef enum { MIDDLE_ON_CONVEX_HULL, MIDDLE_NOT_ON_CONVEX_HULL } Middle_t;

template<class AG>
int number_of_infinite_edges(const AG& ag, typename AG::Vertex_handle v)
{
  int count(0);

  typename AG::Edge_circulator ec_start = ag.incident_edges(v);
  typename AG::Edge_circulator ec = ec_start;

  do {
    if ( ag.is_infinite(*ec) ) { ++count; }
    ++ec;
  } while ( ec != ec_start );

  return count;
}

template<class AG>
int number_of_infinite_faces(const AG& ag, typename AG::Vertex_handle v)
{
  int count(0);

  typename AG::Face_circulator fc_start = ag.incident_faces(v);
  typename AG::Face_circulator fc = fc_start;

  do {
    if ( ag.is_infinite(fc) ) { ++count; }
    ++fc;
  } while ( fc != fc_start );

  return count;
}

template<class AG>
int number_of_infinite_vertices(const AG& ag, typename AG::Vertex_handle v)
{
  int count(0);

  typename AG::Vertex_circulator vc_start = ag.incident_vertices(v);
  typename AG::Vertex_circulator vc = vc_start;

  do {
    if ( ag.is_infinite(vc) ) { ++count; }
    ++vc;
  } while ( vc != vc_start );

  return count;
}


template<class AG>
bool test_three_sites(const AG& ag,
		      typename AG::Vertex_handle v1,
		      typename AG::Vertex_handle v2,
		      typename AG::Vertex_handle v3,
		      Middle_t mid_type)
{
  std::cout << "Site 1: " << v1->site() << std::endl;
  std::cout << "Site 2: " << v2->site() << std::endl;
  std::cout << "Site 3: " << v3->site() << std::endl;
  std::cout << std::endl;

  if ( !ag.is_valid() ) { return false; }

  if ( number_of_infinite_vertices(ag, v1) != 1 ) { return false; }
  if ( number_of_infinite_vertices(ag, v3) != 1 ) { return false; }

  if ( number_of_infinite_edges(ag, v1) != 1 ) { return false; }
  if ( number_of_infinite_edges(ag, v3) != 1 ) { return false; }

  if ( number_of_infinite_faces(ag, v1) != 2 ) { return false; }
  if ( number_of_infinite_faces(ag, v3) != 2 ) { return false; }

  if ( mid_type == MIDDLE_ON_CONVEX_HULL ) {
    if ( number_of_infinite_vertices(ag, v2) != 2 ) { return false; }
    if ( number_of_infinite_edges(ag, v2) != 2 ) { return false; }
    if ( number_of_infinite_faces(ag, v2) != 4 ) { return false; }

    if ( ag.tds().degree(v1) != 2 ) { return false; }
    if ( ag.tds().degree(v3) != 2 ) { return false; }
    if ( ag.tds().degree(v2) != 4 ) { return false; }

    if ( ag.tds().degree(ag.infinite_vertex()) != 4 ) { return false; }
  } else {
    if ( number_of_infinite_vertices(ag, v2) != 0 ) { return false; }
    if ( number_of_infinite_edges(ag, v2) != 0 ) { return false; }
    if ( number_of_infinite_faces(ag, v2) != 0 ) { return false; }

    if ( ag.tds().degree(v1) != 4 ) { return false; }
    if ( ag.tds().degree(v3) != 4 ) { return false; }
    if ( ag.tds().degree(v2) != 2 ) { return false; }

    if ( ag.tds().degree(ag.infinite_vertex()) != 2 ) { return false; }
  }

  return true;
}

#endif // TEST_AG_THREE_SITES_H
