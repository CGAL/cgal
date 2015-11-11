

static CGAL_THREAD_LOCAL_VARIABLE(int, maxd,30);
static CGAL_THREAD_LOCAL_VARIABLE(std::vector<Face_handle> , f, maxd);
static CGAL_THREAD_LOCAL_VARIABLE(std::vector<int>, i, maxd);
static CGAL_THREAD_LOCAL_VARIABLE(std::vector<Vertex_handle>, w, maxd);

