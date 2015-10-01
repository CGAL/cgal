

static CGAL_THREAD_LOCAL_VARIABLE(int, maxd,30);
static CGAL_THREAD_LOCAL_VARIABLE(std::vector<Face_handle> , f, std::vector<Face_handle>(maxd));
static CGAL_THREAD_LOCAL_VARIABLE(std::vector<int>, i, std::vector<int>(maxd));
static CGAL_THREAD_LOCAL_VARIABLE(std::vector<Vertex_handle>, w, std::vector<Vertex_handle>(maxd));

