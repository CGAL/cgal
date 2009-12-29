----------------------------------------------------------------------
-- CGAL multi delaunay ipelet description
----------------------------------------------------------------------

label = "k order Delaunay"

about = [[
This ipelet is part of the CGAL_ipelet package. See www.cgal.org.
]]

-- this variable will store the C++ ipelet when it has been loaded
ipelet = false

function run(model, num)
  if not ipelet then ipelet = assert(ipe.Ipelet(dllname)) end
  model:runIpelet(methods[num].label, ipelet, num)
end

methods = {
  { label= "Delaunay" },
  { label= "Delaunay 2" },
  { label= "Delaunay 3" },
  { label= "Delaunay n-1" },
  { label= "Delaunay k" },
  { label= "Voronoi" },
  { label= "Voronoi 2" },
  { label= "Voronoi 3" },
  { label= "Voronoi n-1" },
  { label= "Voronoi k" },
  { label="Help" },
}

----------------------------------------------------------------------
