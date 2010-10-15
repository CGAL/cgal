----------------------------------------------------------------------
-- CGAL Triangulations ipelet description
----------------------------------------------------------------------

label = "Triangulations"

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
  { label="Delaunay" },
  { label="Constrained Delaunay" },
  { label="Conforming Delaunay" },
  { label="Conforming Gabriel" },
  { label="Regular" },
  { label="Help" },
}

----------------------------------------------------------------------
