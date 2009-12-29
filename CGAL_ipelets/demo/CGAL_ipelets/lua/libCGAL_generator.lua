----------------------------------------------------------------------
-- CGAL generators ipelet description
----------------------------------------------------------------------

label = "Generators"

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
  { label="Points in a disk" },
  { label="Points on a grid" },
  { label="Points in a square" },
  { label="Points on a convex hull" },
  { label="Polygon" },
  { label="Segments in a square" },
  { label="Circles (center in a square)" },
  { label="Help" },
}

----------------------------------------------------------------------
