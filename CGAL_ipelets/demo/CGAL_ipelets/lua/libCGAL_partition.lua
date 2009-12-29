----------------------------------------------------------------------
-- CGAL Polygon Partition ipelet description
----------------------------------------------------------------------

label = "Polygon Partition"

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
  { label="Y monotone partition" },
  { label="Greene's approx Convex Partition"},
  { label="Approx Convex Partition"},
  { label="Optimal Convex Partition"},
  { label="Help" },
}

----------------------------------------------------------------------
