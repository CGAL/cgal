----------------------------------------------------------------------
-- WSPD ipelet description
----------------------------------------------------------------------

label = "Cone Spanners"

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
  { label="Theta-k-graph" },
  { label="Yao-k-graph" },
  { label="Half-theta-k-graph with even cones" },
  { label="Half-Yao-k-graph with even cones" },
  { label="Half-theta-k-graph with odd cones" },
  { label="Half-Yao-k-grap with odd cones" },
  { label="k cones" },
  { label="Help" },
}

----------------------------------------------------------------------
