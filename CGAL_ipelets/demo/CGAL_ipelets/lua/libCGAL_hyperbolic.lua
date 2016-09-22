----------------------------------------------------------------------
-- CGAL hyperbolic geometry ipelet description
----------------------------------------------------------------------

label = "Hyperbolic"

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
  { label= "Line through two points" },
  { label= "Segment through two points" },
  { label= "Bisector of two points" },
  { label= "Circle by center(prim. sel.) and point" },
  { label= "Circle center" },
  { label= "Octagon group translations" },
  { label= "Translation a" },
  { label= "Translation inv-b" },
  { label= "Translation c" },
  { label= "Translation inv-d" },
  { label= "Translation inv-a" },
  { label= "Translation b" },
  { label= "Translation inv-c" },
  { label= "Translation d" },
  { label= "Test" },
  { label= "Help" },
}

----------------------------------------------------------------------
