----------------------------------------------------------------------
-- CGAL multi regular ipelet description
----------------------------------------------------------------------

label = "k-order Regular"

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
  { label="Regular" },
  { label="Regular 2" },
  { label="Regular 3" },
  { label="Regular n-1" },
  { label="Regular k" },
  { label="Power Diagram" },
  { label="Power Diagram 2" },
  { label="Power Diagram 3" },
  { label="Power Diagram n-1" },
  { label="Power Diagram k" },
  { label="Help" },
}

----------------------------------------------------------------------
