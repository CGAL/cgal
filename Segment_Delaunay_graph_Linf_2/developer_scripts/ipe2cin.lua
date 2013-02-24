-- ipe2cin.lua
-- convert points and segments in ipe file to cin format
-- (only for first page of ipe file)

if #argv ~= 1 then
  io.stderr:write("Usage: ipescript ipe2cin <file>\n")
  return
end
fname = argv[1]

doc = assert(ipe.Document(fname))

-- io.write("# number of pages = ", #doc, "\n")
-- print ("# of pages =", #doc)

p = assert(doc[1])

-- print ("# of objects in p1 =", #p)

-- iterate over objects of page:
for i, obj, sel, layer in p:objects() do
  --print(i, obj, sel, layer)
  if obj:type() == "reference" then
    if ((obj:get("markshape")):find("mark")) then
      print ("p", obj:position().x, obj:position().y)
    end
  end
  if obj:type() == "path" then
    allsegments = true
    shape = obj:shape()
    for i, subpath in pairs(shape) do
      if subpath["type"] == "curve" then
        manycomponents = false
        for i, seg in ipairs(subpath) do
          if i == 1 then
            segf = seg
          end
          if seg["type"] == "segment" then
            print("s", seg[1].x, seg[1].y, seg[2].x, seg[2].y)
          else
            allsegments = false
          end
          if i > 1 then
            manycomponents = true
            segl = seg
          end
        end
        if subpath["closed"] and allsegments and manycomponents then
          -- connect last point of last segment to
          -- first point of first segment
         print("s", segl[2].x, segl[2].y, segf[1].x, segf[1].y)
        end
      end
    end
    --print(obj:xml())
  end
end

