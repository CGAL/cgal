-- ipe2cin.lua
-- convert points and segments in ipe file to cin format,
-- only for first page and its first view of ipe file
-- and only for visible layers in the first view

if #argv ~= 1 then
  io.stderr:write("Usage: ipescript ipe2cin <file>\n")
  return
end
fname = argv[1]

doc = assert(ipe.Document(fname))

-- io.write("# number of pages = ", #doc, "\n")
-- print("# of pages =", #doc)

p = assert(doc[1])

-- print("# of views=", p:countViews())

-- print("# of objects in p1 =", #p)

function process_object(obj)
  -- print(obj:matrix())
  m = obj:matrix()
  is_id = m:isIdentity()
  if obj:type() == "reference" then
    if (obj:get("markshape")):find("mark") then
      if is_id then
        vec = obj:position()
      else
        vec = m * obj:position()
      end
      print("p", vec.x, vec.y)
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
            if is_id then
              endp1 = seg[1]
              endp2 = seg[2]
            else
              endp1 = m * seg[1]
              endp2 = m * seg[2]
            end
            print("s", endp1.x, endp1.y, endp2.x, endp2.y)
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
          if is_id then
            endp1 = segl[2]
            endp2 = segf[1]
          else
            endp1 = m * segl[2]
            endp2 = m * segf[1]
          end
          print("s", endp1.x, endp1.y, endp2.x, endp2.y)
        end
      end
    end
    --print(obj:xml())
  end
end

for i, obj, sel, layer in p:objects() do
  --print(i, obj, sel, layer)
  if p:visible(1, layer) then
    process_object(obj)
  end
end
