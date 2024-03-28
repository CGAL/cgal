TODO:

- Try to add a custom exact queue that supports EPECK, because now we convert to doubles there.

- What about accumulating errors in EPICK? Maybe we could use hybrid mode doing only important computations exactly like intersections e.g. Another idea is to compute intersections with respect to the original input data instead of the data obtained from the previous iteration.

- Fix case with touching along an edge polygons at the initialization step, they should propagate while now they do not. See e.g. the edge-case-test/test-same-time.off data set.

- Can we avoid any inexact computations such as sqrt completely?

- Fix stress-test-6 cases - they are random, their results may depend on an initial unlucky configuration so we probably need to implement random perturbation before running them as in the original paper.

- Should we also add random perturbation that can be on or off depending on the input.

- EPICK fails: 40 polygons k = 1 and coplanarity = 0.1 or 40 polygons k = 6 and coplanarity = 0.5; 40 polygons - the one from real-data-test; coplanarity is from Support_plane.h operator==().

- EPECK stress-test-5/test-2-rnd-polygons-20-4 runs extremely slow, does it actually finish? Btw real data with the same number of polygons work much faster! Maybe we better test it from now on only on real data because this is what is going to be most-likely used in real life. Can this speed issue be related to conversion from exact to inexact when creating a queue? It seems to hang at the end of each event that is when we check and continue execution with the next event. Maybe the event queue access is so slow? Do we utilize the technique from the paper, storing and accessing only the last three events? The values maybe grow so much that each intersection takes very long time if we use exact type everywhere.

- Try to avoid errors by using a smaller step.

- Make intersections a part of kinetic traits that is exact.

- Try to find out where we lose precision the most.

- In naive hybrid mode, we manage to succeed for more events, which means that making more exact computations improves the algorithm, that is a good sign,
however doing them all exactly is very slow, so the trade-off must be found.