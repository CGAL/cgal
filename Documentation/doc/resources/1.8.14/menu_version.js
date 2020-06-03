Head:     master Merge branch 'releases/CGAL-5.0-branch'
Merge:    cgal/master Merge branch 'releases/CGAL-5.0-branch'
Push:     cgal/master Merge branch 'releases/CGAL-5.0-branch'
Tag:      master_before_no_tws_nor_tabs (1749)

Untracked files (1)
Documentation/doc/resources/1.8.14/1.8.13

Unstaged changes (1)
modified   Documentation/doc/resources/1.8.13/menu_version.js
@@ -3,10 +3,11 @@

   var url_re =  /(cgal\.geometryfactory\.com\/CGAL\/doc\/|doc\.cgal\.org\/)(master|latest|(\d\.\d+|\d\.\d+\.\d+))\//;
   var url_local =  /.*\/doc_output\//;
-  var current_version_local = '5.1-dev'
+  var current_version_local = '5.1-beta1'
   var all_versions = [
       'master',
       'latest',
+      '5.1-beta1',
       '5.0.2',
       '4.14.3',
       '4.13.2',
@@ -44,7 +45,7 @@
   }

   function patch_url(url, new_version) {
-    if(url.includes("doc.cgal.org")||url.includes("cgal.geometryfactory.com")){
+    if(url.includes("doc.cgal.org")||url.includes("cgal.geometryfactory.com")){
       return url.replace(url_re, 'doc.cgal.org/' + new_version + '/');
     }
     else{
@@ -65,7 +66,7 @@
       var motherNode=$("#back-nav ul")[0];
       var node = document.createElement("LI");
       var spanNode = document.createElement("SPAN");
-      var titleNode =document.createTextNode("CGAL Version: ");
+      var titleNode =document.createTextNode("CGAL Version: ");
       var textNode = document.createTextNode("x.y");
       spanNode.setAttribute("class", "version_menu");
       spanNode.appendChild(textNode);
@@ -90,4 +91,4 @@
         }
      }
   });
-})();
+})();

Stashes (104)
stash@{0} WIP on Mesh_3-fix_optimisers_parallel-jtournois-WIP: 46820eda8bf debugging perturber...
stash@{1} WIP on CGAL-5.0-branch: dc12dea7766 Merge branch 'releases/CGAL-4.14-branch' into releases/CGAL-5.0-branch
stash@{2} On Installation-add_CGALConfigVersion-GF: RELOCATABLE?
stash@{3} On Tetrahedral_remeshing-new-jtournois: DEBUG code
stash@{4} WIP on nurbs-viewer-mesher: 495bc2b6022 Fix the link
stash@{5} On master: point_cloud_to_inr
stash@{6} Patch ExxonMobile
stash@{7} WIP on master: b6d5129364a Merge branch 'releases/CGAL-5.0-branch'
stash@{8} Scene_surface_mesh_item::setAllPatchIds
stash@{9} On master: Replaced by "Update warning macro usages #4474"
stash@{10} On master: Cleanup CDT_plus_2
stash@{11} On master: Fix issue "Warning in Triangulation_2" #4371
stash@{12} On master: f0c82986576 updated crontab (automated commit)
stash@{13} Surface_mesh: use override instead of virtual
stash@{14} On T3_accelerate_insert_in_hole: TDS::reserve
stash@{15} On T3_accelerate_insert_in_hole: Essai d'utilisation des time stamps
stash@{16} On heads/Triangulation_segment_traverser_3-tvanlank__rewrote_history-GF: 7e0f93f4c9e Add a case (assertion)
stash@{17} On master: cleanup Polyline_constraint_hierarchy_2.h
stash@{18} On Sweep_2-bug_57-GF: 02bfdcf828a Better way to measure the recursion depth
stash@{19} On master: 3d1450b71fb Merge branch 'releases/CGAL-4.14-branch'
stash@{20} WIP on master: 3d1450b71fb Merge branch 'releases/CGAL-4.14-branch'
stash@{21} On Polyhedron-demo__add_qtscript_support_to_Mesh_3_plugin-GF: e068363fbfa Add an API to replace QMessageBox static functions
stash@{22} On CGAL-clang_tidy__nullptr_on_Mesh_2-GF: Mesh_3 plugin
stash@{23} On integration: CORE MemoryPool DEBUG
stash@{24} Try to fix check headers warnings about incorrect flags
stash@{25} WIP on releases/CGAL-4.13-branch: a34c09084c6 Merge pull request #3855 from sgiraudot/Intersections_3-Fix_almost_collinear_segments_bug-GF
stash@{26} WIP on releases/CGAL-4.14-branch: 490589b48f6 Merge branch 'releases/CGAL-4.13-branch' into releases/CGAL-4.14-branch
stash@{27} On integration: 32f063f25f0 Merge pull request #3886 from lrineau/CGAL-Adapt_to_Boost_1.70-GF
stash@{28} On integration: GMPXX by default
stash@{29} On integration: Simplify Mesh_3, and MPZF
stash@{30} On Number_types-intervals3-glisse: Include all header twice
stash@{31} On Mesh_3-tricubic-GF: HACKS
stash@{32} On Mesh_3-fix_polyhedral_complex_domain-GF: 9964b18243c Fix Polyhedral_complex_mesh_domain_3 when detect_features() is not called
stash@{33} On Mesh_3-fix_Index-GF: Mesh_3 with union instead of boost::variant
stash@{34} On Mesh_3-fix_Index-GF: Testsuite/test/parse-ctest-dashboard-xml.py with zlib
stash@{35} On master: Tests with -funsafe-math-optimizations
stash@{36} On Triangulation_2-Debug_CDT2-lrineau: Debug CDT_2
stash@{37} On releases/CGAL-4.12-branch: Scene_polygon_soup_item: remove duplicated triangles
stash@{38} On integration: Mesh_3: Ident nested #ifdef
stash@{39} On Travis-Check_including_all_headers-GF: Try to improve "check_headers"
stash@{40} On Polyhedron-demo_offset_plugin_Mesh_3_triangle_soup-GF: attempt of a traversal traits
stash@{41} On Polyhedron-demo_offset_plugin_Mesh_3_triangle_soup-GF: Variable offset distance
stash@{42} On master: Problems with SEP image
stash@{43} WIP on master: 441170768df Prepare CGAL-4.12-beta1
stash@{44} On integration: CSS!!
stash@{45} WIP on Installation-improve_CGAL_DEV_MODE-lrineau: 4811fae85c2 Do not document the CGAL pure header-only mode
stash@{46} Mesh_3: debug trick with CGAL_assertion_msg and a lambda
stash@{47} On integration: CTest
stash@{48} On (no branch): a8704de Merge branch 'releases/CGAL-4.10-branch'
stash@{49} On Mesh_3-fix_NaN-lrineau: DEBUG Sliver_perturber
stash@{50} On Mesh_2-restore_Qt3_demo-lrineau: Doxygen warnings
stash@{51} On integration: Debug Mesh_2
stash@{52} On Mesh_3-API_with_incidences-GF: e17da6445d1 Write the documentation
stash@{53} On Mesh_3-test_polyhedral_complex_with_surface_mesh-GF: TESTS
stash@{54} ?
stash@{55} On Mesh_3-fix_bug_1944-GF: 4325f22 Merge pull request #2145 from gdamiand/patch-1
stash@{56} On integration: a5ea993 Merge remote-tracking branch 'cgal/master' into integration
stash@{57} On Polyhedron_demo-Use_sm_in_Deformation-GF: Selection plugin, pb ODR?
stash@{58} On Polyhedron-clipping_snapping_new_snapping-GF-wip: 23e0762 fixup! Add a safety check
stash@{59} On Polyhedron-clipping_snapping_new_snapping-GF-wip: cbfc043 Revert "wip"
stash@{60} On Polyhedron-clipping_snapping_new_snapping-GF-wip: 720086f Add snap_corners_to_curves_when_possible
stash@{61} On (no branch): Debug Mesh_3
stash@{62} WIP on Polyhedron-clipping_snapping_new_snapping-GF-wip: 50d11a6 Move get_curve_id to Clipping_snapping_tool_details
stash@{63} On Polyhedron-clipping_snapping_new_snapping-GF-wip: WIP on fix corner-to-curve
stash@{64} WIP on Polyhedron-clipping_snapping_new_snapping-GF-wip: 26e2d29 Add missing `#include`
stash@{65} On integration: c540a7e Merge pull request #1808 from MaelRL/Spatial_searching-Fix_fuzzy_query_item_border
stash@{66} WIP on integration: c540a7e Merge pull request #1808 from MaelRL/Spatial_searching-Fix_fuzzy_query_item_border
stash@{67} WIP on integration: c540a7e Merge pull request #1808 from MaelRL/Spatial_searching-Fix_fuzzy_query_item_border
stash@{68} WIP on integration: c540a7e Merge pull request #1808 from MaelRL/Spatial_searching-Fix_fuzzy_query_item_border
stash@{69} WIP on integration: c540a7e Merge pull request #1808 from MaelRL/Spatial_searching-Fix_fuzzy_query_item_border
stash@{70} More assertion in multi-thread Handle.h
stash@{71} CMake: if(NOT CMAKE_CROSSCOMPILING)
stash@{72} Bug Mesh_3 TROU!
stash@{73} WIP on Polyhedron-clipping_snapping_new_snapping-GF-wip: 4f429fa Merge remote-tracking branch 'cgal/master' into Polyhedron-clipping_snapping_new_snapping-GF-wip
stash@{74} On master: convert an image from unsigned short to float
stash@{75} On master: Mesh_3 Robust_weighted_circumcenter_filtered_traits_3 use certainly
stash@{76} On CGAL-remove_support_for_LEDA_5_and_before-GF: 0c69001 Remove all usage of CGAL_LEDA_VERSION
stash@{77} WIP on Mesh_3-new_facet_criterion_with_normals-lrineau: c4b81cf Fix a warning in Sizing_field_with_aabb_tree
stash@{78} On CGAL-license_check-GF: 366976b fix header
stash@{79} Polyhedron demo with keywords for plugins
stash@{80} Mesh_3: save, load a C3t3 and refine it
stash@{81} On master: 2153e65 Fix a typo in doc: remove extra "`"
stash@{82} On Mesh_3-hybrid_mesh_domain-GF: Shifted_sphere_implicit_function
stash@{83} On CGAL_headers_only_step1-gdamiand_cjamin: Fix an error in headers-only
stash@{84} Bug Intel 2017
stash@{85} WIP on Mesh_3-improve_polylines_to_protect-GF: ee0fb3b Fix the header guard macro, and copyright years
stash@{86} On master: Restore CGAL_Qt3
stash@{87} On Mesh_2-fix_issue_781-GF: Mesh_2: Display queues sizes
stash@{88} On Mesh_3-improve_images-GF: "Fix" C3t3_io_plugin to deal with INT_MIN
stash@{89} WIP on Polyhedron-clipping_snapping_new_snapping-GF-wip: 7787959 Re-test `modified_features` after a corner-snapping
stash@{90} On Mesh_3-experimental-GF: Experiments on Mesh-3
stash@{91} On master: Fix orient polygon soup with PMP, temp
stash@{92} On CGAL-Qt5_support-GF: Replace gluErrorString
stash@{93} On master: Mesh_3, for pipeNotWorking
stash@{94} On master: Mesh_3 for Medicim/Nobelbiocare
stash@{95} WIP on Mesh_3-experimental-GF: a1b342e Allow to open a binary CDT_3
stash@{96} WIP on master: 0df4095 Merge pull request #30 from afabri/Documentation-addHome-GF
stash@{97} WIP on CGAL-Qt5_support-GF: 0d2f838 Allow to find QGLViewer-qt5_
stash@{98} On master: LICENSE (Dijsktra), Memory_size.h (near)
stash@{99} On Triangulation_2-Fix_CDT_plus-GF: wip
stash@{100} On master: Pretty printers, et Sliver_perturber.h
stash@{101} On master: bench mesh_3
stash@{102} On Intersection_3-fix-do_intersect_Iso_cuboid_3_Segment_3-lrineau: Fix intersection(Iso_cuboid_3, Segment_3)
stash@{103} On master: Try to fix warnings for CMap

Recent commits
97123fc3c72 cgal/master Merge branch 'releases/CGAL-5.0-branch'
871c97273af Merge pull request #4496 from lrineau/CGAL-move_semantic_for_triangulations-GF
a828cb0d066 Merge pull request #4620 from janetournois/Tetrahedral_remeshing-new-jtournois
8db45039044 Merge pull request #4710 from danston/CGAL-clangmp_bug_fix-danston
c15030bf39b Merge pull request #4740 from afabri/T2-low_dimensional-GF
814689552b1 Merge pull request #4752 from lrineau/CGAL-fix_cpp20-mglisse_GF
863b1decf69 Merge pull request #4754 from maxGimeno/PMP-Fix_parallel_haussdorf_dist-maxGimeno
4354b2c87f0 releases/CGAL-5.0-branch cgal/releases/CGAL-5.0-branch Merge pull request #4710 from danston/CGAL-clangmp_bug_fix-danston
dc2ae1614c2 Merge remote-tracking branch 'cgal/releases/CGAL-4.14-branch' into releases/CGAL-5.0-branch
520fbf7c4b7 refs/pull/4754/head Add missing include
