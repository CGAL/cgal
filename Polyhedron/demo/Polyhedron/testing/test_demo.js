var testItem = "./testing/data/sphere.off";
main_window.open(testItem, 'off_plugin');
main_window.test_all_actions();
testItem = "./testing/data/oni.pwn";
main_window.open(testItem, 'xyz_plugin');
main_window.test_all_actions();

var itemA = "./testing/data/itema.off";
var itemB = "./testing/data/itemb.off";
main_window.open(itemA, 'off_to_nef_plugin');
main_window.open(itemB, 'off_to_nef_plugin');
scene.setItemA(0);
scene.setItemB(1);
main_window.test_all_actions();

testItem = "./testing/data/poly.polylines.txt";
main_window.open(testItem, 'polylines_io_plugin');
main_window.test_all_actions();

testItem = "./testing/data/sphere.mesh";
main_window.open(testItem, 'C3t3_io_plugin');
main_window.test_all_actions();

testItem = "./testing/data/camera.camera.txt";
main_window.open(testItem, 'camera_positions_plugin');

testItem = "./testing/data/sphere.ts";
main_window.open(testItem, 'gocad_plugin');
scene.erase(0);

testItem = "./testing/data/sphere.nef3";
main_window.open(testItem, 'io_nef_plugin');
scene.erase(0);

testItem = "./testing/data/sphere.ply";
main_window.open(testItem, 'ply_plugin');
scene.erase(0);

testItem = "./testing/data/oni.ply";
main_window.open(testItem, 'ply_plugin');
scene.erase(0);

testItem = "./testing/data/sphere.stl";
main_window.open(testItem, 'stl_plugin');
scene.erase(0);

if (typeof vtk_plugin !== "undefined") {
  testItem = "./testing/data/sphere.vtk";
  main_window.open(testItem, 'vtk_plugin');
  scene.erase(0);

  testItem = "./testing/data/sphere.vtp";
  main_window.open(testItem, 'vtk_plugin');
  scene.erase(0);

  testItem = "./testing/data/sphere.vtu";
  main_window.open(testItem, 'vtk_plugin');
  scene.erase(0);
}

testItem = "./testing/data/sphere.off";
main_window.open(testItem, 'off_plugin');
testItem = "./testing/data/sphere.selection.txt";
main_window.open(testItem, 'selection_sm_plugin');
scene.erase(1);
scene.erase(0);

testItem = "./testing/data/mini.surf";
main_window.open(testItem, 'surf_io_plugin');
scene.erase(3); //id of the group contaning the items.

testItem = "./testing/data/sphere.inr";
main_window.open(testItem, 'segmented images');
scene.erase(0);

testItem = "./testing/data/sphere.inr.gz";
main_window.open(testItem, 'segmented images');
scene.erase(0);

testItem = "./testing/data/sphere.inr.gz";
main_window.open(testItem, 'segmented images');
scene.erase(0);

if (typeof las_plugin !== "undefined") {
    testItem = "./testing/data/oni.las";
    main_window.open(testItem, 'las_plugin');
    scene.erase(0);
}
