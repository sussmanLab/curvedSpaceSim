#include "submesher.h"
#include <set>
#include <stack>

triangleMesh submesher::constructSubmeshFromSourceAndTargets(triangleMesh &mesh, faceLocation &source, std::vector<faceLocation> &targets, double &maximumDistanceFromSource)
    {
    triangleMesh submesh;
	/*
  Triangle_mesh submesh; 
  Point_3 source_r3 = PMP::construct_point(source,mesh); 

  Face_index source_face = source.first;
   
  std::set<Face_index> goal_faces; 
  
  std::set<Face_index> visited; 
  visited.insert(source_face);
  
  // first: add to the list of goal faces all the faces which contain 
  // a target but are not the source face
  for (Face_location target: targets) {
    //std::cout << PMP::construct_point(target, mesh) << " in face " << target.first << std::endl;  
    if (target.first != source_face) goal_faces.insert(target.first); 
  }
  if (goal_faces.empty()) return create_submesh_from_visited(visited,mesh); // early return if all our targets are in one face

  std::vector<Face_index> exploration_stack; 
  std::vector<Point_3> face_verts; 
  face_verts.reserve(3);
  Face_index neighboring_face;
  
  Halfedge_index hf = mesh.halfedge(source_face);

  for(Halfedge_index hi : halfedges_around_face(hf, mesh)){
    // Access the neighboring face through the opposite halfedge of the current halfedge
    neighboring_face = mesh.face(mesh.opposite(hi));
    face_verts = getVertexPositions(mesh, neighboring_face);
    visited.insert(neighboring_face); 
    exploration_stack.push_back(neighboring_face);
    
    if(goal_faces.count(neighboring_face)) {
      goal_faces.erase(neighboring_face);
      // since we know there can't be holes in a mesh with only faces that are directly connected 
      // to the source face: 
      if (goal_faces.empty()) {
        return create_submesh_from_visited(visited, mesh); 
      } 
    } 
    // push to queue to investigate neighbors so long as we still have targets to find
  }

  Face_index current_face;
  Triangle_3 comp_triangle;

  // now: breadth-first search, checking every connected face and trying to make a complete patch
  // while (!exploration_stack.empty() && !goal_faces.empty()) { // early stop condition -- faster for small # of targets but possible to get holes 
  while (!exploration_stack.empty()) {   
    // get the current face we're "standing" in    
    current_face = exploration_stack.back(); 
    exploration_stack.pop_back(); 
    // redefining hf since it will be serving the same purpose
    hf = mesh.halfedge(current_face); 

    for (Halfedge_index hi : halfedges_around_face(hf, mesh)) {
      neighboring_face = mesh.face(mesh.opposite(hi));
      
      // if we've seen this face before, continue so we don't create a closed loop 
      if(visited.count(neighboring_face)) { 
        continue;
      }
      
      face_verts = getVertexPositions(mesh,neighboring_face);
      // if all our vertices are outside the cutoff radius, we're out of bounds
      // and shouldn't add this face to the exploration stack
      bool out = false;
      comp_triangle = Triangle_3(face_verts[0],face_verts[1],face_verts[2]);

      if (sqrt(CGAL::squared_distance(source_r3, comp_triangle)) > cutoff_dist) {
       	continue;
      }
       
      // the two statements above guarantee closure, although closed surfaces 
      // with small internal "radii"/large cutoff radii can lead to delayed closure/
      // submeshes that are too large becuase they will require the first condition

      visited.insert(neighboring_face); 
      
      // unit test -- in the future, would be cool to make this only happen in debug mode 
      if(goal_faces.count(neighboring_face)) {
        goal_faces.erase(neighboring_face); 
      }
      
      //finally, if we absolutely need to keep looking down this path (
      //it is not past the neighbor range, we haven't seen it before), we do. 
      exploration_stack.push_back(neighboring_face);
    }    

  }
  if (!goal_faces.empty()) {
    std::cout << "All goal faces not found, debug!" << std::endl;
    printf("remaining goal faces: \n");
    for (Face_index face: goal_faces) std::cout << face << " ";
  }
  printf("\n");


  submesh = create_submesh_from_visited(visited,mesh); 
  */
    return submesh;
    };
