#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
using namespace OpenMesh;
using namespace Eigen;

// Returns area of given face on mesh. It's assumed it's a triangle.
double area(Mesh &mesh, FaceHandle fh) {
  //http://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
  Mesh::FaceVertexIter fvit = mesh.fv_iter(fh);
  Vec3f v1 = mesh.point(fvit.handle());
  Vec3f v2 = mesh.point((++fvit).handle());
  Vec3f v3 = mesh.point((++fvit).handle());
  assert (!(++fvit));
  Vec3f ab = v2 - v1;
  Vec3f ac = v3 - v1;
  return (ab % ac).length() / 2.0; //% = cross product
}

void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature) {

  for (Mesh::VertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
    // WRITE CODE HERE TO COMPUTE THE CURVATURE AT THE CURRENT VERTEX ----------------------------------------------
    Vec3f normal = mesh.normal(it.handle());
    Vector3d N(normal[0],normal[1],normal[2]); // example of converting to Eigen's vector class for easier math

    // In the end you need to fill in this struct

    // Find Nvi: normal vector at this vertex (vi). First need to find total area of adjacent faces.
    Vec3f Nvi = mesh.normal(it.handle());

    Vec3f vi = mesh.point(it.handle());
    for(Mesh::VertexVertexIter vvit = mesh.vv_iter(it.handle()); vvit; ++vvit){
      Vec3f vj = mesh.point(vvit.handle());
      Nvi * Nvi.transpose();
      Nvi.identity();
    }
          
    CurvatureInfo info;
    info.curvatures[0] = 0;
    info.curvatures[1] = 0;
    info.directions[0] = Vec3f();
    info.directions[1] = Vec3f();

    mesh.property(curvature,it) = info;
    // -------------------------------------------------------------------------------------------------------------
  }
}
