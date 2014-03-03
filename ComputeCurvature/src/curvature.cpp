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

    // In the end you need to fill in this struct

    // Find Nvi: normal vector at this vertex (vi). First need to find total area of adjacent faces.
    Vec3f mesh_Nvi = mesh.normal(it.handle());
    Vec3f mesh_vi = mesh.point(it.handle());
    Vector3d Nvi(mesh_Nvi[0], mesh_Nvi[1], mesh_Nvi[2]);
    Vector3d vi(mesh_vi[0], mesh_vi[1], mesh_vi[2]);

    //compute total area of triangles around vi
    double areaSum = 0;
    for(Mesh::VertexFaceIter vfit = mesh.vf_iter(it.handle()); vfit; ++vfit){
      areaSum += area(mesh, vfit.handle());
    }

    //compute the matrix Mvi
    // Get the vertex-outgoing halfedges circulator of vertex _vh
    Matrix3d Mvi = Matrix3d::Zero();
    for(Mesh::VertexOHalfedgeIter vohit = mesh.voh_iter(it.handle()); vohit; ++vohit){
      Mesh::VertexHandle vj_handle = mesh.to_vertex_handle(vohit.handle());
      Vec3f mesh_vj = mesh.point(vj_handle);
      Vector3d vj(mesh_vj[0], mesh_vj[1], mesh_vj[2]);

      // Calculate Kij
      Vector3d edge = vj - vi;
      double Kij = 2.0 * Nvi.transpose() * edge;
      Kij /= edge.squaredNorm();

      // Calculate Tij
      Matrix3d I = Matrix3d::Identity();
      Vector3d Tij = (I - Nvi * Nvi.transpose()) * (vi - vj);
      Tij.normalize();

      // Calculate wij
      //faces on both sides of halfedge
      Mesh::FaceHandle fh1 = mesh.face_handle(vohit.handle());
      Mesh::FaceHandle fh2 = mesh.opposite_face_handle(vohit.handle());
      double wij = (area(mesh, fh1) + area(mesh, fh2)) / (2 * areaSum);

      Mvi += wij * Kij * Tij * Tij.transpose();
    }

    //TODO: get eigenstuff and put it in the right place
    EigenSolver<Matrix3d> es(Mvi, true);

    Vector3d eVec1 = es.pseudoEigenvectors().block(0,0,3,1),
      eVec2 = es.pseudoEigenvectors().block(0,1,3,1),
      eVec3 = es.pseudoEigenvectors().block(0,2,3,1);

    double eVal1 = real(es.eigenvalues()(0)),
      eVal2 = real(es.eigenvalues()(1)),
      eVal3 = real(es.eigenvalues()(2));

    // Find which eigenvector is Nvi
    float thresh = 0.0001;
    if (abs(Nvi.dot(eVec1.normalized())) > thresh) {
      // eVec1 == Nvi
      std::cout << "eVec1 = Nvi" << std::endl;
    } else if (abs(Nvi.dot(eVec2.normalized())) > thresh) {
      // eVec2 == Nvi
      std::cout << "eVec2 = Nvi" << std::endl;
    } else if (abs(Nvi.dot(eVec3.normalized())) > thresh) {
      // eVec3 == Nvi
      std::cout << "eVec3 = Nvi" << std::endl;
    } else {
      // Nvi not an eigenvector => too high threshold
      std::cout << "Nvi not found in eigenvectors" << std::endl;
    }

    //TODO: Assign curvatures and directions
    CurvatureInfo info;
    info.curvatures[0] = 1;
    info.curvatures[1] = 1;
    info.directions[0] = Vec3f(1,0,0);
    info.directions[1] = Vec3f(0,0,1);

    mesh.property(curvature,it) = info;
    // -------------------------------------------------------------------------------------------------------------
  }
}
