#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
using namespace OpenMesh;
using namespace Eigen;

using namespace std;

const double EPSILON = 0.00001;

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
    Mesh::VertexHandle vi_handle = it.handle();

    // In the end you need to fill in this struct

    // Find Nvi: normal vector at this vertex (vi). First need to find total area of adjacent faces.
    Vec3f mesh_Nvi = mesh.normal(vi_handle);
    Vec3f mesh_vi = mesh.point(vi_handle);
    Vector3d Nvi(mesh_Nvi[0], mesh_Nvi[1], mesh_Nvi[2]);
    Vector3d vi(mesh_vi[0], mesh_vi[1], mesh_vi[2]);

    //compute total area of triangles around vi
    double areaSum = 0;
    for(Mesh::VertexFaceIter vfit = mesh.vf_iter(it.handle()); vfit; ++vfit){
      areaSum += area(mesh, vfit.handle());
    }

    //compute the matrix Mvi
    // Get the vertex-outgoing halfedges circulator of vertex _vh
    double wijSum = 0;
    Matrix3d Mvi = Matrix3d::Zero();
    for(Mesh::VertexOHalfedgeIter vohit = mesh.voh_iter(it.handle()); vohit; ++vohit){
      Mesh::VertexHandle vj_handle = mesh.to_vertex_handle(vohit.handle());
      assert(vj_handle != vi_handle);
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
      assert (fh1 != fh2);
      double wij = (area(mesh, fh1) + area(mesh, fh2)) / (2 * areaSum);
      wijSum += wij;

      Mvi += Tij * Tij.transpose() * wij * Kij;
    }
    assert(wijSum - 1.0 < EPSILON);

    //get eigenstuff
    EigenSolver<Matrix3d> es(Mvi, true);

    //TODO: check that there are no complex parts
    assert (abs(imag(es.eigenvalues()[0]) < EPSILON));

    //Find the bigger, smaller, and 0 eigenvalue, and make those T1, T2, N
    int i_t1 = 0, i_t2 = 1, i_n = 2;
    double eVals[3] = {abs(real(es.eigenvalues()[0])),
                       abs(real(es.eigenvalues()[1])),
                       abs(real(es.eigenvalues()[2]))};;
    if (eVals[0] < eVals[1] && eVals[0] < eVals[2]){
      assert (eVals[0] < EPSILON); //It's the zero
      i_n = 0;

      if (eVals[1] > eVals[2]) {
        i_t1 = 1;
        i_t2 = 2;
      }
      else {
        i_t1 = 2;
        i_t2 = 1;
      }
    }
    else if(eVals[1] < eVals[0] && eVals[1] < eVals[2]) {
      assert (eVals[1] < EPSILON); //It's the zero
      i_n = 1;

      if (eVals[0] > eVals[2]) {
        i_t1 = 0;
        i_t2 = 2;
      }
      else {
        i_t1 = 2;
        i_t2 = 0;
      }
    }
    else {
      assert (eVals[2] < EPSILON); //It's the zero
      i_n = 2;

      if (eVals[1] > eVals[0]) {
        i_t1 = 1;
        i_t2 = 0;
      }
      else {
        i_t1 = 0;
        i_t2 = 1;
      }
    }

    //Assign curvatures and directions
    CurvatureInfo info;
    info.curvatures[0] = eVals[i_t1];
    info.curvatures[1] = eVals[i_t2];
    info.directions[0] = Vec3f(real((es.eigenvectors().col(i_t1))[0]),
                               real((es.eigenvectors().col(i_t1))[1]),
                               real((es.eigenvectors().col(i_t1))[2]));
    info.directions[1] = Vec3f(real((es.eigenvectors().col(i_t2))[0]),
                               real((es.eigenvectors().col(i_t2))[1]),
                               real((es.eigenvectors().col(i_t2))[2]));

    mesh.property(curvature,it) = info;
    // -------------------------------------------------------------------------------------------------------------
  }
}
