#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
using namespace OpenMesh;
using namespace Eigen;

void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature) {

const float threshold_parallel = 0.7;

for (Mesh::VertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
// WRITE CODE HERE TO COMPUTE THE CURVATURE AT THE CURRENT VERTEX ----------------------------------------------
Vec3f normal = mesh.normal(it.handle());
Vector3d N(normal[0],normal[1],normal[2]); // example of converting to Eigen's vector class for easier math

// Modified by forestj@stanford

// Sum of the area of surrounding triangles, used in weight denominator
float totalArea = 0.0;
Mesh::VertexOHalfedgeIter voh_iter;
for (voh_iter = mesh.voh_iter(it.handle()); voh_iter; ++voh_iter)
totalArea += mesh.calc_sector_area(voh_iter.handle());

// Prepare the matrix
Matrix3d M = Matrix3d::Zero();

// Neighbor access
for (voh_iter = mesh.voh_iter(it.handle()); voh_iter; ++voh_iter)
{
// Opposite HE
HalfedgeHandle oppHEH = mesh.opposite_halfedge_handle(voh_iter.handle());
// assume the mesh is without boundary
assert(!is_boundary(voh_iter.handle()));
assert(!is_boundary(oppHEH));

// Weight is the sum of area of neighboring triangles
float weight = 0.0;
weight += mesh.calc_sector_area(voh_iter.handle());
weight += mesh.calc_sector_area(oppHEH);
weight = weight / (2 * totalArea);

// Edge vector
Vec3f edge;
mesh.calc_edge_vector(voh_iter.handle(),edge);
Vector3d vecE(edge[0],edge[1],edge[2]);
// kappa
float kappa = (2/(vecE.norm()*vecE.norm())) * N.transpose() * vecE ;
// T vector
Vector3d vecT;
vecT = vecE - N * N.transpose() * vecE;
vecT.normalize();
// sum into matrix M
M += weight * kappa * vecT * vecT.transpose();
}


// ------- Eigen decomposition for the matrix M
EigenSolver<Matrix3d> solver(M);
Vector3d v1 = solver.pseudoEigenvectors().block(0,0,3,1);
Vector3d v2 = solver.pseudoEigenvectors().block(0,1,3,1);
Vector3d v3 = solver.pseudoEigenvectors().block(0,2,3,1);

double eig1 = real(solver.eigenvalues()(0));
double eig2 = real(solver.eigenvalues()(1));
double eig3 = real(solver.eigenvalues()(2));


// Principal curvatures and directions. Consider different possible situations.
double curv1,curv2; // curv1 >= curv2
Vector3d pDir1,pDir2; // principal directions
if (abs(N.dot(v1.normalized())) > threshold_parallel)
{
if (eig2 >= eig3){
curv1 = 3*eig2 - eig3;curv2 = 3*eig3 - eig2;
pDir1 = v2;pDir2 = v3;
}else{
curv1 = 3*eig3 - eig2;curv2 = 3*eig2 - eig3;
pDir1 = v3;pDir2 = v2;
}
}else if (abs(N.dot(v2.normalized())) > threshold_parallel)
{
if (eig1 >= eig3){
curv1 = 3*eig1 - eig3;curv2 = 3*eig3 - eig1;
pDir1 = v1;pDir2 = v3;
}else{
curv1 = 3*eig3 - eig1;curv2 = 3*eig1 - eig3;
pDir1 = v3;pDir2 = v1;
}
}else if (abs(N.dot(v3.normalized())) > threshold_parallel){
if (eig1 >= eig2){
curv1 = 3*eig1 - eig2;curv2 = 3*eig2 - eig1;
pDir1 = v1;pDir2 = v2;
}else{
curv1 = 3*eig2 - eig1;curv2 = 3*eig1 - eig2;
pDir1 = v2;pDir2 = v1;
}
}else {assert(false);}// shouldn't happen

// In the end you need to fill in this struct

CurvatureInfo info;
info.curvatures[0] = curv1;
info.curvatures[1] = curv2;
info.directions[0] = Vec3f(pDir1[0],pDir1[1],pDir1[2]);
info.directions[1] = Vec3f(pDir2[0],pDir2[1],pDir2[2]);

mesh.property(curvature,it) = info;
// -------------------------------------------------------------------------------------------------------------
}
}
