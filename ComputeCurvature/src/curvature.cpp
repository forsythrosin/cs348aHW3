#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include "curvature.h"
using namespace OpenMesh;
using namespace Eigen;

void computeCurvature(Mesh &mesh, OpenMesh::VPropHandleT<CurvatureInfo> &curvature) {

	for (Mesh::VertexIter it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
		// WRITE CODE HERE TO COMPUTE THE CURVATURE AT THE CURRENT VERTEX ----------------------------------------------
		Vec3f normal = mesh.normal(it.handle());
		Vector3d N(normal[0],normal[1],normal[2]); // example of converting to Eigen's vector class for easier math

		// In the end you need to fill in this struct
		
		CurvatureInfo info;
		info.curvatures[0] = 0;
		info.curvatures[1] = 0;
		info.directions[0] = Vec3f();
		info.directions[1] = Vec3f();

		mesh.property(curvature,it) = info;
		// -------------------------------------------------------------------------------------------------------------
	}
}
