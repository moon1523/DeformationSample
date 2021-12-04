#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <functional>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#include <igl/readPLY.h>
#include <igl/readTGF.h>
#include <igl/directed_edge_parents.h>
#include <igl/directed_edge_orientations.h>
#include <igl/forward_kinematics.h>
#include <igl/deform_skeleton.h>
#include <igl/lbs_matrix.h>
#include <igl/dqs.h>
#include <igl/opengl/glfw/Viewer.h>

using namespace std;
using namespace Eigen;

typedef vector<Quaterniond, aligned_allocator<Quaterniond>>
RotationList;

int main(int argc, char** argv)
{
	double angle = 60;
	MatrixXd V, U, C, W, M;
	MatrixXi F, BE;
	VectorXi P;

	if (argc > 1)
		angle = atof(argv[1]);

	igl::readPLY("./cuboid.ply",V,F);  // V: Vertices, F: Face
	U=V;

	// Linear weights
	double minZ = V.col(2).minCoeff();
	double maxZ = V.col(2).maxCoeff();
	W = MatrixXd::Ones(V.rows(), 2);
	W.col(1) = (V.col(2).array() - minZ) * ( 1 / (maxZ - minZ));
	W.col(0) -= W.col(1);

	igl::opengl::glfw::Viewer viewer;

	MatrixXd axis(3, 9);
	axis << 0, 0, 0, 1, 0, 0, 1, 0, 0, // axisX: starting point, end point, RGB
			0, 0, 0, 0, 1, 0, 0, 1, 0, // axisY: starting point, end point, RGB
			0, 0, 0, 0, 0, 1, 0, 0, 1; // axisZ: starting point, end point, RGB
	viewer.data().lines = axis;

	viewer.data().set_mesh(U, F);
	viewer.data().show_faces = true;
	viewer.data().show_lines = true;
	viewer.data().line_width = 0.5;

	// lambda function
	auto pre_draw = [&](igl::opengl::glfw::Viewer &viewer)->bool {
		RotationList vQ;
		Quaterniond quat0 = Quaterniond::Identity();
		Quaterniond quat1(AngleAxisd(angle*M_PI/180, Vector3d::UnitY()));
		vQ.push_back(quat0);
		vQ.push_back(quat1);
		vector<Vector3d> vT;
		vT.push_back(Vector3d::Zero());
		vT.push_back(Vector3d::Zero());
		igl::dqs(V,W,vQ,vT,U);

		viewer.data().set_vertices(U);
		viewer.data().set_data(W.col(0), 0, 1, igl::COLOR_MAP_TYPE_PARULA);
		return false;
	};

	auto key_down = [&](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)->bool {

		switch(key) {
		case 'D':
		case 'd':
			cout << "'d' key down" << endl;
			return true;
		}
		return false;
	};

	viewer.callback_pre_draw = pre_draw;
	viewer.callback_key_down = key_down;

	viewer.launch();

	return EXIT_SUCCESS;
}
