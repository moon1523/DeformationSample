#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

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
using namespace igl;

typedef vector<Quaterniond, aligned_allocator<Quaterniond>>
RotationList;

const Eigen::RowVector3d sea_green(70./255.,252./255.,167./255.);
MatrixXd V, U, C, W, M;
MatrixXi F, BE;
VectorXi P;
double angle;

bool use_dqs = false;
bool recompute = true;

bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
	if (recompute)
	{
		RotationList vQ;
		Quaterniond quat0 = Quaterniond::Identity();
		Quaterniond quat1(AngleAxisd(angle*M_PI/180, Vector3d::UnitY()));
		vQ.push_back(quat0);
		vQ.push_back(quat1);
		vector<Vector3d> vT;
		vT.push_back(Vector3d::Zero());
		vT.push_back(Vector3d::Zero());
		igl::forward_kinematics(C,BE,P,vQ,vQ,vT);

		// Set transformation matrix
		const int dim = C.cols();
		MatrixXd T(BE.rows()*(dim+1),dim);
		for (int i=0; i<BE.rows(); i++)
		{
			Affine3d a = Affine3d::Identity();
			a.translate(vT[i]);
			a.rotate(vQ[i]);
			T.block(i*(dim+1), 0, dim+1, dim) = a.matrix().transpose().block(0,0,dim+1,dim);
		}

		if (use_dqs) {
			cout << "DQS" << endl;
			igl::dqs(V,W,vQ,vT,U);
		}
		else
		{
			cout << "LBS" << endl;
			igl::lbs_matrix(V,W,M);
			U = M * T;
		}
		// Also deform skeleton edges
		MatrixXd CT;
		MatrixXi BET;
		igl::deform_skeleton(C,BE,T,CT,BET);

		viewer.data().set_vertices(U);
		viewer.data().set_points(CT, sea_green);
		viewer.data().set_edges(CT, BET, sea_green);

		recompute = false;
	}

	return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
{
	recompute = true;
	switch(key)
	{
	case 'D':
	case 'd':
		use_dqs = !use_dqs;
		return true;
	}
	return false;
}

int main(int argc, char** argv)
{
	angle = 60;
	if (argc > 1)
		angle = stod(argv[1]);

	RotationList rest_pose;

	igl::readPLY("./cuboid.ply",V,F);  // V: Vertices, F: Face
	U=V;
	igl::readTGF("./cuboid.tgf",C,BE); // C: joint points, BE: BonE
	igl::directed_edge_parents(BE, P); // P: Parents bone ID
	igl::directed_edge_orientations(C, BE, rest_pose); // rest_pose: set bone rotation

	// Linear weights
	double minZ = V.col(2).minCoeff();
	double maxZ = V.col(2).maxCoeff();
	W = MatrixXd::Ones(V.rows(), 2);
	W.col(1) = (V.col(2).array() - minZ) * ( 1 / (maxZ - minZ));
	W.col(0) -= W.col(1);


	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(U, F);
	viewer.data().set_edges(C, BE, sea_green);
	viewer.data().set_points(C, sea_green);
	viewer.data().show_faces = true;
	viewer.data().show_lines = true;
	viewer.data().line_width = 0.5;
	viewer.data().point_size = 10;
//	viewer.core().trackball_angle.normalize();

	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	cout << "Press [d] to toggle between LBS and DQS" << endl;
	viewer.launch();

	return EXIT_SUCCESS;
}

























