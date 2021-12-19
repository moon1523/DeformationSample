#include <igl/boundary_conditions.h>
#include <igl/directed_edge_parents.h>
#include <igl/deform_skeleton.h>
#include <igl/normalize_row_sums.h>
#include <igl/readPLY.h>
#include <igl/readTGF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/forward_kinematics.h>
#include <igl/bbw.h>
#include <igl/dqs.h>

#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <vector>
#include <algorithm>
#include <iostream>

typedef std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond>>
    RotationList;

using namespace Eigen;
using namespace std;

int main(int argc, char *argv[])
{
	MatrixXd V, C, U;
	MatrixXi F, BE;
	VectorXi P;

	igl::readPLY("outer.ply", V, F);
	igl::readTGF("patient.tgf",C,BE);
	igl::directed_edge_parents(BE, P);

	V.conservativeResize(V.rows() + C.rows(), 3);
	V.bottomRows(C.rows()) = C;

	MatrixXd VT;
	MatrixXi TT, FT;
	igl::copyleft::tetgen::tetrahedralize(V, F, "pYq", VT, TT, FT);

	VectorXi b; b.resize(C.rows());
	MatrixXd bc = MatrixXd::Identity(C.rows(),C.rows());
	for (int i=0;i<C.rows();i++) {
		b(i) = V.rows() - C.rows() + i;
	}

	// compute BBW weights matrix
	igl::BBWData bbw_data;
	// only a few iterations for sake of demo
	bbw_data.active_set_params.max_iter = 8;
	bbw_data.verbosity = 2;
	MatrixXd W;
	if (!igl::bbw(VT, TT, b, bc, bbw_data, W))
		return EXIT_FAILURE;
	igl::normalize_row_sums(W, W);


	MatrixXd CT, jointTrans;
	CT.resize(C.rows(), 3);
	CT = C;
	jointTrans.resize(C.rows(),3);

	double scalingFactor = -0.1;
	for (int i=0;i<BE.rows();i++) {
		if ( i == 12 || i == 13 || i == 14 || i == 15 ) { // legs
			Vector3d c0 = C.row(BE(i,0)).transpose();
			Vector3d c1 = C.row(BE(i,1)).transpose();
			Vector3d dir = (c0-c1).normalized();
			double length = (c0-c1).norm();
			CT.row(BE(i,1)) = (c1 + dir*length*scalingFactor).transpose();
		}
		else if ( i == 8 || i == 9 || i == 10 || i == 11 ) { // arm
			Vector3d c0 = C.row(BE(i,0)).transpose();
			Vector3d c1 = C.row(BE(i,1)).transpose();
			Vector3d dir = (c0-c1).normalized();
			double length = (c0-c1).norm();
			CT.row(BE(i,1)) = (c1 - dir*length*scalingFactor).transpose();
		}
	}

	jointTrans = CT-C;
	U = VT + W*jointTrans;

	// Plot the mesh with pseudocolors
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(U, F);
	viewer.data().set_points(C, RowVector3d(70. / 255., 252. / 255., 167. / 255.));
	viewer.data().set_edges(C, BE, RowVector3d(70. / 255., 252. / 255., 167. / 255.));
	viewer.data().show_lines = false;
	viewer.data().show_overlay_depth = false;
	viewer.data().line_width = 1;

	int wgt(0);
	viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &_viewer, unsigned char key, int mods) -> bool
	{

		switch (key)
		{
		case 'W':
		  wgt ++;
		  if (wgt==BE.rows()+1) wgt = 0;
		  break;
		case 'Q':
		  wgt --;
		  if (wgt<0) wgt = BE.rows()-1;
		  break;
		case 'A':
			U = VT;
			_viewer.data().set_points(C, RowVector3d(70. / 255., 252. / 255., 167. / 255.));
			_viewer.data().set_edges(C, BE, RowVector3d(70. / 255., 252. / 255., 167. / 255.));
			break;
		case 'S':
			U = VT + W*jointTrans;
			_viewer.data().set_points(CT, RowVector3d(70. / 255., 252. / 255., 167. / 255.));
			_viewer.data().set_edges(CT, BE, RowVector3d(70. / 255., 252. / 255., 167. / 255.));

			break;
		}
		_viewer.data().set_vertices(U);
		_viewer.data().set_data(W.topRows(U.rows()).col(wgt));

		return true;
	};
	int v1 = viewer.data_list[0].id;
	viewer.data(v1).point_size = 8;
	viewer.core().is_animating = false;
	viewer.launch();

	return EXIT_SUCCESS;
}







