#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <float.h>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "dual_quat_cu.hpp"

using namespace std;
using namespace Eigen;
using namespace Tbx;

void LBS(double angle, MatrixXd V, MatrixXi F, MatrixXd W);
void DQS(double angle, MatrixXd V, MatrixXi F, MatrixXd W);
void PrintPLY(string fileName, MatrixXd U, MatrixXi F, MatrixXd W);

int main(int argc, char** argv)
{
	double angle(0);
	if (argc == 2)
	{
		angle = stod(argv[1]);
	}
	ifstream ifs("./cuboid.ply");
	if (!ifs.is_open()) {
		cerr << "file is not opened" << endl; exit(1);
	}

	MatrixXd V;
	MatrixXi F;

	int vertNo, faceNo;
	string dump;
	map<double, vector<int>> zMap;

	while (ifs >> dump)
	{
		if (dump == "vertex")
		{
			ifs >> dump;
			vertNo = stoi(dump);
			V.resize(vertNo, 3);
		}
		else if (dump == "face")
		{
			ifs >> dump;
			faceNo = stoi(dump);
			F.resize(faceNo, 4);
		}
		else if (dump == "end_header")
		{
			double x,y,z;
			for (int i=0; i<vertNo; i++)
			{
				ifs >> x >> y >> z;
				V(i,0) = x;
				V(i,1) = y;
				V(i,2) = z;
				zMap[z].push_back(i);
			}

			int a,b,c;
			for (int i=0; i<faceNo; i++)
			{
				ifs >> dump >> a >> b >> c;
				F(i,0) = 3;
				F(i,1) = a;
				F(i,2) = b;
				F(i,3) = c;
			}
		}
	}
	ifs.close();

	double minZ = V.col(2).minCoeff();
	double maxZ = V.col(2).maxCoeff();

	// Linear weights
	MatrixXd W = MatrixXd::Ones(V.rows(), 2);
	W.col(1) = (V.col(2).array() - minZ) * ( 1 / (maxZ - minZ));
	W.col(0) -= W.col(1);

	// MatrixXd W;
	// W.resize(vertNo, 2);
	// //	// Linear weights
	// //	for (auto itr:zMap) {
	// //		for (auto itr2:itr.second) {
	// //			W(itr2,0) =  -1/6. * itr.first + 0.5;
	// //			W(itr2,1) =  1 - W(itr2,0);
	// //		}
	// //	}
	
	// // Quadratic weights
	// for (auto itr:zMap) {
	// 	for (auto itr2:itr.second) {
	// 		if (itr.first < 0) {
	// 			W(itr2,0) = -1/18. * (itr.first + 3) * (itr.first + 3) + 1;
	// 			W(itr2,1) = 1 - W(itr2,0);
	// 		}
	// 		else {
	// 			W(itr2,0) = 1/18. * (itr.first - 3) * (itr.first - 3);
	// 			W(itr2,1) = 1 - W(itr2,0);
	// 		}
	// 	}
	// }

	LBS(angle, V, F, W);
	DQS(angle, V, F, W);

	return 0;
}


void DQS(double angle, MatrixXd V, MatrixXi F, MatrixXd W)
{
	AngleAxisd yAA1(0, Vector3d::UnitY());
	AngleAxisd yAA2(angle * M_PI/180, Vector3d::UnitY());
	MatrixXd T1 = yAA1.matrix();
	MatrixXd T2 = yAA2.matrix();

	MatrixXd U;
	U.resize(V.rows(),3);

	// 4x4 Matrix <=> Dual Quaternion (DQ)
	// They include the rotation and translation data.
	Point3 center(0, 0, 0);
	Transfo mat1( T1(0,0), T1(0,1), T1(0,2), center.x,
			      T1(1,0), T1(1,1), T1(1,2), center.y,
				  T1(2,0), T1(2,1), T1(2,2), center.z,
				        0,       0,       0, 1  );
	Transfo mat2( T2(0,0), T2(0,1), T2(0,2), center.x,
				  T2(1,0), T2(1,1), T2(1,2), center.y,
				  T2(2,0), T2(2,1), T2(2,2), center.z,
						0,       0,       0, 1  );

	// Push the DQ data for the joints
	// Joint 1: Box bottom
	// Joint 2: Box center
	vector<Dual_quat_cu> dual_quat;
	dual_quat.push_back(Dual_quat_cu(mat1));
	dual_quat.push_back(Dual_quat_cu(mat2));

	for (int i=0; i<V.rows(); i++)
	{
		Dual_quat_cu dq_blend;
		// DQ = a + epsilon * b;
		// a = q0: quaternion which contains rotation data (non-dual part)
		// b = qe: quaternion which contains translation data (dual part)

		bool first(true);
		Quat_cu q0;

		for (int j=0; j<2; j++)
		{
			if (first)
			{
				dq_blend = dual_quat[j] * W(i,j); // translation is not considered in root joint
				q0 = dual_quat[j].rotation();
				first = false;
				continue;
			}

			if (dual_quat[j].rotation().dot( q0 ) < 0.f )
				dq_blend = dq_blend + dual_quat[j] * (-W(i,j));
			else
				dq_blend = dq_blend + dual_quat[j] * W(i,j);
		}
		Point3 v(V(i,0), V(i,1), V(i,2));
		Point3 u = dq_blend.transform(v);

		U(i,0) = u.x;
		U(i,1) = u.y;
		U(i,2) = u.z;
	}
	PrintPLY("DQS", U, F, W);
}

void LBS(double angle, MatrixXd V, MatrixXi F, MatrixXd W)
{
	MatrixXd U;
	// LBS
	AngleAxisd yAA1(0, Vector3d::UnitY());
	AngleAxisd yAA2(angle * M_PI/180, Vector3d::UnitY());
	MatrixXd T1 = yAA1.matrix();
	MatrixXd T2 = yAA2.matrix();

	U.resize(V.rows(), 3);
	U = V;

	// for (int i=0; i<U.rows(); i++) {
	// 	Vector3d u = W(i,0) * T1 * U.row(i).transpose() + W(i,1) * T2 * U.row(i).transpose();
	// 	U.row(i) = u.transpose();
	// }

	// U.col(0) = W.array().col(0) * (V * T1.transpose()).array().col(0) + W.array().col(1) * (V * T2.transpose()).array().col(0);
	// U.col(1) = W.array().col(0) * (V * T1.transpose()).array().col(1) + W.array().col(1) * (V * T2.transpose()).array().col(1);
	// U.col(2) = W.array().col(0) * (V * T1.transpose()).array().col(2) + W.array().col(1) * (V * T2.transpose()).array().col(2);

	U = (V * T1.transpose()).array().colwise() * W.col(0).array() + (V * T2.transpose()).array().colwise() * W.col(1).array();

	PrintPLY("LBS", U, F, W);
}

void PrintPLY(string fileName, MatrixXd U, MatrixXi F, MatrixXd W)
{
	ofstream ofs (fileName + ".ply");
	double vMax = W.col(0).maxCoeff();
	double vMin = W.col(0).minCoeff();

	double dv(vMax-vMin);

	MatrixXi RGB;
	RGB.resize(W.rows(), 3);
	for (int i=0; i<W.rows(); i++)
	{
		double v = W(i,0);
		int r,g,b;

		if ( v < ( vMin + 0.25 * dv) ) {
			r = 0;
			g = (255 * (4 * (v - vMin) / dv));
			b = 255;
		}
		else if ( v < vMin + 0.5 * dv) {
			r = 0;
			g = 255;
			b = (255 * (1 + 4 * (vMin + 0.25 * dv - v) / dv ));
		}
		else if ( v < ( vMin + 0.75 * dv) ) {
			r = (255 * (4 * ( v - vMin - 0.5 * dv) / dv ));
			g = 255;
			b = 0;
		}
		else {
			r = 255;
			g = (255 * (1 + 4 * (vMin + 0.75 * dv - v) / dv ));
			b = 0;
		}
		RGB(i,0) = r;
		RGB(i,1) = g;
		RGB(i,2) = b;
	}

	ofs << "ply" << endl;
	ofs << "format ascii 1.0" << endl;
	ofs << "comment Exported by RapidForm" << endl;
	ofs << "element vertex " << U.rows() << endl;
	ofs << "property float x" << endl;
	ofs << "property float y" << endl;
	ofs << "property float z" << endl;
	ofs << "property uchar red" << endl;
	ofs << "property uchar green" << endl;
	ofs << "property uchar blue" << endl;
	ofs << "element face " << F.rows() << endl;
	ofs << "property list uchar int vertex_indices" << endl;
	ofs << "end_header" << endl;

	for (int i=0;i<U.rows();i++) {
		ofs << U.row(i) << " " << RGB.row(i) << endl;
	}
	ofs << F << endl;
	ofs.close();
}














