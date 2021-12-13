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

MatrixXd GenerateBonePoints(MatrixXd C, MatrixXi BE, double interval)
{
  MatrixXd output = C;
  for (int i = 0; i < BE.rows(); i++)
  {
    Vector3d c0 = C.row(BE(i, 0)).transpose();
    Vector3d c1 = C.row(BE(i, 1)).transpose();
    Vector3d dir = (c1 - c0).normalized();
    double l = (c1 - c0).norm();
    int num = floor(l / interval);
    double interval1 = l / (double)num;
    int prevSize = output.rows();
    output.conservativeResize(prevSize + num - 1, NoChange);
    for (double intvl = 1; intvl < num; intvl++)
      output.row(prevSize + intvl - 1) = (c0 + dir * interval1 * intvl).transpose();
  }
  return output;
}

int main(int argc, char *argv[])
{
	MatrixXd V;
	MatrixXi F;
	igl::readPLY("outer.ply", V, F);

  // option 1>>
  //- modify code to obtain C and BE matrix by reading proper tgf file.
  //- generate tgf file by your own.
	MatrixXd C;
	MatrixXi BE;
	VectorXi P;
	igl::readTGF("patient.tgf",C,BE);

  // option 2>>
  // extend C to include all necessary joint points, and generate BE matrix manually.

//	C.resize(17,3);
//	C << 0, 4.48, 0,
//		-11, 4.48, 0,
//		8.8, 4.5, 0,
//		0, 37, 0,
//		0, 66.6, 0,
//		8.14, 57.61, 0,
//		-8.6, 57.6, 0,
//		-19, 56, 0,
//		19, 56, 0,
//		21, 27, -6.7,
//		-23, 27, -6.7,
//		25, 0.6, -6.7,
//		-25, 0.6, -6.7,
//		8.8,-37,0,
//		-11,-37,0,
//		6.7,-76,-5.6,
//		-7.6,-76,-5.6;
//
//
//	BE.resize(16, 2);
//	BE << 0, 1,
//		0, 2,
//		0, 3,
//		3, 4,
//		3, 5,
//		3, 6,
//		6, 7,
//		5, 8,
//		8, 9,   // 8  left upper A
//		7, 10,  // 9  right upper A
//		9, 11,  // 10 left lower A
//		10, 12, // 11 right lower A
//		2, 13, // 12 left thigh
//		1, 14, // 13 right thigh
//		13, 15, // 14 left calf
//		14, 16; // 15 right calf

  igl::directed_edge_parents(BE, P);

  MatrixXd boneV = GenerateBonePoints(C, BE, 2);

  V.conservativeResize(V.rows() + boneV.rows(), 3);
  V.bottomRows(boneV.rows()) = boneV;

  MatrixXd VT;
  MatrixXi TT, FT;
  igl::copyleft::tetgen::tetrahedralize(V, F, "pYq", VT, TT, FT);

  MatrixXd U = V;

  // List of boundary indices (aka fixed value indices into VV)
  VectorXi b;
  // List of boundary conditions of each weight function
  MatrixXd bc;
  igl::boundary_conditions(VT, TT, C, VectorXi(), BE, MatrixXi(), b, bc);

  // compute BBW weights matrix
  igl::BBWData bbw_data;
  // only a few iterations for sake of demo
  bbw_data.active_set_params.max_iter = 8;
  bbw_data.verbosity = 2;
  MatrixXd W;
  if (!igl::bbw(VT, TT, b, bc, bbw_data, W))
  {
    return EXIT_FAILURE;
  }

  igl::normalize_row_sums(W, W);

  // Plot the mesh with pseudocolors
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(U, F);
  viewer.data().set_edges(C, BE, RowVector3d(70. / 255., 252. / 255., 167. / 255.));
  viewer.data().show_lines = false;
  viewer.data().show_overlay_depth = false;
  viewer.data().line_width = 1;

  double deg(0), deg1(0);
  int wgt(0);
  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &_viewer, unsigned char key, int mods) -> bool
  {
    switch (key)
    {
    case '.':
      deg ++;
      break;
    case ',':
      deg --;
      break;
    case ']':
      deg1 ++;
      break;
    case '[':
      deg1 --;
      break;
    case 'W':
	  wgt ++;
	  if (wgt==BE.rows()) wgt = 0;
	  break;
    case 'Q':
	  wgt --;
	  if (wgt<0) wgt = BE.rows()-1;
	  break;
    }

    RotationList vQ(BE.rows(), Quaterniond::Identity());
    vQ[8] = AngleAxisd(deg/180.*M_PI, Vector3d::UnitZ());
    vQ[9] = AngleAxisd(deg/180.*M_PI, -Vector3d::UnitZ());
    vQ[10] = AngleAxisd(deg1/180.*M_PI, Vector3d::UnitZ());
    vQ[11] = AngleAxisd(deg1/180.*M_PI, -Vector3d::UnitZ());
    vQ[12] = AngleAxisd(deg/180.*M_PI, Vector3d::UnitZ());
    vQ[13] = AngleAxisd(deg/180.*M_PI, -Vector3d::UnitZ());
    vQ[14] = AngleAxisd(deg1/180.*M_PI, Vector3d::UnitZ());
    vQ[15] = AngleAxisd(deg1/180.*M_PI, -Vector3d::UnitZ());

    vector<Vector3d> vT;
    // propagate relative rotations via FK to retrieve absolute transformations
    igl::forward_kinematics(C, BE, P, vQ, vQ, vT);
    igl::dqs(V, W.topRows(V.rows()), vQ, vT, U);
    _viewer.data().set_vertices(U);

    MatrixXd C1 = C;
    for(int i=0;i<BE.rows();i++)
    {
      C1.row(BE(i, 1)) = vQ[i] * C.row(BE(i, 1)) + vT[i];
    }


    _viewer.data().set_points(C1, RowVector3d(70. / 255., 252. / 255., 167. / 255.));
    _viewer.data().set_edges(C1, BE, RowVector3d(70. / 255., 252. / 255., 167. / 255.));

    _viewer.data().set_data(W.topRows(V.rows()).col(wgt));
	int v1 = viewer.data_list[0].id;
	_viewer.data(v1).point_size = 8;
    return true;
  };


  viewer.core().is_animating = false;
  cout << "Press '.' to show next weight function." << endl
       << "Press ',' to show previous weight function." << endl
       << "Press [space] to toggle animation." << endl;
  viewer.launch();
  return EXIT_SUCCESS;
}

// Inputs: C, BE, P, dQ, dT
// Outputs: vQ, vT

//    IGL_INLINE void igl::forward_kinematics(
//      const Eigen::MatrixXd & C,
//      const Eigen::MatrixXi & BE,
//      const Eigen::VectorXi & P,
//      const std::vector<
//        Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > & dQ,
//      const std::vector<Eigen::Vector3d> & dT,
//      std::vector<
//        Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > & vQ,
//      std::vector<Eigen::Vector3d> & vT)
//    {
//      using namespace std;
//      using namespace Eigen;
//      const int m = BE.rows();
//      assert(m == P.rows());
//      assert(m == (int)dQ.size());
//      assert(m == (int)dT.size());
//      vector<bool> computed(m,false);
//      vQ.resize(m);
//      vT.resize(m);
//      // Dynamic programming
//      function<void (int) > fk_helper = [&] (int b)
//      {
//        if(!computed[b])
//        {
//          if(P(b) < 0)
//          {
//            // base case for roots
//            vQ[b] = dQ[b];
//            const Vector3d r = C.row(BE(b,0)).transpose();
//            vT[b] = r-dQ[b]*r + dT[b];
//          }else
//          {
//            // Otherwise first compute parent's
//            const int p = P(b);
//            fk_helper(p);
//            vQ[b] = vQ[p] * dQ[b];
//            const Vector3d r = C.row(BE(b,0)).transpose();
//            vT[b] = vT[p] - vQ[b]*r + vQ[p]*(r + dT[b]);
//          }
//          computed[b] = true;
//        }
//      };
//      for(int b = 0;b<m;b++)
//      {
//        fk_helper(b);
//      }
//    }
//
//    IGL_INLINE void igl::forward_kinematics(
//      const Eigen::MatrixXd & C,
//      const Eigen::MatrixXi & BE,
//      const Eigen::VectorXi & P,
//      const std::vector<
//        Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > & dQ,
//      std::vector<
//        Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> > & vQ,
//      std::vector<Eigen::Vector3d> & vT)
//    {
//      std::vector<Eigen::Vector3d> dT(BE.rows(),Eigen::Vector3d(0,0,0));
//      return forward_kinematics(C,BE,P,dQ,dT,vQ,vT);
//    }
