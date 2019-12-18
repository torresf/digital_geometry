
#include "Find2DAffineTransform.hpp"

// This code is released in public domain

// Given two sets of 3D points, find the rotation + translation + scale
// which best maps the first set to the second.
// Source: http://en.wikipedia.org/wiki/Kabsch_algorithm

//https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp

// The input 2D points are stored as columns.
Eigen::Affine2d Find2DAffineTransform(const Eigen::Matrix2Xd &_in, const Eigen::Matrix2Xd &_out)
{

  // Default output
  Eigen::Affine2d A;
  A.linear() = Eigen::Matrix2d::Identity(2, 2);
  A.translation() = Eigen::Vector2d::Zero();

  Eigen::Matrix2Xd in = _in;
  Eigen::Matrix2Xd out = _out;

  if (in.cols() > out.cols())
  {
    in.resize(2, out.cols());
  }
  else
  {
    out.resize(2, in.cols());
  }

  std::cout << in.cols() << " " << out.cols() << std::endl;

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.
  double dist_in = 0, dist_out = 0;
  for (int col = 0; col < in.cols() - 1; col++)
  {
    dist_in += (in.col(col + 1) - in.col(col)).norm();
    dist_out += (out.col(col + 1) - out.col(col)).norm();
  }
  if (dist_in <= 0 || dist_out <= 0)
    return A;
  double scale = dist_out / dist_in;
  out /= scale;

  // Find the centroids then shift to the origin
  Eigen::Vector2d in_ctr = Eigen::Vector2d::Zero();
  Eigen::Vector2d out_ctr = Eigen::Vector2d::Zero();
  for (int col = 0; col < in.cols(); col++)
  {
    in_ctr += in.col(col);
    out_ctr += out.col(col);
  }
  in_ctr /= in.cols();
  out_ctr /= out.cols();
  for (int col = 0; col < in.cols(); col++)
  {
    in.col(col) -= in_ctr;
    out.col(col) -= out_ctr;
  }

  // SVD
  Eigen::MatrixXd Cov = in * out.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Find the rotation
  double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
  if (d > 0)
    d = 1.0;
  else
    d = -1.0;
  Eigen::Matrix2d I = Eigen::Matrix2d::Identity(2, 2);
  I(1, 1) = d;
  Eigen::Matrix2d R = svd.matrixV() * I * svd.matrixU().transpose();

  // The final transform
  A.linear() = scale * R;
  A.translation() = scale * (out_ctr - R * in_ctr);

  return A;
}
