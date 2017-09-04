#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

// -----------------------------------------------------------------------------
// Calculate the RMSE
// -----------------------------------------------------------------------------
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // The estimation vector size should not be zero
  // The estimation vector size should equal ground truth vector size
  if (estimations.size() == 0 || (estimations.size() != ground_truth.size())) {
    cout << "estimation vector size is 0 or not matching ground_truth vector size\n" << endl;
    return rmse;
  }

  // Accumulate squared residuals
  for (int i=0; i < estimations.size(); ++i) {
    // Residual
    VectorXd res = estimations[i] - ground_truth[i];

    // Multiply element-wise and accumulate
    //res = res.array() * res.array();
    //rmse += res;
    rmse = rmse.array() + res.array() * res.array();
  }

  // Calculate the mean
  rmse = rmse / estimations.size();

  // Calculate the squared root to get RMSE
  rmse = rmse.array().sqrt();

  return rmse;
}

// -----------------------------------------------------------------------------
// Normalize angle in radians to the range [-pi, pi)
// -----------------------------------------------------------------------------
const float M_2PI = 2 * M_PI;
void Tools::NormalizeAngle1(double &r) {
  r = fmod(r + M_PI, M_2PI);
  if (r < 0)
    r += M_2PI;
  r -= M_PI;
}

// -----------------------------------------------------------------------------
// Normalize angle in radians to the range [-pi, pi]
// -----------------------------------------------------------------------------
void Tools::NormalizeAngle2(double &r) {
  float y = fmod(r + M_PI, M_2PI);
  if (y < 0)
    y += M_2PI;

  if (r > 0.0 && y == 0.0)
     r = y + M_PI;
  else
     r = y - M_PI;
}
