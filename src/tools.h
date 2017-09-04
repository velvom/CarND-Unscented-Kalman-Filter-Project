#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * Normalize angle in radians to the range [-pi, pi)
  */
  static void NormalizeAngle1(double &r);

  /**
  * Normalize angle in radians to the range [-pi, pi]
  */
  static void NormalizeAngle2(double &r);
};

#endif /* TOOLS_H_ */
