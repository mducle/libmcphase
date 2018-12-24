/* eigen.hpp
 * 
 * (C) 2018 Duc Le - duc.le@stfc.ac.uk
 * This program is licensed under the GNU General Purpose License, version 3. Please see the LICENSE file
 */

#ifndef EIGEN_H
#define EIGEN_H

#include <Eigen/Dense>

// Use row major matrices for compatibility with Numpy
using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

#endif
