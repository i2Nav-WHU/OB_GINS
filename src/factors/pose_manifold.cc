/*
 * OB_GINS: An Optimization-Based GNSS/INS Integrated Navigation System
 *
 * Copyright (C) 2022 i2Nav Group, Wuhan University
 *
 *     Author : Hailiang Tang
 *    Contact : thl@whu.edu.cn
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "src/factors/pose_manifold.h"
#include "src/common/rotation.h"

int PoseManifold::AmbientSize() const {
    return 7;
}

int PoseManifold::TangentSize() const {
    return 6;
}

bool PoseManifold::Plus(const double *x, const double *delta, double *x_plus_delta) const {
    Eigen::Map<const Eigen::Vector3d> p_x(x);
    Eigen::Map<const Eigen::Quaterniond> q_x(x + 3);

    Eigen::Map<const Eigen::Vector3d> p_delta(delta);

    Eigen::Quaterniond q_delta = Rotation::rotvec2quaternion(Eigen::Map<const Eigen::Vector3d>(delta + 3));

    Eigen::Map<Eigen::Vector3d> p_x_plus(x_plus_delta);
    Eigen::Map<Eigen::Quaterniond> q_x_plus(x_plus_delta + 3);

    p_x_plus = p_x + p_delta;
    q_x_plus = (q_x * q_delta).normalized();

    return true;
}

bool PoseManifold::PlusJacobian(const double *x, double *jacobian) const {
    Eigen::Map<Eigen::Matrix<double, 7, 6, Eigen::RowMajor>> jaco(jacobian);

    jaco.topRows<6>().setIdentity();
    jaco.bottomRows<1>().setZero();

    return true;
}

bool PoseManifold::Minus(const double *y, const double *x, double *y_minus_x) const {
    Eigen::Map<const Eigen::Vector3d> p_y(y);
    Eigen::Map<const Eigen::Quaterniond> q_y(y + 3);

    Eigen::Map<const Eigen::Vector3d> p_x(x);
    Eigen::Map<const Eigen::Quaterniond> q_x(x + 3);

    Eigen::Map<Eigen::Vector3d> p_y_minus_x(y_minus_x);
    Eigen::Map<Eigen::Vector3d> q_y_minus_x(y_minus_x + 3);

    p_y_minus_x = p_y - p_x;
    q_y_minus_x = Rotation::quaternion2vector((q_x.inverse() * q_y).normalized());

    return true;
}

bool PoseManifold::MinusJacobian(const double *x, double *jacobian) const {
    Eigen::Map<Eigen::Matrix<double, 6, 7, Eigen::RowMajor>> jaco(jacobian);

    jaco.rightCols<6>().setIdentity();
    jaco.leftCols<1>().setZero();

    return true;
}