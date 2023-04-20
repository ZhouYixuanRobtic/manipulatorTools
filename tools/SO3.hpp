/**
 * Some functions are copied from <Modern Robotics: Mechanics, Planning, and Control> Code library,
 * see https://github.com/NxRLab/ModernRobotics for more information.
 * These functions may renamed, if so will be denoted
 * Author: YX.E.Z
 * Date: Oct 1 2022
 */

#ifndef GCOPTER_SO3_HPP
#define GCOPTER_SO3_HPP

#include "Eigen/Eigen"
#include "iostream"
namespace SO3{
    constexpr double eps=10e-16;

    /**
     * # From MR Code Library with "RotInv" #
     *  Inverts a rotation matrix
     * @param R A rotation matrix
     * @return The inverse of R
     */
    inline static Eigen::Matrix3d inv(const Eigen::Matrix3d& R){
        return R.transpose();
    }

    /**
     * # From MR Code Library with "RotInv" #
     *  Inverts a rotation matrix in place
     * @param R A rotation matrix
     * @return The inverse of R
     */
    inline static void invInPlace(Eigen::Matrix3d& R){
        R.transposeInPlace();
        return;
    }

    /**
     * # From MR Code Library with "VecToso3" #
     * Converts a 3-vector to an so(3) representation
     * @param omg A 3-vector
     * @return The skew symmetric representation of omg
     */
    inline static Eigen::Matrix3d skew(const Eigen::Vector3d& omg){
        Eigen::Matrix3d R;
        R<<0, -omg[2], omg[1],
           omg[2], 0, -omg[0],
           -omg[1], omg[0], 0;
        return R;
    }

    /**
     * # From MR Code Library with "so3ToVec" #
     * Converts an so(3) representation to a 3-vector
     * @param so3mat A 3x3 skew-symmetric matrix
     * @return The 3-vector corresponding to so3mat
     */
    inline static Eigen::Vector3d vex(const Eigen::Matrix3d& so3mat){
        return Eigen::Vector3d{so3mat(2, 1), so3mat(0, 2), so3mat(1, 0)};
    }

    /**
     * # From MR Code Library with "MatrixExp3" #
     * Computes the matrix exponential of a matrix in so(3)
     * @param so3mat A 3x3 skew-symmetric matrix
     * @return The matrix exponential of so3mat
     */
    inline static Eigen::Matrix3d exp(const Eigen::Matrix3d& so3mat){
        double theta = vex(so3mat).norm();
        if(theta < eps){
            return Eigen::Matrix3d::Identity();
        }else{
            Eigen::Matrix3d omgMat = so3mat / theta;
            Eigen::Matrix3d ret; ret.setIdentity();
            ret.noalias() += std::sin(theta) * omgMat;
            ret.noalias() += (1. - std::cos(theta)) * (omgMat * omgMat);
            return ret;
        }
    }

    /**
     * # From MR Code Library with "MatrixLog3" #
     * Computes the matrix logarithm of a rotation matrix
     * @param R A 3x3 rotation matrix
     * @return The matrix logarithm of R
     */
    inline static Eigen::Matrix3d log(const Eigen::Matrix3d& R){
        double acosInput = (R.trace() - 1) / 2.;
        if(acosInput >= 1.){
            return Eigen::Matrix3d::Zero();
        }else if(acosInput <= -1.){
            Eigen::Vector3d omg;
            double constant;
            for(int i = 2; i >=0; --i){
                constant = 1 + R(i, i);
                if(constant >= eps or i==0){
                    omg = R.col(i);
                    omg[i] = constant;
                    omg /= std::sqrt(2. * constant);
                    omg *= M_PI;
                    return skew(omg);
                }
            }
        }else{
            //warning std::acos only gives one result, so take angleAxis to get angle.
            double theta = Eigen::AngleAxisd(R).angle();
            if(theta > std::numeric_limits<double>::epsilon()){
                theta /= 2. * std::sin(theta);
                return theta * (R - R.transpose());
            }    
        }
        return Eigen::Matrix3d::Zero();
    }

    /**
     * # From MR Code Library with "AxisAng3" #
     * Converts a 3-vector of exponential coordinates for rotation into axis-angle form
     * @param omg A 3-vector of exponential coordinates for rotation
     * @return <A unit rotation axis, The corresponding rotation angle>
     */
    inline static std::pair<Eigen::Vector3d, double> axisAngle(const Eigen::Vector3d& omg){
        return std::make_pair(omg.normalized(), omg.norm());
    }

    /**
     * A convenient function for directly compute a SO3 element
     * @param omega The exponential coordinates
     * @return The corresponding SO3 element
     */
    inline static Eigen::Matrix3d Exp(const Eigen::Vector3d& omega){
        return exp(skew(omega));
    }

    /**
     * A convenient function for directly compute the exponential coordinates
     * @param R The SO3 element
     * @return The corresponding exponential coordinates
     */
    inline static Eigen::Vector3d Log(const Eigen::Matrix3d& R){
        return vex(log(R));
    }

    /**
     * Returns the left Jacobian on lie group.
     * Originally, this is copied from Sophus so3.hpp \leftJacobian
     * @param omega The exponential coordinates
     * @param theta_o A precomputed `theta` can be optionally passed in
     * @return the left Jacobian on lie group
     */
    inline static Eigen::Matrix3d leftJacobian(const Eigen::Vector3d& omega, double* theta_o= nullptr){
        double const theta_sq = theta_o ? *theta_o * *theta_o : omega.squaredNorm();
        Eigen::Matrix3d const Omega = skew(omega);
        Eigen::Matrix3d V = Eigen::Matrix3d::Identity();

        if (theta_sq < eps * eps) {
            V.noalias() += 0.5 * Omega;
        } else {
            double theta = theta_o ? *theta_o : std::sqrt(theta_sq);
            V.noalias() += (1. - std::cos(theta)) / theta_sq * Omega;
            V.noalias() += (theta - std::sin(theta)) / (theta_sq * theta) * (Omega * Omega);
        }
        return V;
    }

    /**
     * Returns the left Jacobian inverse on lie group.
     * Originally, this is copied from Sophus so3.hpp \leftJacobianInverse
     * @param omega The exponential coordinates
     * @param theta_o A precomputed `theta` can be optionally passed in
     * @return the left Jacobian inverse on lie group
     */
    inline static Eigen::Matrix3d leftJacobianInverse(const Eigen::Vector3d& omega, double* theta_o= nullptr){
        double const theta_sq = theta_o ? *theta_o * *theta_o : omega.squaredNorm();
        Eigen::Matrix3d const Omega = SO3::skew(omega);

        Eigen::Matrix3d V_inv = Eigen::Matrix3d::Identity();
        V_inv.noalias() -= 0.5 * Omega;

        if (theta_sq < eps * eps) {
            V_inv.noalias() += double (1. / 12.) * (Omega * Omega);
        } else {
            double const theta = theta_o ? *theta_o : std::sqrt(theta_sq);
            double const half_theta = 0.5 * theta;
            V_inv.noalias() += (1. - half_theta * std::cos(half_theta) / std::sin(half_theta)) / (theta * theta)
                    * (Omega * Omega);
        }
        return V_inv;
    }
    /**
     * # From MR Code Library with "ProjectToSO3" #
     * Returns a projection of mat into SO(3)
     * @param mat A matrix near SO(3) to project to SO(3)
     * @return The closest matrix to R that is in SO(3)
     * Projects a matrix mat to the closest matrix in SO(3) using singular-value
     * decomposition (see
     * https://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review).
     * This function is only appropriate for matrices close to SO(3).
     */
    inline static Eigen::Matrix3d project(const Eigen::Matrix3d& mat){
        Eigen::JacobiSVD<Eigen::Matrix3d> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix3d R = svd.matrixU() * svd.matrixV().transpose();
        if (R.determinant() < 0)
            // In this case the result may be far from M; reverse sign of 3rd column
            R.col(2) *= -1;
        return R;
    }

    /**
     * # From MR Code Library with "DistanceToSO3" #
     * Returns the Frobenius norm to describe the distance of mat from the
     * SO(3) manifold
     * @param mat A 3x3 matrix
     * @return A quantity describing the distance of mat from the SO(3) manifold
     * Computes the distance from mat to the SO(3) manifold using the following method:
     * If det(mat) <= 0, return a large number.
     * If det(mat) > 0, return norm(mat^T.mat - I).
     */
    inline static double distanceToManifold(const Eigen::Matrix3d& mat){
        if(mat.determinant() > 0.){
            return (mat.transpose() * mat - Eigen::Matrix3d::Identity()).norm();
        }else{
            return 1e+9;
        }
    }

    /**
     * Returns the geodesic distance of two SO3 elements from start to the end
     * @param start a SO3 mat
     * @param end another SO3 mat
     * @return geodesic distance
     * Computes the Euclidean norm of the exponential coordinates of the interstate
     * rho(R,R') = ||log(R^T R')||_F
     */
    inline static double distance(const Eigen::Matrix3d& start, const Eigen::Matrix3d& end){
        return log(start.transpose() * end).norm();
    }

    /**
     * # From MR Code Library with "TestIfSO3" #
     * Returns true if mat is close to or on the manifold SO(3)
     * @param mat A 3x3 matrix
     * @return True if mat is very close to or in SO(3), false otherwise
     * Computes the distance d from mat to the SO(3) manifold using the
     * following method:
     * If det(mat) <= 0, d = a large number.
     * If det(mat) > 0, d = norm(mat^T.mat - I).
     * If d is close to zero, return true. Otherwise, return false.
     */
    inline static bool isSO3(const Eigen::Matrix3d& mat){
        return std::abs(distanceToManifold(mat)) < eps;
    }

    /**
     *  Compute arithmetic mean of a set of SO3 elements
     * @param list A list of all SO3 elements
     * @return The arithmetic mean of all orientation in the form of SO3
     * Computes the mean using the following method:
     * \bar{R} = 1/n * \sum_{i=1}^n R_i
     * \bar{R} = SO3(\bar{R}) # ensure the mean is a SO3 matrix
     */
    inline static Eigen::Matrix3d arithmeticMean(const std::vector<Eigen::Matrix3d>& list){
        assert(!list.empty());
        Eigen::Matrix3d pose = Eigen::Matrix3d::Zero();
        for(const auto& item : list){
            pose.noalias() += item;
        }
        pose /= list.size();
        return project(pose);
    }

    /**
     *  Compute tangent mean of a set of SO3 elements
     * @param list A list of all SO3 elements
     * @return The tangent mean of all orientation in the form of SO3
     * Computes the mean using the following method:
     * [omega_i] = logm(R_i)
     * \bar{omega} = 1/n \sum_{i=1}^n omega_i
     * \bar{R} = expm([\bar{omega}])
     */
    inline static Eigen::Matrix3d tangentMean(const std::vector<Eigen::Matrix3d>& list){
        assert(!list.empty());
        Eigen::Vector3d omg = Eigen::Vector3d::Zero();
        for(const auto& item : list){
            omg += Log(item);
        }
        omg /= list.size();
        return Exp(omg);
    }

    /**
     *  Compute geodesic mean of a set of SO3 elements
     * @param list A list of all SO3 elements
     * @return The geodesic mean of all orientation in the form of SO3
     * Computes the SO3 mean on the SO3 manifold using the following method:
     * \bar{R} = \mathop{argmin}_X \sum_{i=1}^n ||log(X^{-1} R_i||^2
     */
    inline static Eigen::Matrix3d GeodesicMean(const std::vector<Eigen::Matrix3d>& list, const double& error_thresh = 1e-3){
        assert(!list.empty());
        Eigen::Matrix3d mean = list[0];
        while(true){
            Eigen::Vector3d omg = Eigen::Vector3d::Zero();
            for(const auto& item : list){
                omg += Log(mean.transpose() * item);
            }
            omg /= list.size();
            if(omg.norm() < error_thresh){
                return mean;
            }else{
                mean *= Exp(omg);
            }
        }
    }

    /**
     * A convenient function for return a SO3 identity matrix
     * @return A 3X3 identity matrix
     */
    inline static Eigen::Matrix3d Identity(){
        return Eigen::Matrix3d::Identity();
    }
}


#endif //GCOPTER_SO3_HPP
