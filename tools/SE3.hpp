/**
 * Some functions are copied from <Modern Robotics: Mechanics, Planning, and Control> Code library,
 * see https://github.com/NxRLab/ModernRobotics for more information.
 * These functions may renamed, if so will be denoted
 * Author: YX.E.Z
 * Date: Oct 1 2022
 */

#ifndef GCOPTER_SE3_HPP
#define GCOPTER_SE3_HPP

#include "SO3.hpp"

namespace SE3{
    constexpr double eps = SO3::eps;
    typedef Eigen::Matrix4d TMat;
    typedef Eigen::Matrix<double, 6, 1> TVec;
    typedef Eigen::Matrix<double, 6, 6> TAdj;

    /**
     * # From MR Code Library with "TransInv" #
     * Inverts a homogeneous transformation matrix
     * @param T A homogeneous transformation matrix
     * @return The inverse of T
     * Uses the structure of transformation matrices to avoid taking a matrix
     * inverse, for efficiency.
     */
    inline static TMat inv(const TMat& T){
        TMat ret = T;
        ret.topLeftCorner<3, 3>().transposeInPlace();
        ret.topRightCorner<3, 1>() = -1. * ret.topLeftCorner<3, 3>() * T.topRightCorner<3, 1>();
        return ret;
    }

    /**
     * # From MR Code Library with "VecTOse3" #
     * Converts a spatial velocity vector into a 4x4 matrix in se3
     * @param V A 6-vector representing a spatial velocity
     * @return The 4x4 se3 representation of V
     */
    inline static TMat skew(const TVec& V){
        TMat ret = TMat::Zero();
        ret.topLeftCorner<3, 3>() = SO3::skew(V.head(3));
        ret.topRightCorner<3, 1>() = V.tail(3);
        return ret;
    }

    /**
     * # From MR Code Library with "se3ToVec" #
     * Converts an se3 matrix into a spatial velocity vector
     * @param se3Mat A 4x4 matrix in se3
     * @return The spatial velocity 6-vector corresponding to se3mat
     */
    inline static TVec vex(const TMat& se3Mat){
        TVec v;
        v<<se3Mat(2, 1), se3Mat(0, 2), se3Mat(1, 0),
           se3Mat(0, 3), se3Mat(1, 3), se3Mat(2, 3);
        return v;
    }

    /**
     * Computes the adjoint representation of a homogeneous transformation matrix
     *
     * @param T A homogeneous transformation matrix
     * @param inverse for the sake of compute efficiency, when true directly compute adjoint for the inverse matrix
     * @return The 6x6 adjoint representation [AdT] of T
     */
    inline static TAdj adjoint(const TMat& T, bool inverse = false){
        TAdj ret; ret.setZero();
        const auto& R = T.topLeftCorner<3, 3>();
        const auto& p = T.topRightCorner<3, 1>();
        if(inverse){
            ret.topLeftCorner<3, 3>() = R.transpose();
            ret.bottomRightCorner<3, 3>() = ret.topLeftCorner<3, 3>();
            ret.bottomLeftCorner<3, 3>() = -1. * ret.topLeftCorner<3, 3>() * SO3::skew(p);
        }else{
            ret.topLeftCorner<3, 3>() = R;
            ret.bottomRightCorner<3, 3>() = ret.topLeftCorner<3, 3>();
            ret.bottomLeftCorner<3, 3>() = SO3::skew(p) * ret.topLeftCorner<3, 3>();
        }
        return ret;
    }

    /**
     * # From MR Code Library with "AxisAng6" #
     * Converts a 6-vector of exponential coordinates into screw axis-angle form
     * @param expc6 A 6-vector of exponential coordinates for rigid-body motion S*theta
     * @return <The corresponding normalized screw axis, The distance traveled along/about S>
     */
    inline static std::pair<TVec, double> axisAngle(const TVec& expc6){
        double theta = expc6.head(3).norm();
        if(theta < eps){
            theta = expc6.tail(3).norm();
        }
        return theta > eps ? std::make_pair(TVec{expc6 / theta}, theta) : std::make_pair(TVec::Zero(), 0.);
    }

    /**
     * # From MR Code Library with "MatrixExp6" #
     * Computes the matrix exponential of an se3 representation of exponential coordinates
     * @param se3Mat A matrix in se3
     * @return The matrix exponential of se3mat
     */
    inline static TMat exp(const TMat& se3Mat){
        const auto& so3 = se3Mat.topLeftCorner<3, 3>();
        const auto& t = se3Mat.topRightCorner<3, 1>();
        const auto omega = SO3::vex(so3);

        TMat ret = TMat::Identity();
        ret.topLeftCorner<3, 3>() = SO3::exp(so3);
        ret.topRightCorner<3, 1>() = SO3::leftJacobian(omega) * t;
        return ret;
    }

    /**
     * # From MR Code Library with "MatrixLog6" #
     * Computes the matrix logarithm of a homogeneous transformation matrix
     * @param T A matrix in SE3
     * @return The matrix logarithm of R
     */
    inline static TMat log(const TMat& T){
        const auto& R = T.topLeftCorner<3, 3>();
        const auto& t = T.topRightCorner<3, 1>();
        const Eigen::Matrix3d omgMat = SO3::log(R);
        TMat ret = TMat::Zero();
        ret.topLeftCorner<3, 3>() = omgMat;
        ret.topRightCorner<3, 1>() = SO3::leftJacobianInverse(SO3::vex(omgMat)) * t;
        return ret;
    }

    /**
     * A convenient function for directly compute a SE3 element
     * @param V The exponential coordinates
     * @return The corresponding SE3 element
     */
    inline static TMat Exp(const TVec & V){
        return exp(skew(V));
    }

    /**
     * A convenient function for directly compute the exponential coordinates
     * @param T The SE3 element
     * @return The corresponding exponential coordinates
     */
    inline static TVec Log(const TMat & T){
        return {vex(log(T))};
    }

    /**
     * # From MR Code Library with "ProjectToSE3" #
     * Returns a projection of mat into SE(3)
     * @param mat A 4x4 matrix to project to SE(3)
     * @return The closest matrix to T that is in SE(3)
     * Projects a matrix mat to the closest matrix in SE(3) using singular-value
     * decomposition (see
     * http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review).
     * This function is only appropriate for matrices close to SE(3).
     */
    inline static TMat project(const TMat& mat){
        TMat ret = mat;
        ret.topLeftCorner<3, 3>() = SO3::project(mat.topLeftCorner<3, 3>());
        return ret;
    }

    /**
     * # From MR Code Library with "DistanceToSE3" #
     * Returns the Frobenius norm to describe the distance of mat from the
     * SE(3) manifold
     * @param mat A 4x4 matrix
     * @return A quantity describing the distance of mat from the SE(3) manifold
     * Computes the distance from mat to the SE(3) manifold using the following
     * method:
     * Compute the determinant of matR, the top 3x3 submatrix of mat.
     * If det(matR) <= 0, return a large number.
     * If det(matR) > 0, replace the top 3x3 submatrix of mat with matR^T.matR,
     * and set the first three entries of the fourth column of mat to zero. Then
     * return norm(mat - I).
     */
    inline static double distanceToManifold(const TMat& mat){
        const auto& R = mat.topLeftCorner<3, 3>();
        if(R.determinant() > 0.){
            TMat ret = mat;
            ret.topLeftCorner<3, 3>() = R.transpose() * R;
            ret.topRightCorner<3, 1>().setZero();
            ret -= TMat::Identity();
            return ret.norm();
        }else{
            return 1e+9;
        }
    }

    /**
     * Returns the geodesic distance of two SE3 elements from start to the end
     * @param start a SE3 mat
     * @param end another SE3 mat
     * @return geodesic distance
     * Computes the Euclidean norm of the exponential coordinates of the interstate
     * rho(T,T') = ||log(R^{-1} R')||_F
     */
    inline static double distance(const TMat& start, const TMat& end){
        return log(inv(start) * end).norm();
    }

    /**
     * # From MR Code Library with "TestIfSE3" #
     * Returns true if mat is close to or on the manifold SE(3)
     * @param mat A 4x4 matrix
     * @return True if mat is very close to or in SE(3), false otherwise
     * Computes the distance d from mat to the SE(3) manifold using the
     * following method:
     * Compute the determinant of the top 3x3 submatrix of mat.
     * If det(mat) <= 0, d = a large number.
     * If det(mat) > 0, replace the top 3x3 submatrix of mat with mat^T.mat, and
     * set the first three entries of the fourth column of mat to zero.
     * Then d = norm(T - I).
     * If d is close to zero, return true. Otherwise, return false.
     */
    inline static bool isSE3(const TMat& mat){
        return std::abs(distanceToManifold(mat)<eps);
    }

    /**
     *  Compute matrix converting jacobian matrix w.r.t. B frame into A frame according to pose matrix T
     * @param T Tab, pose of B frame w.r.t. A frame
     * @param inverse for computing efficiency take inverse as condition.
     * @return matrix converts jacobian
     */
    inline static TAdj changeJacobianMatrix(const TMat &T, bool inverse = false){
        TAdj ret = TAdj::Zero();
        const auto& R = T.topLeftCorner<3, 3>();
        if(inverse){
            ret.topLeftCorner<3, 3>() = R.transpose();
        }else{
            ret.topLeftCorner<3, 3>() = R;
        }
        ret.bottomRightCorner<3, 3>() = ret.topLeftCorner<3, 3>();
        return ret;
    }


    /**
     *  Compute arithmetic mean of a set of SE3 elements
     * @param list A list of all SE3 elements
     * @return The arithmetic mean of all orientation in the form of SE3
     * Computes the mean using the following method:
     * \bar{T} = 1/n * \sum_{i=1}^n T_i
     * \bar{T} = SE3(\bar{T}) # ensure the mean is a SE3 matrix
     */
    inline static TMat arithmeticMean(const std::vector<TMat>& list){
        assert(!list.empty());
        TMat pose = TMat ::Zero();
        for(const auto& item : list){
            pose.noalias() += item;
        }
        pose /= list.size();
        return project(pose);
    }

    /**
     *  Compute tangent mean of a set of SE3 elements
     * @param list A list of all SE3 elements
     * @return The tangent mean of all orientation in the form of SE3
     * Computes the mean using the following method:
     * [V_i] = logm(T_i)
     * \bar{V} = 1/n \sum_{i=1}^n V_i
     * \bar{T} = expm([\bar{V}])
     */
    inline static TMat tangentMean(const std::vector<TMat>& list){
        assert(!list.empty());
        TVec omg = TVec::Zero();
        for(const auto& item : list){
            omg += Log(item);
        }
        omg /= list.size();
        return Exp(omg);
    }

    /**
     *  Compute geodesic mean of a set of SE3 elements
     * @param list A list of all SE3 elements
     * @return The geodesic mean of all orientation in the form of SE3
     * Computes the SE3 mean on the SE3 manifold using the following method:
     * \bar{T} = \mathop{argmin}_X \sum_{i=1}^n ||log(X^{-1} T_i||^2
     */
    inline static TMat GeodesicMean(const std::vector<TMat>& list, const double& error_thresh = 1e-3){
        assert(!list.empty());
        TMat mean = list[0];
        while(true){
            TVec omg = TVec::Zero();
            for(const auto& item : list){
                omg += Log(inv(mean) * item);
            }
            omg /= list.size();
            if(omg.norm() < error_thresh){
                return mean;
            }else{
                mean *= Exp(omg);
            }
        }
    }
}

#endif //GCOPTER_SE3_HPP
