/**
 * Some functions are copied from <Modern Robotics: Mechanics, Planning, and Control> Code library,
 * see https://github.com/NxRLab/ModernRobotics for more information.
 * These functions may renamed, if so will be denoted
 * Common Interface of screw-based robot functions
 * Author: YX.E.Z
 * Date: 2.12.2022
 */

#ifndef GCOPTER_ROBOT_HPP
#define GCOPTER_ROBOT_HPP

#include "SE3.hpp"

namespace robot{
    using SE3::TMat;
    typedef Eigen::Matrix<double, 6, -1> ScrewList;
    typedef Eigen::Matrix<double, 6, -1> Jacobian;
    typedef Eigen::Matrix<double, -1, 1> ThetaList;
    typedef TMat (*lfk_t)(const ThetaList& thetaList);
    typedef std::function<TMat(const ThetaList& thetaList)> lfk_func;


    /**
     * # From Mr library with "FKinSpace"
     * Computes forward kinematics in the space frame for an open chain robot
     * @param M The home configuration (position and orientation) of the end-effector
     * @param SList The joint screw axes in the space frame when the manipulator is at the home position,
     * in the format of a matrix with axes as the columns
     * @param thetaList A list of joint coordinates
     * @return A homogeneous transformation matrix representing the end-
     *        effector frame when the joints are at the specified coordinates
     *        (i.t.o Space Frame)
     */
    inline static TMat fkInSpace(const TMat& M, const ScrewList& SList, const ThetaList& thetaList){
        assert(thetaList.size() == SList.cols());
        TMat ret = M;
        for(int i = thetaList.size() - 1; i >= 0; --i){
            ret = SE3::Exp(SList.col(i) * thetaList[i]) * ret;
        }
        return ret;
    }

    /**
     * # From MR library with "FKInBody"
     * Computes forward kinematics in the body frame for an open chain robot
     * @param M The home configuration (position and orientation) of the end-
              effector
     * @param BList The joint screw axes in the end-effector frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
     * @param thetaList  A list of joint coordinates
     * @return  A homogeneous transformation matrix representing the end-
             effector frame when the joints are at the specified coordinates
             (i.t.o Body Frame)
     */
    inline static TMat fkInBody(const TMat& M, const ScrewList& BList, const ThetaList& thetaList){
        assert(thetaList.size() == BList.cols());
        TMat ret = M;
        for(int i = 0; i < thetaList.size(); ++i){
            ret *= SE3::Exp(BList.col(i) * thetaList[i]);
        }
        return ret;
    }

    /**
     * # From MR library
     * Computes the body Jacobian for an open chain robot
     * @param BList The joint screw axes in the end-effector frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
     * @param thetaList A list of joint coordinates
     * @return The body Jacobian corresponding to the inputs (6xn real
             numbers)
     */
    inline static Jacobian jacobianBody(const ScrewList& BList, const ThetaList& thetaList){
        assert(thetaList.size() == BList.cols());
        Jacobian ret(6, BList.cols());
        ret.rightCols(1) = BList.rightCols(1);
        TMat T = TMat::Identity();
        for(int i = thetaList.size() - 2; i > -1; --i){
            T *= SE3::Exp(BList.col(i + 1) * -thetaList[i + 1]);
            ret.col(i) = SE3::adjoint(T, false) * BList.col(i);
        }
        return ret;
    }

    /**
     * # From MR library
     * Computes the space Jacobian for an open chain robot
     * @param SList The joint screw axes in the space frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
     * @param thetaList A list of joint coordinates
     * @return The space Jacobian corresponding to the inputs (6xn real
             numbers)
     */
    inline static Jacobian jacobianSpace(const ScrewList& SList, const ThetaList& thetaList){
        assert(thetaList.size() == SList.cols());
        Jacobian ret(6, SList.cols());
        ret.leftCols(1) = SList.leftCols(1);
        TMat T = TMat::Identity();
        for(int i = 1; i < SList.cols(); ++i){
            T *= SE3::Exp(SList.col(i - 1) * thetaList[i - 1]);
            ret.col(i) = SE3::adjoint(T, false) * SList.col(i);
        }
        return ret;
    }

    /**
     * Convert a MR-form space jacobian to a norm form jacobian
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param space_jacobian jacobian in MR form
     * @param in_ee True for change into the end frame
     * @return Norm form jacobian
     */
    inline static Jacobian fromMRSpaceJacobian(const TMat& T, const Jacobian& space_jacobian, bool in_ee=false){
        Jacobian ret = SE3::adjoint(T, true) * space_jacobian;
        ret.topRows(3).swap(ret.bottomRows(3));
        if(!in_ee){
            ret = SE3::changeJacobianMatrix(T, false) * ret;
        }
        return ret;
    }

    /**
     * Inplace version of \fromMRSpaceJacobian
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param space_jacobian jacobian in MR form
     * @param in_ee True for change into the end frame
     */
    inline static void fromMRSpaceJacobianInPlace(const TMat& T, Jacobian& space_jacobian, bool in_ee=false){
        space_jacobian = SE3::adjoint(T, true) * space_jacobian;
        space_jacobian.topRows(3).swap(space_jacobian.bottomRows(3));
        if(!in_ee){
            space_jacobian = SE3::changeJacobianMatrix(T, false) * space_jacobian;
        }
    }

    /**
     * Convert a MR-form space jacobian to a norm form jacobian
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param body_jacobian jacobian in MR form
     * @param in_ee True for change into the end frame
     * @return Norm form jacobian
     */
    inline static Jacobian fromMRBodyJacobian(const TMat& T, const Jacobian& body_jacobian, bool in_ee=true){
        Jacobian ret = body_jacobian;
        ret.topRows(3).swap(ret.bottomRows(3));
        if(!in_ee){
            ret = SE3::changeJacobianMatrix(T, false) * ret;
        }
        return ret;
    }

    /**
     * Inplace version of \fromMRBodyJacobian
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param space_jacobian jacobian in MR form
     * @param in_ee True for change into the end frame
     */
    inline static void fromMRBodyJacobianInPlace(const TMat& T, Jacobian& body_jacobian, bool in_ee=true){
        body_jacobian.topRows(3).swap(body_jacobian.bottomRows(3));
        if(!in_ee){
            body_jacobian = SE3::changeJacobianMatrix(T, false) * body_jacobian;
        }
    }

    /**
     * Convert a Norm-form space jacobian to a MR form jacobian
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param space_jacobian jacobian in MR form
     * @param in_ee True for change into the end frame
     * @return MR form jacobian
     */
    inline static Jacobian fromNormSpaceJacobian(const TMat& T, const Jacobian& space_jacobian, bool in_ee=false){
        Jacobian ret = SE3::changeJacobianMatrix(T, true) * space_jacobian;
        ret.topRows(3).swap(ret.bottomRows(3));
        if(!in_ee){
            ret = SE3::adjoint(T, false) * ret;
        }
        return ret;
    }

    /**
     * Inplace version of \fromNormSpaceJacobian
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param space_jacobian jacobian in MR form
     * @param in_ee True for change into the end frame
     * @return MR form jacobian
     */
    inline static void fromNormSpaceJacobianInPlace(const TMat& T, Jacobian& space_jacobian, bool in_ee=false){
        space_jacobian = SE3::changeJacobianMatrix(T, true) * space_jacobian;
        space_jacobian.topRows(3).swap(space_jacobian.bottomRows(3));
        if(!in_ee){
            space_jacobian = SE3::adjoint(T, false) * space_jacobian;
        }
    }

    /**
     * Convert a Norm-form space jacobian to a MR form jacobian
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param body_jacobian jacobian in MR form
     * @param in_ee True for change into the end frame
     * @return MR form jacobian
     */
    inline static Jacobian fromNormBodyJacobian(const TMat& T, const Jacobian& body_jacobian, bool in_ee= false){
        Jacobian ret = body_jacobian;
        ret.topRows(3).swap(ret.bottomRows(3));
        if(!in_ee){
            ret = SE3::adjoint(T, false) * ret;
        }
        return ret;
    }

    /**
     * Inplace version of \fromNormBodyJacobian
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param body_jacobian jacobian in MR form
     * @param in_ee True for change into the end frame
    */
    inline static void fromNormBodyJacobianInPlace(const TMat& T, Jacobian& body_jacobian, bool in_ee= false){
        body_jacobian.topRows(3).swap(body_jacobian.bottomRows(3));
        if(!in_ee){
            body_jacobian = SE3::adjoint(T, false) * body_jacobian;
        }
    }

    /**
     * Compute Screw list through fk function
     * @param fk_function function handle of forward kinematics
     * @param SList Output Screw list
     * @param body True for computing screw w.r.t. body
     * @param joint_numbers the total joint numbers
     * @return Initial condition
     */
    inline static TMat computeScrewList(lfk_t fk_function, ScrewList& SList, bool body=false, const int joint_numbers=7){
        SList.resize(6, joint_numbers);
        ThetaList q(joint_numbers); q.setZero();
        TMat M = SE3::inv(fk_function(q)); TMat Ti;
        for(int i = 0; i < joint_numbers; ++i){
            q.setZero(); q[i] = 1.;
            Ti = fk_function(q);
            SList.col(i) = SE3::Log(Ti * M);
        }
        M = SE3::inv(M);
        if(body){
            SList = SE3::adjoint(M, true) * SList;
        }
        return M;
    }

    /**
     * Compute Screw list through fk function
     * @param fk_function function handle of forward kinematics
     * @param SList Output Screw list
     * @param body True for computing screw w.r.t. body
     * @param joint_numbers the total joint numbers
     * @return Initial condition
     */
    inline static TMat computeScrewList(lfk_func fk_function, ScrewList& SList, bool body=false, const int joint_numbers=7){
        SList.resize(6, joint_numbers);
        ThetaList q(joint_numbers); q.setZero();
        TMat M = SE3::inv(fk_function(q)); TMat Ti;
        for(int i = 0; i < joint_numbers; ++i){
            q.setZero(); q[i] = 1.;
            Ti = fk_function(q);
            SList.col(i) = SE3::Log(Ti * M);
        }
        M = SE3::inv(M);
        if(body){
            SList = SE3::adjoint(M, true) * SList;
        }
        return M;
    }

    /**
     * Compute analytical jacobian d([r, x])/d(q)
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param BList The joint screw axes in the end-effector frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
     * @param thetaList  A list of joint coordinates
     * @return analytical jacobian d([r, x])/d(q)
     */
    inline static Jacobian analyticalJacobianBody(const TMat& T, const ScrewList& BList, const ThetaList& thetaList){
        Jacobian ret = jacobianBody(BList, thetaList);
        const auto &R = T.topLeftCorner<3, 3>();
        Eigen::Matrix3d LJI = SO3::leftJacobianInverse(SO3::Log(R));
        ret.topRows(3) = LJI * ret.topRows(3);
        ret.bottomRows(3) = R * ret.bottomRows(3);
        return ret;
    }

    /**
     * Compute analytical jacobian d([r, x])/d(q)
     * @param T The pose of the jacobian defined link w.r.t. world
     * @param BList The joint screw axes in the space frame when the
                  manipulator is at the home position, in the format of a
                  matrix with axes as the columns
     * @param thetaList  A list of joint coordinates
     * @return analytical jacobian d([r, x])/d(q)
     */
    inline static Jacobian analyticalJacobianSpace(const TMat& T, const ScrewList& SList, const ThetaList& thetaList){
        Jacobian ret = jacobianSpace(SList, thetaList);
        const auto &R = T.topLeftCorner<3, 3>();
        const auto &t = T.topRightCorner<3, 1>();
        Eigen::Matrix3d LJI = SO3::leftJacobianInverse(SO3::Log(R)) * R.transpose();
        ret.topRows(3) = LJI * ret.topRows(3);
        ret.bottomRows(3) = -SO3::skew(t) * ret.topRows(3) + ret.bottomRows(3);
        return ret;
    }

    /**
     * Compute numerical solution of inverse kinematics
     * @param T The query end-effector pose in SE3;
     * @param angles The joint angles
     * @param M The home configuration
     * @param SList The screw list in the body frame
     * @param eomg error thresh for omega
     * @param ev error thresh for v
     * @return The theta list whose fk result equals to T
     */
    inline static bool numericalIKInBody(const TMat& T, ThetaList& angles, const TMat& M, const ScrewList& BList, double eomg, double ev){
        int i = 0;
        int max_iterations = 50;
        TMat Tfk = fkInBody(M, BList, angles);
        TMat Tdiff = SE3::inv(Tfk) * T;
        SE3::TVec Vb = SE3::Log(Tdiff);
        bool err = (Vb.head(3).norm() > eomg || Vb.tail(3).norm() > ev);
        Jacobian Jb;
        while (err && i++ < max_iterations) {
            Jb = jacobianBody(BList, angles);
            angles += Jb.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Vb);
            // iterate
            Tfk = fkInBody(M, BList, angles);
            Tdiff = SE3::inv(Tfk) * T;
            Vb = SE3::Log(Tdiff);
            err = (Vb.head(3).norm() > eomg || Vb.tail(3).norm() > ev);
        }
        for (int i =0; i < angles.size(); ++i) {
            auto& angle = angles[i];
            angle = angle < SE3::eps ? 0. : angle;
        }
        return !err;
    }

    /**
     * Compute numerical solution of inverse kinematics
     * @param T The query end-effector pose in SE3;
     * @param angles The joint angles
     * @param M The home configuration
     * @param SList The screw list in the space space
     * @param eomg error thresh for omega
     * @param ev error thresh for v
     * @return The theta list whose fk result equals to T
     */
    inline static bool numericalIKInSpace(const TMat& T, ThetaList& angles, const TMat& M, const ScrewList& SList, double eomg, double ev){
        int i = 0;
        int max_iterations = 50;
        TMat Tfk = fkInSpace(M, SList, angles);
        TMat Tdiff = SE3::inv(Tfk) * T;
        SE3::TVec Vs = SE3::adjoint(Tfk, false) * SE3::Log(Tdiff);
        bool err = (Vs.head(3).norm() > eomg || Vs.tail(3).norm() > ev);
        Jacobian Js;
        while (err && i++ < max_iterations) {
            Js = jacobianSpace(SList, angles);
            angles += Js.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(Vs);
            // iterate
            Tfk = fkInSpace(M, SList, angles);
            Tdiff = SE3::inv(Tfk) * T;
            Vs = SE3::adjoint(Tfk, false) * SE3::Log(Tdiff);
            err = (Vs.head(3).norm() > eomg || Vs.tail(3).norm() > ev);
        }
        for (int i =0; i < angles.size(); ++i) {
            auto& angle = angles[i];
            angle = angle < SE3::eps ? 0. : angle;
        }
        return !err;
    }
}




#endif //GCOPTER_ROBOT_HPP
