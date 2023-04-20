/**
 * @brief A wrapper class for a cpp-style class-based interface.
 * Methods to obtain Hessian and Jacoiban of R/T w.r.t. to q can be found in
 * "A Systematic Approach to Computing the Manipulator Jacobian and Hessian using the Elementary Transform Sequence"
 * "Jesse Haviland, Peter Corke"
 * Thanks to their excellent work.
 */
#ifndef MANIPULATORTOOLS_WRAPPER_HPP
#define MANIPULATORTOOLS_WRAPPER_HPP

#include "robot.hpp"

namespace robot {
    typedef std::vector<Eigen::Matrix<double, 6, -1>> Hessian; ///< Every element denotes that \f(dJ / dq_i\f)
    typedef std::vector<Eigen::Matrix<double, 3, 3>> RJacobian; ///< Every element denotes that \f(dR / dq_i\f)
    typedef std::vector<Eigen::Matrix<double, 4, 4>> TJacobian; ///< Every element denotes that \f(dT / dq_i\f)

    /**
     * @brief A wrapper class for a cpp-style interface
     * @tparam m The number of active joints
     */
    template<int m>
    class ManipulatorTool {
        ManipulatorTool() = default;

        ~ManipulatorTool() = default;

    public:
        /**
         * @brief set screw list through a c-style fk function
         * @param fk a c-style fk function
         */
        void setScrew(lfk_t fk) {
            M = robot::computeScrewList(fk, sList, false, m);
            screwReady = true;
        }

        /**
         * @brief set screw list through a cpp-style fk function
         * @param fk a cpp-style fk function
         */
        void setScrew(const lfk_func &fk) {
            M = robot::computeScrewList(fk, sList, false, m);
            screwReady = true;
        }

        /**
         * @brief set screw list through a pre-computed result
         * @param screwList A list of pre-computed screw axes
         * @param T The initial transform M
         */
        void setScrew(const ScrewList &screwList, const TMat &T) {
            if (screwList.cols() != m) {
                throw (std::invalid_argument("The screw list provided only has " + std::to_string(screwList.cols())
                                             + "columns but " + std::to_string(m) + "needed"));
            } else {
                sList = screwList;
                M = SE3::project(T);
            }
        }

        /**
         * @brief A interface of forward kinematics function
         * @param q The current joint values
         * @return The transformation T of end effector w.r.t. the base
         */
        TMat getT(const ThetaList &q) {
            if (screwReady)
                return fkInSpace(sList, q);
            else {
                throw (std::invalid_argument("No screw lists have been set, call setScrew() at first"));
            }
        }

        /**
         * @brief To get jacobian \f(J = dg/dq, g = [t, \omega]^T\f), where g is defined as the general coordinates of
         * end effector w.r.t. the base.
         * @param q The current joint values
         * @return The jacobian \f(6\times n\f)
         */
        Jacobian getJacobian(const ThetaList &q) {
            if (screwReady) {
                Jacobian J = jacobianSpace(sList, q);
                TMat T = fkInSpace(M, sList, q);
                return fromMRSpaceJacobian(T, J, false);
            } else {
                throw (std::invalid_argument("No screw lists have been set, call setScrew() at first"));
            }
        }

        /**
         * @brief To get Hessian \f(H = dJ / dq\f), note that this function computes jacobian, if you have a pre-computed
         * jacobian use another interface.
         * @param q The current joint values
         * @param H The output Hessian
         */
        void getHessian(const ThetaList &q, Hessian &H) {
            Jacobian J = getJacobian(q);
            getHessian(J, H);
        }

        /**
         * @brief To get Hessian \f(H = dJ/ dq\f) thorough pre-computed jacobian.
         * @param J The pre-computed jacobain, must be in the format of \f(J_a, J_w\f). The rotation part is located downside.
         * @param H The output Hessian
         */
        static void getHessian(const Jacobian &J, Hessian &H) {
            const Eigen::Matrix<double, 3, m> &Ja = J.topRows(3);
            const Eigen::Matrix<double, 3, m> &Jw = J.bottomRows(3);

            H.resize(m);
            for (int i = 0; i < m; ++i) {
                auto &Hi = H[i];
                Hi.resize(6, m);
                Hi.setZero();
                for (int j = 0; j < m; ++j) {
                    Hi.col(j).tail(3) = SO3::skew(Jw.col(j)) * Jw.col(i);
                    if (j >= i) {
                        Hi.col(j).head(3) = SO3::skew(Jw.col(i)) * Ja.col(j);
                        if (j != i)
                            H[j].col(i).head(3) = Hi.col(j).head(3);
                    }
                }
            }
        }

        /**
         * @brief To get the jacobian of R w.r.t. q, which should be a tensor but store in a vecotr format
         * @param q The current joint values. Notice that this function computes R and J. Use another static function
         * with pre-computed R and J in case of efficiency
         * @param JR The output jacobian of R.
         */
        void getRJacobian(const ThetaList &q, RJacobian &JR) {
            TMat T = getT(q);
            const Eigen::Matrix<double, 3, 3> &R = T.topLeftCorner<3, 3>();

            Jacobian J = fromMRSpaceJacobian(T, jacobianSpace(sList, q), false);
            getRJacobian(J, R, JR);
        }

        /**
         * @brief To get the jacobian of R w.r.t. q, which should be a tensor but store in a vecotr format
         * @param J The pre-computed jacobian, must in the format of \f(J_a, J_w\f)
         * @param R The pre-computed rotation matrix
         * @param JR The output Jacobian vector.
         */
        static void getRJacobian(const Jacobian &J, const Eigen::Matrix3d &R, RJacobian &JR) {
            const Eigen::Matrix<double, 3, m> &Jw = J.bottomRows(3);

            JR.resize(m);
            for (int i = 0; i < m; ++i) {
                JR[i] = SO3::skew(Jw.col(i)) * R;
            }
        }

        /**
         * @brief To get the jacobain of T w.r.t. q, which should be a tensor but stored in vector;
         * @param q The current joint values, notice that this function computes T and J. Use another static function
         * with pre-computed T and J in case of efficiency
         * @param JT The output jacobain vector
         */
        void getTJacobian(const ThetaList &q, TJacobian &JT) {
            TMat T = fkInSpace(q, sList);
            Jacobian J = fromMRSpaceJacobian(T, jacobianSpace(sList, q), false);

            getTJacobian(T, J, JT);
        }

        /**
         * @brief To get the jacobain of T w.r.t. q, which should be a tensor but stored in vector
         * @param T The pre-computed transformation T
         * @param J The pre-computed jacobian, must in the format of \f(J_a, J_w\f)
         * @param JT The output jacobain
         */
        static void getTJacobian(const TMat &T, const Jacobian &J, TJacobian &JT) {
            const Eigen::Matrix<double, 3, m> &Jw = J.bottomRows(3);
            const Eigen::Matrix<double, 3, m> &Ja = J.topRows(3);
            const Eigen::Matrix<double, 3, 3> &R = T.topLeftCorner<3, 3>();

            JT.resize(m);
            for (int i = 0; i < m; ++i) {
                auto &JTi = JT[i];
                JTi.setZero();
                JTi.topLeftCorner<3, 3>() = SO3::skew(Jw.col(i)) * R;
                JTi.topRightCorner<3, 1>() = Ja.col(i);
            }
        }

        /**
         * @brief To compute end-effector general velocity \f(v = [\dot{t}, \theta \hat{\omega}]\f)
         * @param q The current joint values, notice that this function computes jacoiban, use another static function
         * with pre-computed Jacobian in case of efficiency
         * @param qv The current joint velocity
         * @return The general end effector velocity
         */
        Eigen::Matrix<double, 6, 1> getVelocity(const ThetaList &q, const ThetaList &qv) {
            assert(qv.size() == q.size() && qv.size() == m);
            return getJacobian(q) * qv;
        }

        /**
         * @brief To compute end-effector general acceleration \f(a = [\ddot{a}, \alpha]\f)
         * @param q Th current joint values, notice that this function computes jaociban and hessian, use another static
         * fucntion with pre-computed jacobian and hessian in case of efficiency
         * @param qv The current joint velocity
         * @param qa The current joint acceleration
         * @return The general end effector acceleration
         */
        Eigen::Matrix<double, 6, 1> getAcceleration(const ThetaList &q, const ThetaList &qv, const ThetaList &qa) {
            assert(q.size() == m);
            assert(qv.size() == qa.size() && qv.size() == m);
            Jacobian J = getJacobian(q);

            Hessian H;
            getHessian(J, H);

            getAcceleration(J, H, qv, qa);
        }

        /**
         * @brief To compute end-effector general acceleration \f(a = [\ddot{a}, \alpha]\f)
         * @param J The pre-computed jacobian
         * @param H The pre-computed hessian
         * @param qv The current joint velocity
         * @param qa The current joint acceleration
         * @return The general end effector acceleration
         */
        static Eigen::Matrix<double, 6, 1>
        getAcceleration(const Jacobian &J, const Hessian &H, const ThetaList &qv, const ThetaList &qa) {
            assert(H.size() == m && qv.size() == qa.size() && qv.size() == m);
            assert(J.cols() == m && H[0].cols() == m);
            Eigen::Matrix<double, 6, m> Hq;
            Hq.setZero();
            for (int i = 0; i < m; ++i) {
                Hq += H[i] * qv[i];
            }
            return Hq * qv + J * qa;
        }


    protected:
        bool screwReady{false};
        Eigen::Matrix<double, 6, m> sList;
        TMat M;
    };
}


#endif //MANIPULATORTOOLS_WRAPPER_HPP
