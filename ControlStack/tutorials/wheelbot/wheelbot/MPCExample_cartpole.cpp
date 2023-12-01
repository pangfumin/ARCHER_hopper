/**
 * @file MPCExample.cpp
 * @author Giulio Romualdi
 * @copyright Released under the terms of the BSD 3-Clause License
 * @date 2018
 */

// osqp-eigen
#include "OsqpEigen/OsqpEigen.h"

// eigen
#include <Eigen/Dense>

#include <iostream>

#define STATE_DIM 4
#define INPUT_DIM 1

    // M = 0.5
    // m = 0.2
    // b = 0.1
    // ftheta = 0.1
    // l = 0.3
    // g = 9.81

    // Ts = 50e-3
const double M = 0.5;
const double m = 0.2;
const double b = 0.1;
const double ftheta = 0.1;
const double l = 0.3;
const double g = 9.81;
const double Ts = 50e-3;


void setDynamicsMatrices(Eigen::Matrix<double, STATE_DIM, STATE_DIM>& Ad, Eigen::Matrix<double, STATE_DIM, INPUT_DIM>& Bd)
{
    // a << 1., 0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.1, 0., 0.,
    //     0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.1, 0., 0., 0., 0.0488, 0., 0., 1., 0., 0., 0.0016,
    //     0., 0., 0.0992, 0., 0., 0., -0.0488, 0., 0., 1., 0., 0., -0.0016, 0., 0., 0.0992, 0., 0.,
    //     0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.0992, 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
    //     0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,
    //     0., 0., 0.9734, 0., 0., 0., 0., 0., 0.0488, 0., 0., 0.9846, 0., 0., 0., -0.9734, 0., 0., 0.,
    //     0., 0., -0.0488, 0., 0., 0.9846, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.9846;

    // b << 0., -0.0726, 0., 0.0726, -0.0726, 0., 0.0726, 0., -0.0152, 0.0152, -0.0152, 0.0152, -0.,
    //     -0.0006, -0., 0.0006, 0.0006, 0., -0.0006, 0.0000, 0.0106, 0.0106, 0.0106, 0.0106, 0,
    //     -1.4512, 0., 1.4512, -1.4512, 0., 1.4512, 0., -0.3049, 0.3049, -0.3049, 0.3049, -0.,
    //     -0.0236, 0., 0.0236, 0.0236, 0., -0.0236, 0., 0.2107, 0.2107, 0.2107, 0.2107;

    //  Ac =np.array([[0,       1,          0,                  0],
    //               [0,       -b/M,       -(g*m)/M,           (ftheta*m)/M],
    //               [0,       0,          0,                  1],
    //               [0,       b/(M*l),    (M*g + g*m)/(M*l),  -(M*ftheta + ftheta*m)/(M*l)]])

    // Bc = np.array([
    //     [0.0],
    //     [1.0/M],
    //     [0.0],
    //     [-1/(M*l)]
    // ])

    // [nx, nu] = Bc.shape # number of states and number or inputs

    // # Simple forward euler discretization
    // Ad = np.eye(nx) + Ac*Ts
    // Bd = Bc*Ts

//     A_: 
//        0        1        0        0
//  25.6246        0        0        0
//        0        0        0        1
// -11.4737        0        0        0
// B_: 
//        0
// -7.56335
//        0
//  9.35673

    Eigen::Matrix<double, 4, 4, Eigen::RowMajor> Ac;
    // Ac <<   0, 1, 0, 0,
    //         0, -b/M, -(g*m)/M, (ftheta*m)/M,
    //         0, 0, 0, 1,
    //         0, b/(M*l), (M*g + g*m)/(M*l), -(M*ftheta + ftheta*m)/(M*l);

      Ac <<         0,       1,        0,        0,
                25.6246 ,       0 ,       0 ,       0,
                    0,        0 ,       0 ,       1,
                -11.4737,        0,        0 ,       0;

    Eigen::Matrix<double, 4, 1> Bc;
    // Bc << 0.0, 1.0/M, 0.0, -1.0/(M * l);
    Bc << 
                0,
        -7.56335,
            0,
        9.35673;


    Ad = Eigen::Matrix4d::Identity() + Ac * Ts;
    Bd = Bc * Ts;



}

void setInequalityConstraints(Eigen::Matrix<double, STATE_DIM, 1>& xMax,
                              Eigen::Matrix<double, STATE_DIM, 1>& xMin,
                              Eigen::Matrix<double, INPUT_DIM, 1>& uMax,
                              Eigen::Matrix<double, INPUT_DIM, 1>& uMin)
{
    //     # Constraints
    // xmin = np.array([-1.0, -100, -100, -100])
    // xmax = np.array([0.3,   100.0, 100, 100])

    // umin = np.array([-20])
    // umax = np.array([20])

    // input inequality constraints
    uMin << -200;

    uMax << 200;

    // state inequality constraints
    xMin << -100.0, -100, -100.0, -100;

    xMax << 100.0, 100, 100.0, 100;
}

void setWeightMatrices(Eigen::DiagonalMatrix<double, STATE_DIM>& Q, Eigen::DiagonalMatrix<double, INPUT_DIM>& R)
{
    Q.diagonal() << 100, 0.1, 0.010, 100;
    R.diagonal() << 1.00;
}

void castMPCToQPHessian(const Eigen::DiagonalMatrix<double, STATE_DIM>& Q,
                        const Eigen::DiagonalMatrix<double, INPUT_DIM>& R,
                        int mpcWindow,
                        Eigen::SparseMatrix<double>& hessianMatrix)
{

    hessianMatrix.resize(STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow,
                         STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow);

    // populate hessian matrix
    for (int i = 0; i < STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow; i++)
    {
        if (i < STATE_DIM * (mpcWindow + 1))
        {
            int posQ = i % STATE_DIM;
            float value = Q.diagonal()[posQ];
            if (value != 0)
                hessianMatrix.insert(i, i) = value;
        } else
        {
            int posR = i % INPUT_DIM;
            float value = R.diagonal()[posR];
            if (value != 0)
                hessianMatrix.insert(i, i) = value;
        }
    }
}

void castMPCToQPGradient(const Eigen::DiagonalMatrix<double, STATE_DIM>& Q,
                         const Eigen::Matrix<double, STATE_DIM, 1>& xRef,
                         int mpcWindow,
                         Eigen::VectorXd& gradient)
{

    Eigen::Matrix<double, STATE_DIM, 1> Qx_ref;
    Qx_ref = Q * (-xRef);

    // populate the gradient vector
    gradient = Eigen::VectorXd::Zero(STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow, 1);
    for (int i = 0; i < STATE_DIM * (mpcWindow + 1); i++)
    {
        int posQ = i % STATE_DIM;
        float value = Qx_ref(posQ, 0);
        gradient(i, 0) = value;
    }
}

void castMPCToQPConstraintMatrix(const Eigen::Matrix<double, STATE_DIM, STATE_DIM>& dynamicMatrix,
                                 const Eigen::Matrix<double, STATE_DIM, INPUT_DIM>& controlMatrix,
                                 int mpcWindow,
                                 Eigen::SparseMatrix<double>& constraintMatrix)
{
    constraintMatrix.resize(STATE_DIM * (mpcWindow + 1) + STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow,
                            STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow);

    // populate linear constraint matrix
    for (int i = 0; i < STATE_DIM * (mpcWindow + 1); i++)
    {
        constraintMatrix.insert(i, i) = -1;
    }

    for (int i = 0; i < mpcWindow; i++)
        for (int j = 0; j < STATE_DIM; j++)
            for (int k = 0; k < STATE_DIM; k++)
            {
                float value = dynamicMatrix(j, k);
                if (value != 0)
                {
                    constraintMatrix.insert(STATE_DIM * (i + 1) + j, STATE_DIM * i + k) = value;
                }
            }

    for (int i = 0; i < mpcWindow; i++)
        for (int j = 0; j < STATE_DIM; j++)
            for (int k = 0; k < INPUT_DIM; k++)
            {
                float value = controlMatrix(j, k);
                if (value != 0)
                {
                    constraintMatrix.insert(STATE_DIM * (i + 1) + j, INPUT_DIM * i + k + STATE_DIM * (mpcWindow + 1))
                        = value;
                }
            }

    for (int i = 0; i < STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow; i++)
    {
        constraintMatrix.insert(i + (mpcWindow + 1) * STATE_DIM, i) = 1;
    }
}

void castMPCToQPConstraintVectors(const Eigen::Matrix<double, STATE_DIM, 1>& xMax,
                                  const Eigen::Matrix<double, STATE_DIM, 1>& xMin,
                                  const Eigen::Matrix<double, INPUT_DIM, 1>& uMax,
                                  const Eigen::Matrix<double, INPUT_DIM, 1>& uMin,
                                  const Eigen::Matrix<double, STATE_DIM, 1>& x0,
                                  int mpcWindow,
                                  Eigen::VectorXd& lowerBound,
                                  Eigen::VectorXd& upperBound)
{
    // evaluate the lower and the upper inequality vectors
    Eigen::VectorXd lowerInequality
        = Eigen::MatrixXd::Zero(STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow, 1);
    Eigen::VectorXd upperInequality
        = Eigen::MatrixXd::Zero(STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow, 1);
    for (int i = 0; i < mpcWindow + 1; i++)
    {
        lowerInequality.block(STATE_DIM * i, 0, STATE_DIM, 1) = xMin;
        upperInequality.block(STATE_DIM * i, 0, STATE_DIM, 1) = xMax;
    }
    for (int i = 0; i < mpcWindow; i++)
    {
        lowerInequality.block(INPUT_DIM * i + STATE_DIM * (mpcWindow + 1), 0, INPUT_DIM, 1) = uMin;
        upperInequality.block(INPUT_DIM * i + STATE_DIM * (mpcWindow + 1), 0, INPUT_DIM, 1) = uMax;
    }

    // evaluate the lower and the upper equality vectors
    Eigen::VectorXd lowerEquality = Eigen::MatrixXd::Zero(STATE_DIM * (mpcWindow + 1), 1);
    Eigen::VectorXd upperEquality;
    lowerEquality.block(0, 0, STATE_DIM, 1) = -x0;
    upperEquality = lowerEquality;
    lowerEquality = lowerEquality;

    // merge inequality and equality vectors
    lowerBound = Eigen::MatrixXd::Zero(2 * STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow, 1);
    lowerBound << lowerEquality, lowerInequality;

    upperBound = Eigen::MatrixXd::Zero(2 * STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow, 1);
    upperBound << upperEquality, upperInequality;
}

void updateConstraintVectors(const Eigen::Matrix<double, STATE_DIM, 1>& x0,
                             Eigen::VectorXd& lowerBound,
                             Eigen::VectorXd& upperBound)
{
    lowerBound.block(0, 0, STATE_DIM, 1) = -x0;
    upperBound.block(0, 0, STATE_DIM, 1) = -x0;
}

double getErrorNorm(const Eigen::Matrix<double, STATE_DIM, 1>& x, const Eigen::Matrix<double, STATE_DIM, 1>& xRef)
{
    // evaluate the error
    Eigen::Matrix<double, STATE_DIM, 1> error = x - xRef;

    // return the norm
    return error.norm();
}

int main()
{
    // set the preview window
    int mpcWindow = 20;

    // allocate the dynamics matrices
    Eigen::Matrix<double, STATE_DIM, STATE_DIM> Ad;
    Eigen::Matrix<double, STATE_DIM, INPUT_DIM> Bd;

    // allocate the constraints vector
    Eigen::Matrix<double, STATE_DIM, 1> xMax;
    Eigen::Matrix<double, STATE_DIM, 1> xMin;
    Eigen::Matrix<double, INPUT_DIM, 1> uMax;
    Eigen::Matrix<double, INPUT_DIM, 1> uMin;

    // allocate the weight matrices
    Eigen::DiagonalMatrix<double, STATE_DIM> Q;
    Eigen::DiagonalMatrix<double, INPUT_DIM> R;

    // allocate the initial and the reference state space
    Eigen::Matrix<double, STATE_DIM, 1> x0;
    Eigen::Matrix<double, STATE_DIM, 1> xRef;

    // set the initial and the desired states
    double phi0 = 20*2*M_PI/360;

    x0 << phi0, 0, 1, 0;// # initial state
    Eigen::Matrix<double, STATE_DIM, 1> x_init = x0;

    xRef << 0.0, 0.0, 0.0, 0.12;

    // set MPC problem quantities
    setDynamicsMatrices(Ad, Bd);
    setInequalityConstraints(xMax, xMin, uMax, uMin);
    setWeightMatrices(Q, R);

    // allocate QP problem matrices and vectores
    Eigen::SparseMatrix<double> hessian;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> linearMatrix;
    Eigen::VectorXd lowerBound;
    Eigen::VectorXd upperBound;

    // cast the MPC problem as QP problem
    castMPCToQPHessian(Q, R, mpcWindow, hessian);
    castMPCToQPGradient(Q, xRef, mpcWindow, gradient);
    castMPCToQPConstraintMatrix(Ad, Bd, mpcWindow, linearMatrix);
    castMPCToQPConstraintVectors(xMax, xMin, uMax, uMin, x0, mpcWindow, lowerBound, upperBound);

    // instantiate the solver
    OsqpEigen::Solver solver;

    // settings
    // solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // set the initial data of the QP solver
    solver.data()->setNumberOfVariables(STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow);
    solver.data()->setNumberOfConstraints(2 * STATE_DIM * (mpcWindow + 1) + INPUT_DIM * mpcWindow);
    if (!solver.data()->setHessianMatrix(hessian))
        return 1;
    if (!solver.data()->setGradient(gradient))
        return 1;
    if (!solver.data()->setLinearConstraintsMatrix(linearMatrix))
        return 1;
    if (!solver.data()->setLowerBound(lowerBound))
        return 1;
    if (!solver.data()->setUpperBound(upperBound))
        return 1;

    // instantiate the solver
    if (!solver.initSolver())
        return 1;

    // len_sim = 20 # simulation length (s)
    // nsim = int(len_sim/Ts) # simulation length(timesteps)
    // xsim = np.zeros((nsim,nx))
    // usim = np.zeros((nsim,nu))
    // tsim = np.arange(0,nsim)*Ts

    double len_sim = 50; // second
    int nsim = (int) (len_sim / Ts);

    // controller input and QPSolution vector
    Eigen::Matrix<double, 1,1> ctr;
    Eigen::VectorXd QPSolution;

    // // number of iteration steps
    // int numberOfSteps = 50;

    for (int i = 0; i < nsim; i++)
    {
        // solve the QP problem
        if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError)
            return 1;

        // get the controller input
        QPSolution = solver.getSolution();
        ctr = QPSolution.block(STATE_DIM * (mpcWindow + 1), 0, INPUT_DIM, 1);

        // save data into file
        auto x0Data = x0.data();

        // propagate the model
        x0 = Ad * x0 + Bd * ctr;

        std::cout << "x_init: " << i << " " << x_init.transpose() << std::endl;
        std::cout << "i:      " << i << " " << x0.transpose()  << " ---- " << ctr  << std::endl;
        



        // update the constraint bound
        updateConstraintVectors(x0, lowerBound, upperBound);
        if (!solver.updateBounds(lowerBound, upperBound))
            return 1;
    }
    return 0;
}
