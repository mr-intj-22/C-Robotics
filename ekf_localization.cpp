#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
// #include <eigen3/unsupported/Eigen>

using namespace Eigen;
using namespace std;

const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ");
const float MAX_TIME = 63, DT = 0.01;

Matrix<double, 4, 4> Q;
Matrix<double, 2, 2> R;

Matrix<double, 2, 2> GPS_NOISE;
Matrix<double, 2, 2> INPUT_NOISE;

struct Vars
{
    Matrix<double, 4, 1> x_true;
    Matrix<double, 2, 1> z;
    Matrix<double, 4, 1> x_dr;
    Matrix<double, 2, 1> ud;

    Matrix<double, 4, 1> x_Est;
    Matrix<double, 4, 4> P_Est;
};

Matrix<double, 2, 1> calc_input()
{
    float v = 1.0;       // [m/s]
    float yawrate = 0.1; // [rad/s]
    Matrix<double, 2, 1> u;
    u << v, yawrate;
    return u;
};

Matrix<double, 4, 1> motion_model(Matrix<double, 4, 1> x, Matrix<double, 2, 1> u)
{
    Matrix<double, 4, 4> F;
    F << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 0;

    Matrix<double, 4, 2> B;
    B << DT * cos(x[2]), 0,
        DT * sin(x[2]), 0,
        0, DT,
        1, 0;

    return F * x + B * u;
};

Matrix<double, 2, 1> observ_model(Matrix<double, 4, 1> x)
{
    Matrix<double, 2, 4> H;
    H << 1, 0, 0, 0,
        0, 1, 0, 0;

    return H * x;
};

void observation(Vars *vars, Matrix<double, 2, 1> u)
{

    vars->x_true << motion_model(vars->x_true, u);

    // add noise to gps x-y
    vars->z << observ_model(vars->x_true) + GPS_NOISE * MatrixXd::Random(2, 1);

    // add noise to input
    vars->ud << u + INPUT_NOISE * MatrixXd::Random(2, 1);

    vars->x_dr << motion_model(vars->x_dr, vars->ud);
};

Matrix<double, 4, 4> jacobF(Vars *vars)
{

    double v = calc_input()[0];
    double yaw = vars->x_Est[2];

    Matrix<double, 4, 4> JF;
    JF << 1, 0, -DT * v * sin(yaw), DT * cos(yaw),
        0, 1, DT * v * cos(yaw), DT * sin(yaw),
        0, 0, 1, 0,
        0, 0, 0, 1;

    return JF;
};

Matrix<double, 2, 4> jacobH()
{

    Matrix<double, 2, 4> JH;
    JH << 1, 0, 0, 0,
        0, 1, 0, 0;

    return JH;
};

void ekf_estimation(Vars *vars)
{
    //  Predict
    Matrix<double, 4, 1> xPred = motion_model(vars->x_Est, vars->ud);
    Matrix<double, 4, 4> jF = jacobF(vars);
    Matrix<double, 4, 4> PPred = jF * vars->P_Est * jF.transpose() + Q;

    // //  Update
    Matrix<double, 2, 4> jH = jacobH();
    Matrix<double, 2, 1> zPred = observ_model(xPred);
    Matrix<double, 2, 1> y = vars->z - zPred;
    Matrix<double, 2, 2> S = jH * PPred * jH.transpose() + R;
    Matrix<double, 4, 2> K = PPred * jH.transpose() * S.inverse();
    vars->x_Est = xPred + K * y;
    vars->P_Est = (Matrix<double, 4, 4>::Identity() - K * jH) * PPred;
}

int main()
{

    // process noise (model imperfections, non-linearities, other contributions)
    Q << 0.1, 0, 0, 0,
        0, 0.1, 0, 0,
        0, 0, 1.0, 0,
        0, 0, 0, 1.0;
    // measurement noise (sensor, environment related)
    R << pow(0.1, 2), 0,
        0, pow(0.1, 2);
    // GPS Noise (for simulation reasons)
    GPS_NOISE << pow(0.5, 2), 0,
        0, pow(0.5, 2);
    // Input Noise (for simulation reasons)
    INPUT_NOISE << pow(1, 2), 0,
        0, pow(30 * M_PI/180.0, 2);

    // initialise simulation data
    Vars vars = Vars();
    vars.x_Est.setZero();
    vars.x_true.setZero();
    vars.x_dr.setZero();
    vars.P_Est << Matrix<double, 4, 4>::Identity();

    float time = 0.0;

    ofstream file;
    file.open("EKF.csv");
    char comma = ',';
    file << time << comma << vars.x_true.transpose().format(CSVFormat) << comma << vars.x_dr.transpose().format(CSVFormat) << comma << vars.x_Est.transpose().format(CSVFormat) << comma << vars.z.transpose().format(CSVFormat) << endl;

    while (time < MAX_TIME)
    {
        time += DT;
        Matrix<double, 2, 1> u = calc_input();
        observation(&vars, u);
        ekf_estimation(&vars);

        file << time << comma << vars.x_true.transpose().format(CSVFormat) << comma << vars.x_dr.transpose().format(CSVFormat) << comma << vars.x_Est.transpose().format(CSVFormat) << comma << vars.z.transpose().format(CSVFormat) << endl;
    };
    file.close();

    return 0;
};