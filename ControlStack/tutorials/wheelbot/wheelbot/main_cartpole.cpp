

#include<stdbool.h> //for bool
//#include<unistd.h> //for usleep
#include<math.h>

#include "mujoco.h"
#include "GLFW/glfw3.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <iostream>
#include "drake/systems/controllers/linear_quadratic_regulator.h"
#include "../..//utility/numdiff.hpp"



char filename[] 
  = "/home/pang/repo/robotics/ARCHER_hopper/ControlStack/tutorials/wheelbot/wheelbot/cartpole.xml";

// MuJoCo data structures
mjModel* m = NULL;                  // MuJoCo model
mjData* d = NULL;                   // MuJoCo data
mjvCamera cam;                      // abstract camera
mjvOption opt;                      // visualization options
mjvScene scn;                       // abstract scene
mjrContext con;                     // custom GPU context

// mouse interaction
bool button_left = false;
bool button_middle = false;
bool button_right =  false;
double lastx = 0;
double lasty = 0;

// holders of one step history of time and position to calculate dertivatives
mjtNum position_history = 0;
mjtNum previous_time = 0;

// controller related variables
float_t ctrl_update_freq = 100;
mjtNum last_update = 0.0;
mjtNum ctrl;


class ForwardModel {
  public:
  ForwardModel()
  {

  }

  bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) {

    // chasis pitch, chasis pitch vel, wheel_q, wheel_vel          
    d->qpos[0] = parameters[0][0];
    d->qvel[0] = parameters[0][1];
    d->qpos[1] = parameters[0][2];
    d->qvel[1] = parameters[0][3];
    d->ctrl[0] = parameters[1][0];
    mj_forward(m,d);

     double q1dot, q2dot;
    q1dot = d->qvel[0];
    q2dot = d->qvel[1];

  //Equation of motion
  //M*qacc + qfrc_bias = ctrl
  //qacc = inv(M)*(ctrl-qfrc_bias) = q1ddot, q2ddot
  int i;
  const int ndof = 2;
  double M[ndof*ndof] = {0};
  mj_fullM(m,M,d->qM);
  //M = {M[0] M[1]; M[2] M[3]}
  double det_M = M[0]*M[3] - M[1]*M[2];
  double Minv[] = {M[3],-M[1],-M[2],M[0]};
  for (i=0;i<4;i++)
    Minv[i] = Minv[i]/det_M;

  //printf("%f %f %f %f \n",M[0],M[1],M[2],M[3]);

  double qacc[ndof]={0};
  double f[ndof]={0};
  //f = (ctrl-qfrc_bias)
  f[0] = d->ctrl[0]-d->qfrc_bias[0]; 
  f[1] = 0-d->qfrc_bias[1];
  //printf("%f %f \n",f[0],f[1]);

  //qacc = inv(M)*(ctrl-qfrc_bias)
  mju_mulMatVec(qacc,Minv,f,2,2);

  double q1ddot = d->qacc[0];
  double q2ddot = d->qacc[1];

    residuals[0] = q1dot;
    residuals[1] = q1ddot;
    residuals[2] = q2dot;
    residuals[3] = q2ddot;

    return true;

  }

  private:
  

};

void getAB(const mjModel* m, mjData* d, 
  Eigen::Matrix<double, 4, 4, Eigen::RowMajor>& A, 
  Eigen::Matrix<double, 4, 1>& B) {

  double x[4] = {0};
  double u[1] = {0};
  double *parameters[2] = {x, u};
  ForwardModel model;
  
  NumDiff<ForwardModel, 2> numdiff(&model);
  numdiff.df_r_xi<4,4>(parameters, 0, A.data());
  std::cout << "A_: \n" << A << std::endl;
  numdiff.df_r_xi<4,1>(parameters, 1 , B.data());
  std::cout << "B_: \n" << B << std::endl;
}

// keyboard callback
void keyboard(GLFWwindow* window, int key, int scancode, int act, int mods)
{
    // backspace: reset simulation
    if( act==GLFW_PRESS && key==GLFW_KEY_BACKSPACE )
    {
        mj_resetData(m, d);
        mj_forward(m, d);
    }
}

// mouse button callback
void mouse_button(GLFWwindow* window, int button, int act, int mods)
{
    // update button state
    button_left =   (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)==GLFW_PRESS);
    button_middle = (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE)==GLFW_PRESS);
    button_right =  (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT)==GLFW_PRESS);

    // update mouse position
    glfwGetCursorPos(window, &lastx, &lasty);
}


// mouse move callback
void mouse_move(GLFWwindow* window, double xpos, double ypos)
{
    // no buttons down: nothing to do
    if( !button_left && !button_middle && !button_right )
        return;

    // compute mouse displacement, save
    double dx = xpos - lastx;
    double dy = ypos - lasty;
    lastx = xpos;
    lasty = ypos;

    // get current window size
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // get shift key state
    bool mod_shift = (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT)==GLFW_PRESS ||
                      glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT)==GLFW_PRESS);

    // determine action based on mouse button
    mjtMouse action;
    if( button_right )
        action = mod_shift ? mjMOUSE_MOVE_H : mjMOUSE_MOVE_V;
    else if( button_left )
        action = mod_shift ? mjMOUSE_ROTATE_H : mjMOUSE_ROTATE_V;
    else
        action = mjMOUSE_ZOOM;

    // move camera
    mjv_moveCamera(m, action, dx/height, dy/height, &scn, &cam);
}


// scroll callback
void scroll(GLFWwindow* window, double xoffset, double yoffset)
{
    // emulate vertical mouse motion = 5% of window height
    mjv_moveCamera(m, mjMOUSE_ZOOM, 0, -0.05*yoffset, &scn, &cam);
}

Eigen::Matrix<double, 4, 4, Eigen::RowMajor> A;
Eigen::Matrix<double, 4, 1> B;
void init_controller(const mjModel* m, mjData* d)
{
   getAB(m, d, A, B);
}

void mycontroller(const mjModel* m, mjData* d)
{

  // get chasis pitch and rate-of-pitch
  int pendulum_id = 1;


  // add noise 
  double noise;
  mju_standardNormal(&noise);
  d->xfrc_applied[6*pendulum_id+0] = 5*noise;

    bool use_LQR = true;

    if (!use_LQR) {
      

    } else {
        Eigen::Matrix4d Q;
        Eigen::Matrix<double, 1,1 > R;
        Q.setZero();
        Q.diagonal() << 1, 1, 1, 1 ;
        // std::cout << "Q: " << Q << std::endl;
        R << 0.1;
        Eigen::Matrix<mjtNum , 4, 1> N = Eigen::Matrix<mjtNum , 4, 1>::Zero();

        drake::systems::controllers::LinearQuadraticRegulatorResult lqr_result =
            drake::systems::controllers::LinearQuadraticRegulator(A, B, Q, R, N);
        // Eigen::Matrix<mjtNum, 2, 2> S1_ = lqr_result.S;
        Eigen::Matrix<mjtNum, 1, 4> K= lqr_result.K;

        // std::cout << "K: " << K.transpose() << std::endl;
        d->ctrl[0] = -K[0]*d->qpos[0]-K[1]*d->qvel[0]-K[2]*d->qpos[1]-K[3]*d->qvel[1];

    }

}

// main function
int main(int argc, const char** argv)
{

    // load and compile model
    char error[1000] = "Could not load binary model";

    // check command-line arguments
    if( argc<2 )
        m = mj_loadXML(filename, 0, error, 1000);

    else
        if( strlen(argv[1])>4 && !strcmp(argv[1]+strlen(argv[1])-4, ".mjb") )
            m = mj_loadModel(argv[1], 0);
        else
            m = mj_loadXML(argv[1], 0, error, 1000);
    if( !m )
        mju_error_s("Load model error: %s", error);

    // make data
    d = mj_makeData(m);


    // init GLFW
    if( !glfwInit() )
        mju_error("Could not initialize GLFW");

    // create window, make OpenGL context current, request v-sync
    GLFWwindow* window = glfwCreateWindow(1244, 700, "Demo", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // initialize visualization data structures
    mjv_defaultCamera(&cam);
    mjv_defaultOption(&opt);
    mjv_defaultScene(&scn);
    mjr_defaultContext(&con);
    mjv_makeScene(m, &scn, 2000);                // space for 2000 objects
    mjr_makeContext(m, &con, mjFONTSCALE_150);   // model-specific context

    // install GLFW mouse and keyboard callbacks
    glfwSetKeyCallback(window, keyboard);
    glfwSetCursorPosCallback(window, mouse_move);
    glfwSetMouseButtonCallback(window, mouse_button);
    glfwSetScrollCallback(window, scroll);

    double arr_view[] = {90,-45, 4, 0.000000, 0.000000, 0.000000};
    cam.azimuth = arr_view[0];
    cam.elevation = arr_view[1];
    cam.distance = arr_view[2];
    cam.lookat[0] = arr_view[3];
    cam.lookat[1] = arr_view[4];
    cam.lookat[2] = arr_view[5];



    init_controller(m, d);

        // install control callback
    mjcb_control = mycontroller;

    //m->opt.gravity[2]=-1;
    //qpos is dim nqx1 = 7x1; 3 translations + 4 quaternions
    d->qpos[1]=0.010;
    // d->qvel[2]=5;
    // d->qvel[0]=2;
    // use the first while condition if you want to simulate for a period.
    while( !glfwWindowShouldClose(window))
    {
        // advance interactive simulation for 1/60 sec
        //  Assuming MuJoCo can simulate faster than real-time, which it usually can,
        //  this loop will finish on time for the next frame to be rendered at 60 fps.
        //  Otherwise add a cpu timer and exit this loop when it is time to render.
        mjtNum simstart = d->time;
        while( d->time - simstart < 1.0/60.0 )
        {
            mj_step(m, d);

            //drag force = -c*v^2*unit_vector(v); v = sqrt(vx^2+vy^2+vz^2)
            // vector (v) = vx i + vy j + vz k
            //unit_vector(v) = vector(v)/v
            //fx = -c*v*vx;
            //fy = -c*v*vy;
            //fz = -c*v*vz;
            // double vx, vy, vz;
            // vx = d->qvel[0]; vy = d->qvel[1]; vz = d->qvel[2];
            // double v;
            // v = sqrt(vx*vx+vy*vy+vz*vz);
            // double fx, fy, fz;
            // double c = 1;
            // fx = -c*v*vx;
            // fy = -c*v*vy;
            // fz = -c*v*vz;
            // d->qfrc_applied[0]=fx;
            // d->qfrc_applied[1]=fy;
            // d->qfrc_applied[2]=fz;
        }

       // get framebuffer viewport
        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

          // update scene and render
        opt.frame = mjFRAME_WORLD;
        // cam.lookat[0] = d->qpos[0];
        mjv_updateScene(m, d, &opt, NULL, &cam, mjCAT_ALL, &scn);
        mjr_render(viewport, &scn, &con);
        //printf("{%f, %f, %f, %f, %f, %f};\n",cam.azimuth,cam.elevation, cam.distance,cam.lookat[0],cam.lookat[1],cam.lookat[2]);

        // swap OpenGL buffers (blocking call due to v-sync)
        glfwSwapBuffers(window);

        // process pending GUI events, call GLFW callbacks
        glfwPollEvents();

    }

    // free visualization storage
    mjv_freeScene(&scn);
    mjr_freeContext(&con);

    // free MuJoCo model and data, deactivate
    mj_deleteData(d);
    mj_deleteModel(m);

    // terminate GLFW (crashes with Linux NVidia drivers)
    #if defined(__APPLE__) || defined(_WIN32)
        glfwTerminate();
    #endif

    return 1;
}
