

#include<stdbool.h> //for bool
//#include<unistd.h> //for usleep
//#include <math.h>

#include "mujoco.h"
#include "GLFW/glfw3.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <Eigen/Core>
#include <iostream>
#include "drake/systems/controllers/linear_quadratic_regulator.h"
#include "../..//utility/numdiff.hpp"

//simul

//simulation end time
double simend = 10;

#define ndof 2

int lqr=0;

//related to writing data to a file
// FILE *fid;
int loop_index = 0;
const int data_frequency = 10; //frequency at which data is written to a file


std::string path = "/home/pang/repo/robotics/ARCHER_hopper/ControlStack/tutorials/dbpendulum_lqr/dbpendulum_lqr/";
std::string xmlfile = "doublependulum.xml";


// char datafile[] = "data.csv";


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
                           
    d->qpos[0] = parameters[0][0];
    d->qvel[0] = parameters[0][1];
    d->qpos[1] = parameters[0][2];
    d->qvel[1] = parameters[0][3];
    d->ctrl[0] = parameters[1][0];
    mj_forward(m,d);

    double q1dot, q2dot;
    q1dot = d->qvel[0];
    q2dot = d->qvel[1];

    // //Equation of motion
    // //M*qacc + qfrc_bias = ctrl
    // //qacc = inv(M)*(ctrl-qfrc_bias) = q1ddot, q2ddot
    // int i;
    // //const int ndof = 2;
    // double M[ndof*ndof] = {0};
    // mj_fullM(m,M,d->qM);
    // //M = {M[0] M[1]; M[2] M[3]}
    // double det_M = M[0]*M[3] - M[1]*M[2];
    // double Minv[] = {M[3],-M[1],-M[2],M[0]};
    // for (i=0;i<4;i++)
    //   Minv[i] = Minv[i]/det_M;

    // //printf("%f %f %f %f \n",M[0],M[1],M[2],M[3]);

    // double qacc[ndof]={0};
    // double f[ndof]={0};
    // //f = (ctrl-qfrc_bias)
    // f[0] = 0-d->qfrc_bias[0]; //no ctrl because there is no control on link 1
    // f[1] = d->ctrl[0]-d->qfrc_bias[1];
    // //printf("%f %f \n",f[0],f[1]);

    // //qacc = inv(M)*(ctrl-qfrc_bias)
    // mju_mulMatVec(qacc,Minv,f,2,2);

    // double q1ddot = qacc[0];
    // double q2ddot = qacc[1];

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


//****************************
//This function is called once and is used to get the headers
void init_save_data()
{
  //write name of the variable here (header)
  //  fprintf(fid,"t, ");

   //Don't remove the newline
  //  fprintf(fid,"\n");
}

//***************************
//This function is called at a set frequency, put data here
void save_data(const mjModel* m, mjData* d)
{
  //data here should correspond to headers in init_save_data()
  //seperate data by a space %f followed by space
  // fprintf(fid,"%f ",d->time);

  //Don't remove the newline
  // fprintf(fid,"\n");
}

/******************************/
void set_torque_control(const mjModel* m,int actuator_no,int flag)
{
  if (flag==0)
    m->actuator_gainprm[10*actuator_no+0]=0;
  else
    m->actuator_gainprm[10*actuator_no+0]=1;
}
/******************************/


/******************************/
void set_position_servo(const mjModel* m,int actuator_no,double kp)
{
  m->actuator_gainprm[10*actuator_no+0]=kp;
  m->actuator_biasprm[10*actuator_no+1]=-kp;
}
/******************************/

/******************************/
void set_velocity_servo(const mjModel* m,int actuator_no,double kv)
{
  m->actuator_gainprm[10*actuator_no+0]=kv;
  m->actuator_biasprm[10*actuator_no+2]=-kv;
}
/******************************/

void f(const mjModel* m, mjData* d,double input[5], double output[4])
{
  //state = q1, q1dot, q2, q2dot
  //inputs = q1, q1dot, q2, q2dot, u
  //outputs = q1dot, q1ddot, q2dot, q2ddot

  d->qpos[0] = input[0];
  d->qvel[0] = input[1];
  d->qpos[1] = input[2];
  d->qvel[1] = input[3];
  d->ctrl[0] = input[4];
  mj_forward(m,d);

  double q1dot, q2dot;
  q1dot = d->qvel[0];
  q2dot = d->qvel[1];

  //Equation of motion
  //M*qacc + qfrc_bias = ctrl
  //qacc = inv(M)*(ctrl-qfrc_bias) = q1ddot, q2ddot
  int i;
  //const int ndof = 2;
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
  f[0] = 0-d->qfrc_bias[0]; //no ctrl because there is no control on link 1
  f[1] = d->ctrl[0]-d->qfrc_bias[1];
  //printf("%f %f \n",f[0],f[1]);

  //qacc = inv(M)*(ctrl-qfrc_bias)
  mju_mulMatVec(qacc,Minv,f,2,2);

  double q1ddot = qacc[0];
  double q2ddot = qacc[1];

  output[0] = q1dot;
  output[1] = q1ddot;
  output[2] = q2dot;
  output[3] = q2ddot;

}


Eigen::Matrix<double, 4, 4, Eigen::RowMajor> A;
Eigen::Matrix<double, 4, 1> B;
//**************************
void init_controller(const mjModel* m, mjData* d)
{
  // int i,j;
  // double input[5]={0};
  // double output[4]={0};
  // double pert = 0.001;

  // //input[0] = 0.1;
  // //input[1] = 0.1;
  // //input[4] = 0.1;
  // double f0[4] = {0};
  // f(m,d,input,output);
  // //printf("%f %f %f %f \n",output[0],output[1],output[2],output[3]);
  // for (i=0;i<4;i++)
  //   f0[i] = output[i];

  // double A[4][4]={0};

  // j = 0;
  // for (i=0;i<5;i++)
  //   input[i]=0;
  // input[j] = pert;
  // f(m,d,input,output);
  // for (i=0;i<4;i++)
  //   A[i][j] = (output[i]-f0[i])/pert;

  // j = 1;
  // for (i=0;i<5;i++)
  //   input[i]=0;
  // input[j] = pert;
  // f(m,d,input,output);
  // for (i=0;i<4;i++)
  //   A[i][j] = (output[i]-f0[i])/pert;

  // j = 2;
  // for (i=0;i<5;i++)
  //   input[i]=0;
  // input[j] = pert;
  // f(m,d,input,output);
  // for (i=0;i<4;i++)
  //   A[i][j] = (output[i]-f0[i])/pert;

  // j = 3;
  // for (i=0;i<5;i++)
  //   input[i]=0;
  // input[j] = pert;
  // f(m,d,input,output);
  // for (i=0;i<4;i++)
  //   A[i][j] = (output[i]-f0[i])/pert;

  // printf("A = [...\n");
  // for (i=0;i<4;i++)
  // {
  //   for (j=0;j<4;j++)
  //   {
  //     printf("%f ",A[i][j]);
  //   }
  //   printf(";\n");
  // }
  // printf(" ];\n\n");

  // double B[4] = {0};
  // j = 4;
  // for (i=0;i<5;i++)
  //   input[i]=0;
  // input[j] = pert;
  // f(m,d,input,output);
  // for (i=0;i<4;i++)
  //   B[i] = (output[i]-f0[i])/pert;

  //   printf("B = [...\n");
  // for (i=0;i<4;i++)
  //   printf("%f ;\n",B[i]);
  // printf(" ];\n\n");


    getAB(m, d, A, B);


}

//**************************
void mycontroller(const mjModel* m, mjData* d)
{
  //write control here


  if (lqr==1)
  {

    //     // System dynamics matrices: https://github.com/RobotLocomotion/drake/blob/master/systems/controllers/zmp_planner.h#L306
    // Eigen::Matrix<mjtNum , 2, 2> A_;
    // Eigen::Matrix<mjtNum , 2, 1> B_;
    // Eigen::Matrix<mjtNum , 2, 2> Q_ = Eigen::Matrix<mjtNum , 2, 2>::Identity();
    // Eigen::Matrix<mjtNum , 1, 1> R_ = Eigen::Matrix<mjtNum , 1, 1>::Identity();
    // Eigen::Matrix<mjtNum , 2, 1> N = Eigen::Matrix<mjtNum , 2, 1>::Zero();


    // // Eqs comming from https://youtu.be/XA6B1IxALhk?t=30m3s
    // A_(0,0) = 0.0;
    // A_(0,1) = 1.0;
    // A_(1,0) = - m->opt.gravity[2]*mju_cos(d->sensordata[0])/(m->geom_size[rod_id*3+1]*2);
    // A_(1,1) = - m->dof_damping[pivot_id]/(m->body_mass[ball_id]*mju_pow(m->geom_size[rod_id*3+1]*2,2.0));
    // B_(0,0) = 0.0;
    // B_(1,0) = 1.0;
    // Q_ *= 10.0;

//     A = [...
// 0.000000 1.000000 0.000000 0.000000 ;
// 12.124717 0.000000 -11.022452 -0.000000 ;
// 0.000000 0.000000 0.000000 1.000000 ;
// -15.431458 -0.000000 40.783088 0.000000 ;
//  ];

// B = [...
// 0.000000 ;
// -3.820225 ;
// 0.000000 ;
// 12.134831 ;
//  ];

    // Eigen::Matrix4d A;
    // A <<  0.000000, 1.000000, 0.000000, 0.000000,
    // 12.124717, 0.000000, -11.022452, -0.000000,
    // 0.000000, 0.000000, 0.000000, 1.000000,
    // -15.431458, -0.000000, 40.783088, 0.000000;
    // Eigen::Vector4d B(
    //   0.000000,
    // -3.820225,
    // 0.000000,
    // 12.134831
    // );


    Eigen::Matrix4d Q;
    Eigen::Matrix<double, 1,1 > R;
    Q.setIdentity();
    R << 0.1;
    Eigen::Matrix<mjtNum , 4, 1> N = Eigen::Matrix<mjtNum , 4, 1>::Zero();


    drake::systems::controllers::LinearQuadraticRegulatorResult lqr_result =
            drake::systems::controllers::LinearQuadraticRegulator(A, B, Q, R, N);
    // Eigen::Matrix<mjtNum, 2, 2> S1_ = lqr_result.S;
    Eigen::Matrix<mjtNum, 1, 4> K= lqr_result.K;

    std::cout << "K: " << K.transpose() << std::endl;

    //double K[4]={-265.4197,  -97.9928 , -66.4967,  -28.8720};
    // double K[4] = { -1.2342*1000,   -0.4575*1000,   -0.3158 *1000,  -0.1330*1000};
    d->ctrl[0] = -K[0]*d->qpos[0]-K[1]*d->qvel[0]-K[2]*d->qpos[1]-K[3]*d->qvel[1];
  }

  if (lqr==1)
  {
    double noise;
    mju_standardNormal(&noise);
    int body = 2;
    d->xfrc_applied[6*body+0] = 2*noise;
    body = 1;
    d->xfrc_applied[6*body+0] = 2*noise;
    d->qfrc_applied[0]=noise;
    d->qfrc_applied[1]=noise;
  }

  //write data here (dont change/dete this function call; instead write what you need to save in save_data)
  if ( loop_index%data_frequency==0)
    {
      save_data(m,d);
    }
  loop_index = loop_index + 1;
}


//************************
// main function
int main(int argc, const char** argv)
{
    // char xmlpath[100]={};
    // char datapath[100]={};

    // strcat(xmlpath,path);
    // strcat(xmlpath,xmlfile);

    // strcat(datapath,path);
    // strcat(datapath,datafile);


    // load and compile model
    char error[1000] = "Could not load binary model";

    // check command-line arguments
    if( argc<2 )
        m = mj_loadXML(std::string(path + xmlfile).c_str(), 0, error, 1000);

    // else
    //     if( strlen(argv[1])>4 && !strcmp(argv[1]+strlen(argv[1])-4, ".mjb") )
    //         m = mj_loadModel(argv[1], 0);
    //     else
    //         m = mj_loadXML(argv[1], 0, error, 1000);
    // if( !m )
    //     mju_error_s("Load model error: %s", error);

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

    double arr_view[] = {89.608063, -11.588379, 5, 0.000000, 0.000000, 2.000000}; //view the left side (for ll, lh, left_side)
    cam.azimuth = arr_view[0];
    cam.elevation = arr_view[1];
    cam.distance = arr_view[2];
    cam.lookat[0] = arr_view[3];
    cam.lookat[1] = arr_view[4];
    cam.lookat[2] = arr_view[5];

    // install control callback
    mjcb_control = mycontroller;

    // fid = fopen(datapath,"w");
    // init_save_data();
    init_controller(m,d);

    //d->qvel[1] = -0.2;
    lqr = 1;

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
        }

        if (d->time>=simend)
        {
          //  fclose(fid);
           break;
         }

       // get framebuffer viewport
        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

          // update scene and render
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
