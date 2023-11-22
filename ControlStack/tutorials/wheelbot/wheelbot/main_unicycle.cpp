

#include<stdbool.h> //for bool
//#include<unistd.h> //for usleep
#include<math.h>

#include "mujoco.h"
#include "GLFW/glfw3.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Geometry>


//simulation end time
double simend = 10;

char filename[] 
  = "/home/pang/repo/robotics/ARCHER_hopper/ControlStack/tutorials/wheelbot/wheelbot/unicycle.xml";

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

template<typename T>
Eigen::Matrix<T, 3, 3> RPY2mat(T roll, T pitch, T yaw)
{
    Eigen::Matrix<T, 3, 3> m;

    T cr = cos(roll);
    T sr = sin(roll);
    T cp = cos(pitch);
    T sp = sin(pitch);
    T cy = cos(yaw);
    T sy = sin(yaw);

    m(0,0) = cy * cp;
    m(0,1) = cy * sp * sr - sy * cr;
    m(0,2) = cy * sp * cr + sy * sr;
    m(1,0) = sy * cp;
    m(1,1) = sy * sp * sr + cy * cr;
    m(1,2) = sy * sp * cr - cy * sr;
    m(2,0) = - sp;
    m(2,1) = cp * sr;
    m(2,2) = cp * cr;
    return m; 
}

template<typename T>
void mat2RPY(const Eigen::Matrix<T, 3, 3>& m, T& roll, T& pitch, T& yaw)
{
    roll = atan2(m(2,1), m(2,2));
    pitch = atan2(-m(2,0), sqrt(m(2,1) * m(2,1) + m(2,2) * m(2,2)));
    yaw = atan2(m(1,0), m(0,0));
}


void estimate_alltitude(const mjModel* m, mjData* d, double* roll, double* roll_vel, double* pitch, double* pitch_vel) {
    int chasis_id = 1;

    double x_angle = 0;
    double y_angle = 0;
    double z_angle = 0;

    Eigen::Quaterniond q(d->qpos[3], d->qpos[4], d->qpos[5], d->qpos[6]);
    Eigen::Matrix3d R = q.toRotationMatrix();
    mat2RPY<double>(R, x_angle, y_angle, z_angle);
    
    // std::cout << "chasis_quat: " << chasis_qx << " " << chasis_qy << " " << chasis_qz << " " << chasis_qw << std::endl;
    // std::cout << "Euler from quaternion in roll, pitch, yaw: " << euler[0] <<" " << euler[1] <<" " << euler[2] << std::endl;
//   // get chasis pitch and rate-of-pitch
    *roll = x_angle;
    *roll_vel = d->qvel[3];

    *pitch = y_angle;
    *pitch_vel = d->qvel[4];
    // std::cout << "vel: " << *roll_vel << " " << *pitch_vel << std::endl;
}

void estimate_wheel_vel(const mjModel* m, mjData* d, double* vel) {
    *vel = d->qvel[6]; //y
} 

double wheel_pos_x_vel_integration = 0;

void wheelcontrol(const mjModel* m, mjData* d, const double& pitch, const double& pitch_vel, const double& wheel_vel) {

  double Kp_chasis = 20;
  double Kd_chasis = 2;

  double Kp_wheel = 0.21;
  double Kd_wheel = 0.0;
  double Ki_wheel = 0.000100;

  double vel_x_setpoint = 1.51; // 0

  double vel_error = wheel_vel - vel_x_setpoint ;

  wheel_pos_x_vel_integration += vel_error;

  std::cout << "wheel_pos_x_vel_integration: " << wheel_pos_x_vel_integration << std::endl;

  d->ctrl[0] = - Kp_chasis * pitch - Kd_chasis * pitch_vel
    + Kp_wheel * vel_error + Ki_wheel * wheel_pos_x_vel_integration;

  std::cout << "d->ctrl[0]: " << d->ctrl[0] << std::endl;
  std::cout << "----------" << std::endl;
    
}

void pendulumcontrol(const mjModel* m, mjData* d, const double& roll, const double& roll_vel) {

  double Kp_chasis = 3;
  double Kd_chasis = 0;

//   d->qvel[6] = 0;
  
  std::cout << "roll " << roll << std::endl;
  d->ctrl[1] =  - (- Kp_chasis * roll - Kd_chasis * roll_vel);//

}



void mycontroller(const mjModel* m, mjData* d)
{

  int chasis_id = 1;


  // add noise 
  double noise;
  mju_standardNormal(&noise);
  d->xfrc_applied[6*chasis_id+0] = 1 * noise;
  d->xfrc_applied[6*chasis_id+1] = 4 * noise;

    double roll, pitch;
    double roll_vel, pitch_vel;
    estimate_alltitude(m, d, &roll, &roll_vel, &pitch, &pitch_vel);
    std::cout << "chasis_quat: " << roll << " " << roll_vel << " " << pitch << " " << pitch_vel << std::endl;
    double wheel_vel;
    estimate_wheel_vel(m,d, &wheel_vel);

    wheelcontrol(m,d, pitch, pitch_vel, wheel_vel);
    pendulumcontrol(m, d, roll, roll_vel);

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

    // install control callback
    mjcb_control = mycontroller;

    //m->opt.gravity[2]=-1;
    //qpos is dim nqx1 = 7x1; 3 translations + 4 quaternions
    // d->qpos[0]=0.21;
    // d->qvel[2]=5;
    // d->qvel[0]=2;
    // use the first while condition if you want to simulate for a period.

    // Eigen::Quaterniond q(1, 2, 0, 0);
    // q.normalized();
    // d->qpos[3] = q.w();
    // d->qpos[4] = q.x();
    // d->qpos[5] = q.y();
    // d->qpos[6] = q.z();

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

    //         double x_angle = 0;
    //         double y_angle = 0;
    // //   for (int i = 0; i < int(simend); i++) {
    //     // for (int j = 0; j < int(simend); j++) {
    //         x_angle = - M_PI / 2 +  d->time / (simend) * M_PI; 
    //         std::cout << "x: " << x_angle << std::endl;
    //         Eigen::Matrix3d R = RPY2mat<double>(x_angle, y_angle, 0);
    //         Eigen::Quaterniond q(R);
    //         d->qpos[3] = q.w();
    //         d->qpos[4] = q.x();
    //         d->qpos[5] = q.y();
    //         d->qpos[6] = q.z();
            // d->xquat[4 * chasis_id + 0] = q.w();
            // d->xquat[4 * chasis_id + 1] = q.x();
            // d->xquat[4 * chasis_id + 2] = q.y();
            // d->xquat[4 * chasis_id + 3] = q.z();
        // }
    //   }


            if (d->time>=simend)
            {
                break;
            }

        }

       // get framebuffer viewport
        mjrRect viewport = {0, 0, 0, 0};
        glfwGetFramebufferSize(window, &viewport.width, &viewport.height);

          // update scene and render
        // opt.frame = mjFRAME_WORLD;
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
