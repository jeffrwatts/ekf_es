
#include <BasicLinearAlgebra.h>
using namespace BLA;

float var_imu_f = 0.10;
float var_imu_w = 0.25;

BLA::Matrix<3,1> g = {0, 0, -9.81};

BLA::Matrix<3,3> I3 = {1, 0, 0,
                       0, 1, 0,
                       0, 0, 1};

BLA::Matrix<4,4> I4 = {1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       0, 0, 0, 1};

BLA::Matrix<9,9> I9 = {1, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 1, 0, 0, 0 ,0, 0, 0, 0,
                       0, 0, 1, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 1, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 1, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 1, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 1, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 1, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 1};

BLA::Matrix<3,3> skew_symmetric(BLA::Matrix<3,1> v){
  return BLA::Matrix<3,3> {0, -v(2), v(1),
                          v(2), 0, -v(0),
                          -v(1), v(0), 0};
}

BLA::Matrix<3,3> quaternion_to_matrix(BLA::Matrix<4,1> q) {
  float w = q(0,0);
  BLA::Matrix<3,1> v = { q(1,0), q(2,0), q(3,0) };
  return (I3 * Matrix<1,1>(Matrix<1,1>(pow(w, 2)) - ~v*v)(0,0)) + ((v * ~v) * 2) + skew_symmetric(v) * w * 2;
}

BLA::Matrix<4,1> euler_to_quaternion(BLA::Matrix<3,1> euler) {
  float roll = euler(0,0);
  float pitch = euler(1,0);
  float yaw = euler(2,0);
  
  float cy = cos(yaw * 0.5);
  float sy = sin(yaw * 0.5);
  float cr = cos(roll * 0.5);
  float sr = sin(roll * 0.5);
  float cp = cos(pitch * 0.5);
  float sp = sin(pitch * 0.5);
  
  return Matrix<4,1> {cr * cp * cy + sr * sp * sy,
                      sr * cp * cy - cr * sp * sy,
                      cr * sp * cy + sr * cp * sy,
                      cr * cp * sy - sr * sp * cy};
}

//BLA::Matrix<4,1> axis_angle_to_quaternion(BLA::Matrix<3,1> axis_angle) {
//  float x = axis_angle(0,0);
//  float y = axis_angle(1,0);
//  float z = axis_angle(2,0);

//  float c1 = cos(y/2);
//  float c2 = cos(z/2);
//  float c3 = cos(x/2);

//  float s1 = sin(y/2);
//  float s2 = sin(z/2);
//  float s3 = sin(x/2);

//  Matrix<4,1> quat = {c1 * c2 * c3 - s1 * s2 * s3,
//                      s1 * s2 * c3 + c1 * c2 * s3,
//                      s1 * c2 * c3 + c1 * s2 * s3,
//                      c1 * s2 * c3 - s1 * c2 * s3};
//  return quat;
//}

BLA::Matrix<4,1> quaternion_multiply(BLA::Matrix<4,1> q_left, BLA::Matrix<4,1> q_right) {
  BLA::Matrix<3,1> v_right = { q_right(1,0), q_right(2,0), q_right(3,0) };
  BLA::Matrix<3,3> v_right_skew = -skew_symmetric(v_right);
  BLA::Matrix<4,4> sum_term = {0, -v_right(0,0), -v_right(1,0), -v_right(2,0),
                               v_right(0,0), v_right_skew(0,0), v_right_skew(0,1), v_right_skew(0,2),
                               v_right(1,0), v_right_skew(1,0), v_right_skew(1,1), v_right_skew(1,2),
                               v_right(2,0), v_right_skew(2,0), v_right_skew(2,1), v_right_skew(2,2)};
  BLA::Matrix<4,4> sigma = I4 * q_right(0,0) + sum_term;
  return sigma * q_left;
}

BLA::Matrix<9,6> L = {0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0,
                      1, 0, 0, 0, 0, 0,
                      0, 1, 0, 0, 0, 0,
                      0, 0, 1, 0, 0, 0,
                      0, 0, 0, 1, 0, 0,
                      0, 0, 0, 0, 1, 0,
                      0, 0, 0, 0, 0, 1};

BLA::Matrix<6,6> Q = {var_imu_f, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, var_imu_f, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, var_imu_f, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, var_imu_w, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, var_imu_w, 0.0, 
                      0.0, 0.0, 0.0, 0.0, 0.0, var_imu_w};


BLA::Matrix<9,9> compute_F(Matrix<3,1> accel, float delta_t) {
  BLA::Matrix<3,3> skew = -skew_symmetric(accel)*delta_t;

  return Matrix<9,9> {1, 0, 0, delta_t, 0, 0, 0, 0, 0,
                      0, 1, 0, 0, delta_t, 0, 0, 0, 0,
                      0, 0, 1, 0, 0, delta_t, 0, 0, 0,
                      0, 0, 0, 1, 0, 0, skew(0,0), skew(0, 1), skew(0, 2),
                      0, 0, 0, 0, 1, 0, skew(1,0), skew(1, 1), skew(1, 2),
                      0, 0, 0, 0, 0, 1, skew(2,0), skew(2, 1), skew(2, 2),
                      0, 0, 0, 0, 0, 0, 1, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 1, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 1};
}

void setup() {
  Serial.begin(115200);
  while ( !Serial ) delay(10);  

  float delta_t = 0.005;
  BLA::Matrix<3,1> pos_prev = {-176.86869371, -107.70275048, 8.13058728};
  BLA::Matrix<3,1> vel_prev = {-22.88698488, -25.38412216, -1.20899051};
  BLA::Matrix<4,1> pose_prev = {0.76350596, 0.05554539, -0.04255976, -0.64199846};
  BLA::Matrix<9,9> p_cov_prev = {1.88392327e+06, -2.28713810e+03,  9.08852141e+04,  9.39861384e+04, -2.17170742e+02,  6.44979567e+03,  0.00000000e+00,  2.55562854e+02, 1.01600236e+01,
                                 -2.28713810e+03,  1.88741617e+06,  4.71358351e+04, -3.24780716e+02, 9.40631853e+04,  6.88874705e+03, -2.55562854e+02,  0.00000000e+00, -1.20143304e+01,
                                 9.08852141e+04,  4.71358351e+04,  8.42123706e+03,  4.49952675e+03, 2.83624926e+03,  6.50058899e+02, -1.01600236e+01,  1.20143304e+01, 0.00000000e+00,
                                 9.39861384e+04, -3.24780716e+02,  4.49952675e+03,  5.00981878e+03, -3.07531265e+01,  3.53706842e+02,  0.00000000e+00,  1.52408509e+01, 1.45193134e+00,
                                 -2.17170742e+02,  9.40631853e+04,  2.83624926e+03, -3.07531265e+01, 5.00028396e+03,  4.03827575e+02, -1.52408509e+01,  0.00000000e+00, -1.20934971e+00,
                                 6.44979567e+03,  6.88874705e+03,  6.50058899e+02,  3.53706842e+02, 4.03827575e+02,  6.74640483e+01, -1.45193134e+00,  1.20934971e+00, 0.00000000e+00,
                                 0.00000000e+00, -2.55562854e+02, -1.01600236e+01,  0.00000000e+00, -1.52408509e+01, -1.45193134e+00,  6.25000000e-02,  0.00000000e+00, 0.00000000e+00,
                                 2.55562854e+02,  0.00000000e+00,  1.20143304e+01,  1.52408509e+01, 0.00000000e+00,  1.20934971e+00,  0.00000000e+00,  6.25000000e-02, 0.00000000e+00,
                                 1.01600236e+01, -1.20143304e+01,  0.00000000e+00,  1.45193134e+00, -1.20934971e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, 6.25000000e-02};

  BLA::Matrix<3,1> imu_f = {0.73546694, 0.53162613, 9.74813043};
  BLA::Matrix<3,1> imu_w = {-0.19374949, -0.085177, 0.24635924};

  // 1. Use motion model to update state with IMU inputs.
  BLA::Matrix<3,1> accel = quaternion_to_matrix(pose_prev) * imu_f + g;   
  BLA::Matrix<3,1> pos = pos_prev + vel_prev * delta_t + accel * pow(delta_t, 2) / 2;
  BLA::Matrix<3,1> vel = vel_prev + accel * delta_t;

  BLA::Matrix<4,1> q_imu = euler_to_quaternion(imu_w * delta_t);
  Matrix<4,1> pose = quaternion_multiply(pose_prev, q_imu); 

  // 2. Linearize the motion model and compute Jacobians.
  BLA::Matrix<9,9> F = compute_F(accel-g, delta_t);
  Serial.println("F");
  Serial << F;
  Serial.println("Q");
  Serial << Q * pow(delta_t,2);

  BLA::Matrix<9,9> p_cov = F * p_cov_prev * ~F + L * (Q * pow(delta_t,2)) * ~L;
  Serial << p_cov;


  
  //Serial << pos;
  //Serial.println("========");
  //Serial << vel;
}

void loop() {

}
