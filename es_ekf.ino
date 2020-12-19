#include <SPI.h>
#include <Wire.h>
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BNO055.h>
#include <utility/imumaths.h>
#include <BasicLinearAlgebra.h>
using namespace BLA;

Adafruit_SSD1306 display = Adafruit_SSD1306(128, 32, &Wire);
Adafruit_BNO055 bno = Adafruit_BNO055(55);

Adafruit_BNO055::adafruit_vector_type_t types[] = {
  Adafruit_BNO055::VECTOR_GRAVITY,
  Adafruit_BNO055::VECTOR_EULER,
  Adafruit_BNO055::VECTOR_LINEARACCEL,
  Adafruit_BNO055::VECTOR_ACCELEROMETER,
  Adafruit_BNO055::VECTOR_MAGNETOMETER,
  Adafruit_BNO055::VECTOR_GYROSCOPE};
 
#define BUTTON_A  9
#define BUTTON_B  6
#define BUTTON_C  5

uint16_t BNO055_SAMPLERATE_DELAY_MS = 50;

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

void displayCalibration(uint8_t sys, uint8_t gyro, uint8_t accel, uint8_t mag) {
  display.setCursor(0,0);
  display.clearDisplay();
  display.print("Sys: ");
  display.println(sys);
  display.print("Gyro: ");
  display.println(gyro);
  display.print("Accel: ");
  display.println(accel);
  display.print("Mag: ");
  display.println(mag);
  display.display();
}

void displayEvent(sensors_event_t* event) {
  display.setCursor(0,0);
  display.clearDisplay();
  double x = -1000000, y = -1000000 , z = -1000000; //dumb values, easy to spot problem

  Serial.print("Version: ");
  Serial.print(event->version);
  Serial.print("; Sensor ID: ");
  Serial.print(event->sensor_id);
  Serial.print("; Type: ");
  Serial.println(event->type);
  
  if (event->type == SENSOR_TYPE_LINEAR_ACCELERATION) {
    // From Adafruit_BNO055::VECTOR_LINEARACCEL: Acceleration not including gravity.
    display.println("LINEAR_ACCELERATION");
    x = event->acceleration.x;
    y = event->acceleration.y;
    z = event->acceleration.z;
  }
  else if (event->type == SENSOR_TYPE_ACCELEROMETER) {
    // From Adafruit_BNO055::VECTOR_ACCELEROMETER or Adafruit_BNO055::VECTOR_GRAVITY
    display.println("ACCELEROMETER");
    x = event->acceleration.x;
    y = event->acceleration.y;
    z = event->acceleration.z;
  }
  else if (event->type == SENSOR_TYPE_ORIENTATION) {
    // From Adafruit_BNO055::VECTOR_EULER
    display.println("ORIENTATION");
    x = event->orientation.x;
    y = event->orientation.y;
    z = event->orientation.z;
  }
  else if (event->type == SENSOR_TYPE_ROTATION_VECTOR) {
    // From Adafruit_BNO055::VECTOR_GYROSCOPE
    display.println("ROTATION_VECTOR");
    x = event->gyro.x;
    y = event->gyro.y;
    z = event->gyro.z;
  }
  else if (event->type == SENSOR_TYPE_MAGNETIC_FIELD) {
    // From Adafruit_BNO055::VECTOR_MAGNETOMETER
    display.println("MAGNETIC_FIELD");
    x = event->magnetic.x;
    y = event->magnetic.y;
    z = event->magnetic.z;
  }
  else {
    display.println("UNKNOWN");
  }

  display.print("x: ");
  display.println(x);
  display.print("y: ");
  display.println(y);
  display.print("z: ");
  display.println(z);
  display.display();
}

boolean calibrate = true;
boolean is_running = false;
unsigned long t = millis();

BLA::Matrix<3,1> pos_prev = {0, 0, 0};
BLA::Matrix<3,1> vel_prev = {0, 0, 0};
BLA::Matrix<4,1> pose_prev = {0, 0, 0};

BLA::Matrix<3,1> pos = {0, 0, 0};
BLA::Matrix<3,1> vel = {0, 0, 0};
BLA::Matrix<4,1> pose = {0, 0, 0};

void setup() {
  Serial.begin(9600);
  Serial.println("Setup");
 
  // SSD1306_SWITCHCAPVCC = generate display voltage from 3.3V internally
  display.begin(SSD1306_SWITCHCAPVCC, 0x3C); // Address 0x3C for 128x32

  pinMode(BUTTON_A, INPUT_PULLUP);
  pinMode(BUTTON_B, INPUT_PULLUP);
  pinMode(BUTTON_C, INPUT_PULLUP);
 
  // Display splash
  display.display();
  delay(1000);

  // Reset for text.
  display.setTextSize(1);
  display.setTextColor(SSD1306_WHITE);

  display.display(); // actually display all of the above

  boolean bno_initialized = bno.begin();
  bno.setExtCrystalUse(true);
}

void loop() {
  if(!digitalRead(BUTTON_A)) {
    calibrate = false;
    display.setCursor(0,0);
    display.clearDisplay();
    if (is_running == true) {
      is_running = false;
      display.println ("Stopping...");
    } else {
      is_running = true;
      pos = {0, 0, 0};
      vel = {0, 0, 0};
      pose = {0, 0, 0, 0};
      pos_prev = {0, 0, 0};
      vel_prev = {0, 0, 0};
      pose_prev = {1, 0, 0, 0};
      display.println ("Starting...");
    }
    display.display();
    delay(1000);
    t = millis();
  }

  if (is_running == true) {
    float delta_t = (millis() - t) / 1000.0;
    t = millis();
    sensors_event_t sensor_imu_f;
    bno.getEvent(&sensor_imu_f, Adafruit_BNO055::VECTOR_LINEARACCEL);
    Matrix<3,1> imu_f = {sensor_imu_f.acceleration.x, sensor_imu_f.acceleration.y, sensor_imu_f.acceleration.z};

    sensors_event_t sensor_imu_w; 
    bno.getEvent(&sensor_imu_w, Adafruit_BNO055::VECTOR_GYROSCOPE);
    Matrix<3,1> imu_w = {sensor_imu_w.gyro.x, sensor_imu_w.gyro.y, sensor_imu_w.gyro.z}; 

    //delta_t = 0.005;
    //pos_prev = {-176.86869371, -107.70275048, 8.13058728};
    //vel_prev = {-22.88698488, -25.38412216, -1.20899051};
    //pose_prev = {0.76350596, 0.05554539, -0.04255976, -0.64199846};
    //imu_f = {0.73546694, 0.53162613, 9.74813043};
    //imu_w = {-0.19374949, -0.085177, 0.24635924};

    // 1. Use motion model to update state with IMU inputs.
    BLA::Matrix<3,1> accel = quaternion_to_matrix(pose_prev) * imu_f;   
    
    pos = pos_prev + vel_prev * delta_t + accel * pow(delta_t, 2) / 2;
    vel = vel_prev + accel * delta_t;
    pose = quaternion_multiply(pose_prev, euler_to_quaternion(imu_w * delta_t)); 

    pos_prev = pos;
    vel_prev = vel;
    pose_prev = pose;
  }

  if (calibrate == false) {
    display.setCursor(0,0);
    display.clearDisplay();
    display.print("x: ");
    display.println(pos(0));
    display.print("y: ");
    display.println(pos(1)); 
    display.print("z: ");
    display.println(pos(2));  
    display.display();
  } else {
    uint8_t sys, gyro, accel, mag = 0;
    bno.getCalibration(&sys, &gyro, &accel, &mag);
    displayCalibration(sys, gyro, accel, mag);
  }
  
  delay(BNO055_SAMPLERATE_DELAY_MS);
}

void TestMotionModel() {
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
  BLA::Matrix<9,9> p_cov = F * p_cov_prev * ~F + L * (Q * pow(delta_t,2)) * ~L;
}
