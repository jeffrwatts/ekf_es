
#include <BasicLinearAlgebra.h>
using namespace BLA;


BLA::Matrix<1,3> imu_f0 = {-0.01996148,  0.03136036,  9.78135591};
BLA::Matrix<1,3> imu_f1 = {-0.01986699,  0.03743271,  9.79679338};
BLA::Matrix<1,4> q = {0.6800922, 0.3398297, 0.2268417, 0.6087144};
BLA::Matrix<1,3> g = {0, 0, -9.81};
float delta_t = 0.004999999999999893;

BLA::Matrix<3,3> I = {1, 0, 0,
                      0, 1, 0,
                      0, 0, 1};

BLA::Matrix<3,3> skew_symmetric(BLA::Matrix<3,1> v){
  return BLA::Matrix<3,3> {0, -v(2), v(1),
                          v(2), 0, -v(0),
                          -v(1), v(0), 0};
}

BLA::Matrix<3,3> quaternion_to_matrix(BLA::Matrix<1,4> q) {
  float w = q(0,0);
  BLA::Matrix<3,1> v = { q(0,1), q(0,2), q(0,3) };
  return (I * Matrix<1,1>(Matrix<1,1>(pow(w, 2)) - ~v*v)(0,0)) + ((v * ~v) * 2) + skew_symmetric(v) * w * 2;
}

void setup() {
  Serial.begin(115200);
  while ( !Serial ) delay(10);   
    
  // put your setup code here, to run once:
  BLA::Matrix<3,3> Cns = quaternion_to_matrix(q);
  Serial << Cns;

  float p_est0 = 0;
  float v_est0 = 0;

  //Serial << p_est0 + delta_t * v_est0 + (pow(delta_t,2) / 2) * (Cns * imu_f0 + g)
}

void loop() {
  // put your main code here, to run repeatedly:

}
