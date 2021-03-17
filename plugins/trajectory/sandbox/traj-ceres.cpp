// Input file is a text file produced by
// fields=PointSourceId:0,GpsTime:6,X:2,Y:2,Z:2,ScanAngleRank:0,NumberOfReturns:0,ReturnNumber:0
// pdal translate
//   --writers.text.order="$fields"
//   --writers.text.write_header=false
//   --writers.text.keep_unspecified=false
//   --writers.text.delimiter=' '
//   C2_L2.las C2_L2.txt

// Run with

// ./traj-ceres < C2_L2.txt > traj.txt

// Columns of traj.txt are
// time xest yest zest xfit yfit zfit

// {x,y,z}est is estimate of plane's position (no smoothing)

// ignore {x,y,z}fit.  This is a polynomial fit to the whole trajectory which
// is not what we want.

#include "lasertraj.hpp"

int main() {
  lasertraj traj(std::cin);
  traj.estimate();
  traj.initialfit();
  if (false) {
    for (double t = 8; t <= 77; ++t) {
      Eigen::Matrix<double, lasertraj::ord, 3> C = traj.ceresfit(t, true);
      for (int i = -10; i <= 10; ++i) {
        double t1 = t + double(i)/10;
        Eigen::Vector3d r = traj.rorg + lasertraj::cubic3(C, t1),
          v = lasertraj::cubic3d(C, t1);
        std::cout << i << " " << std::fixed << std::setprecision(3)
                  << traj.torg + t1 << " "
                  << r(0) << " " << r(1) << " " << r(2) << " "
                  << v(0) << " " << v(1) << " " << v(2) << "\n";
      }
    }
  }
  traj.finalfit();
  if (false) {
    std::cout << std::fixed;
    for (int i = 71; i <= 780; ++i) {
      double t = 0.1 * i;
      auto p = traj.coords(t);
      auto r = p.first + traj.rorg;
      auto v = p.second;
      std::cout << std::setprecision(1) << t + traj.torg << " "
                << std::fixed << std::setprecision(3)
                << r(0) << " " << r(1) << " " << r(2) << " "
                << v(0) << " " << v(1) << " " << v(2) << "\n";
    }
  } else {
    traj.attitudefitfull();
    std::cout << std::fixed;
    for (int i = 71; i <= 780; ++i) {
      double t = 0.1 * i;
      auto p = traj.coords(t);
      auto r = p.first + traj.rorg;
      auto v = p.second;
      auto a = traj.att(t).first;
      std::cout << std::setprecision(1) << t + traj.torg << " "
                << std::fixed << std::setprecision(3)
                << r(0) << " " << r(1) << " " << r(2) << " "
                << v(0) << " " << v(1) << " " << v(2) << " "
                << a(0) << " " << a(1) << " " << a(2) << "\n";
    }
  }
  //  traj.print(std::cout);
}
