#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <LidarTrajectory/PulseCollection.hpp>

int main() {
  try {
    LidarTrajectory::PulseCollection coll(std::cin);
    {
      std::ifstream input(coll.params.lookup<std::string>("input").c_str());
      std::string line;
      std::istringstream str;
      Eigen::Vector3d r;
      while (std::getline(input, line)) {
        str.clear(); str.str(line);
        int p, ret, num;
        double t, ang;
        if (!(str >> p >> t
              >> r(0) >> r(1) >> r(2)
              >> ang >> num >> ret))
          throw std::runtime_error("Short read on " + line);
        coll.Add(t, r, num, ret, ang);
      }
    }
    coll.Finalize();
    {
      double
        tmin = coll.tmin,
        tmax = coll.tmax;
      std::ifstream
        ground(coll.params.lookup<std::string>("groundtruth").c_str());
      std::ofstream
        output(coll.params.lookup<std::string>("output").c_str());
      std::istringstream str;
      std::string line;
      Eigen::Vector3d r;
      int count = 0;
      double t, herr = 0, verr = 0;
      while (std::getline(ground, line)) {
        str.clear(); str.str(line);
        str >> t >> r(0) >> r(1) >> r(2);
        if (t < tmin || t > tmax) continue;
        ++count;
        Eigen::Vector3d v, a, p = coll.Trajectory(t, v, a);
        output << std::fixed << std::setprecision(2)
               << t << " " << std::setprecision(4)
               << r(0) << " " << r(1) << " " << r(2) << " "
               << p(0) << " " << p(1) << " " << p(2) << " "
               << v(0) << " " << v(1) << " " << v(2) << " "
               << a(0) << " " << a(1) << " " << a(2) << "\n";
        r = p - r;
        herr += r(0)*r(0) + r(1)*r(1);
        verr += r(2)*r(2);
      }
      herr = std::sqrt(herr/count);
      verr = std::sqrt(verr/count);
      std::cout << std::fixed << std::setprecision(3)
                << "horiz error " << herr << "; vert error " << verr << "\n";
    }
  }
  catch (const std::exception& e) {
    std::cerr << "Caught exception: " << e.what() << "\n";
    return 1;
  }
  catch (...) {
    std::cerr << "Caught unknown exception\n";
    return 1;
  }

}
