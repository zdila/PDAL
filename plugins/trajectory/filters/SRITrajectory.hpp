#if !defined(SRI_TRAJECTORY_HPP)
#define SRI_TRAJECTORY_HPP 1

#include <pdal/Filter.hpp>

namespace pdal {

  struct PrivateArgs;

  class PDAL_DLL SRITrajectory : public Filter {
  public:
    SRITrajectory();
    std::string getName() const;

  private:
    std::unique_ptr<PrivateArgs> m_args;

    virtual void addArgs(ProgramArgs& args);
    virtual void addDimensions(PointLayoutPtr layout);
    virtual PointViewSet run(PointViewPtr view);
    void prepared(PointTableRef table);

    SRITrajectory& operator=(const SRITrajectory&); // not implemented
    SRITrajectory(const SRITrajectory&); // not implemented
  };

} // namespace pdal

#endif
