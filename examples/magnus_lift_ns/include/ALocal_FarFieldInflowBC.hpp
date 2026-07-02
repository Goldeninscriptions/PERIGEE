#ifndef ALOCAL_FARFIELD_INFLOWBC_HPP
#define ALOCAL_FARFIELD_INFLOWBC_HPP

#include "ALocal_InflowBC.hpp"

class ALocal_FarFieldInflowBC : public ALocal_InflowBC
{
  public:
    ALocal_FarFieldInflowBC( const std::string &fileBaseName, int cpu_rank );

    ALocal_FarFieldInflowBC( const HDF5_Reader * const &h5r );

    virtual ~ALocal_FarFieldInflowBC() = default;

    double get_radius( int nbc_id, const Vector_3 &pt ) const override;
};

#endif
