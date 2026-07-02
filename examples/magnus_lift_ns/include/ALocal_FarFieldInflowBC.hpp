#ifndef ALOCAL_FARFIELD_INFLOWBC_HPP
#define ALOCAL_FARFIELD_INFLOWBC_HPP

#include "Vec_Tools.hpp"
#include "Vector_3.hpp"

class HDF5_Reader;

class ALocal_FarFieldInflowBC
{
  public:
    ALocal_FarFieldInflowBC( const std::string &fileBaseName, int cpu_rank );

    ALocal_FarFieldInflowBC( const HDF5_Reader * const &h5r );

    virtual ~ALocal_FarFieldInflowBC() = default;

    virtual int get_num_nbc() const { return num_nbc; }

    virtual int get_LDN( int nbc_id, int node ) const
    { return LDN[nbc_id][node]; }

    virtual int get_Num_LD( int nbc_id ) const { return Num_LD[nbc_id]; }

    virtual double get_outvec( int nbc_id, int ii ) const
    { return outnormal[nbc_id](ii); }

    virtual Vector_3 get_outvec( int nbc_id ) const
    { return outnormal[nbc_id]; }

    virtual double get_actarea( int nbc_id ) const
    { return act_area[nbc_id]; }

    virtual double get_fularea( int nbc_id ) const
    { return ful_area[nbc_id]; }

    virtual bool is_inLDN( int nbc_id, int ii ) const
    { return VEC_T::is_invec(LDN[nbc_id], ii); }

    virtual double get_radius( int nbc_id, const Vector_3 &pt ) const;

    virtual int get_num_local_node( int nbc_id ) const
    { return num_local_node[nbc_id]; }

    virtual int get_num_local_cell( int nbc_id ) const
    { return num_local_cell[nbc_id]; }

    virtual int get_cell_nLocBas( int nbc_id ) const
    { return cell_nLocBas[nbc_id]; }

    virtual Vector_3 get_local_pt_xyz( int nbc_id, int ii ) const
    { return local_pt_xyz[nbc_id][ii]; }

    virtual int get_local_cell_ien( int nbc_id, int ii ) const
    { return local_cell_ien[nbc_id][ii]; }

    virtual int get_local_cell_ien( int nbc_id, int ee, int ii ) const
    { return local_cell_ien[nbc_id][ee * cell_nLocBas[nbc_id] + ii]; }

    virtual void get_ctrlPts_xyz( int nbc_id, int eindex,
        double * const &ctrl_x, double * const &ctrl_y,
        double * const &ctrl_z ) const;

    virtual void get_SIEN( int nbc_id, int eindex,
        int * const &sien ) const;

    virtual std::vector<int> get_SIEN( int nbc_id, int eindex ) const;

  private:
    int num_nbc;

    std::vector<int> Num_LD;
    std::vector< std::vector<int> > LDN;
    std::vector<Vector_3> outnormal;
    std::vector<double> act_area;
    std::vector<double> ful_area;
    std::vector<int> num_out_bc_pts;
    std::vector< std::vector<double> > outline_pts;
    std::vector<Vector_3> centroid;
    std::vector<int> num_local_node, num_local_cell, cell_nLocBas;
    std::vector< std::vector<Vector_3> > local_pt_xyz;
    std::vector< std::vector<int> > local_cell_ien;
    std::vector< std::vector<int> > local_node_pos;

    ALocal_FarFieldInflowBC() = delete;
};

#endif
