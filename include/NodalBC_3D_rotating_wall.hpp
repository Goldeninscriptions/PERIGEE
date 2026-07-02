#ifndef NODALBC_3D_ROTATING_WALL_HPP
#define NODALBC_3D_ROTATING_WALL_HPP
// ============================================================================
// NodalBC_3D_rotating_wall.hpp
//
// Collect nodal and surface data on a rotating wall for strong velocity
// boundary-condition updates in solver time stepping.
// ============================================================================
#include "INodalBC.hpp"
#include "FEType.hpp"
#include "Tet_Tools.hpp"
#include "Hex_Tools.hpp"

class NodalBC_3D_rotating_wall : public INodalBC
{
  public:
    NodalBC_3D_rotating_wall(
        const std::string &wall_file,
        const int &nFunc,
        const FEType &in_elemtype );

    virtual ~NodalBC_3D_rotating_wall() = default;

    virtual unsigned int get_dir_nodes(const unsigned int &ii) const
    { return dir_nodes[ii]; }

    virtual unsigned int get_num_dir_nodes() const
    { return static_cast<unsigned int>(dir_nodes.size()); }

    virtual unsigned int get_per_master_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: NodalBC_3D_rotating_wall::get_per_master_nodes: periodic nodes are not defined.\n");
      return 0;
    }

    virtual unsigned int get_per_slave_nodes(const unsigned int &ii) const
    {
      SYS_T::print_fatal("Error: NodalBC_3D_rotating_wall::get_per_slave_nodes: periodic nodes are not defined.\n");
      return 0;
    }

    virtual unsigned int get_num_per_nodes() const { return 0; }

    virtual unsigned int get_dir_nodes_on_rotated_surface( const unsigned int &ii ) const
    { return dir_nodes_on_rotated_surface[ii]; }

    virtual unsigned int get_num_dir_nodes_on_rotated_surface() const
    { return static_cast<unsigned int>(dir_nodes_on_rotated_surface.size()); }

    virtual int get_num_node() const { return num_node; }
    virtual int get_num_cell() const { return num_cell; }
    virtual int get_nLocBas() const { return nLocBas; }

    virtual int get_ien(const int &cell, const int &lnode) const
    { return sur_ien[nLocBas * cell + lnode]; }

    virtual double get_pt_xyz(const int &node, const int &dir) const
    { return pt_xyz[3 * node + dir]; }

    virtual int get_global_node(const int &node_idx) const
    { return global_node[node_idx]; }

    virtual int get_global_cell(const int &cell_idx) const
    { return global_cell[cell_idx]; }

  private:
    NodalBC_3D_rotating_wall() = delete;

    std::vector<unsigned int> dir_nodes_on_rotated_surface {};
    std::vector<unsigned int> dir_nodes {};

    const FEType elem_type;

    int num_node {0}, num_cell {0}, nLocBas {0};
    std::vector<int> sur_ien {};
    std::vector<double> pt_xyz {};
    std::vector<int> global_node {};
    std::vector<int> global_cell {};
};

#endif
