#include "NodalBC_3D_rotating_wall.hpp"
#include "Vec_Tools.hpp"
#include "VTK_Tools.hpp"

NodalBC_3D_rotating_wall::NodalBC_3D_rotating_wall(
    const std::string &wall_file,
    const int &nFunc,
    const FEType &in_elemtype )
: elem_type( in_elemtype )
{
  dir_nodes_on_rotated_surface.clear();
  dir_nodes.clear();

  SYS_T::file_check( wall_file );

  if( elem_type == FEType::Tet4 )
    nLocBas = 3;
  else if( elem_type == FEType::Tet10 )
    nLocBas = 6;
  else if( elem_type == FEType::Hex8 )
    nLocBas = 4;
  else if( elem_type == FEType::Hex27 )
    nLocBas = 9;
  else
    SYS_T::print_fatal("Error: NodalBC_3D_rotating_wall unknown element type.\n");

  VTK_T::read_grid( wall_file, num_node, num_cell, pt_xyz, sur_ien );

  global_node = VTK_T::read_int_PointData( wall_file, "GlobalNodeID" );
  global_cell = VTK_T::read_int_CellData( wall_file, "GlobalElementID" );

  for( unsigned int jj=0; jj<global_node.size(); ++jj )
  {
    SYS_T::print_fatal_if( global_node[jj] < 0 || global_node[jj] >= nFunc,
        "Error: NodalBC_3D_rotating_wall nodal index %d is not in [0, %d).\n",
        global_node[jj], nFunc );

    dir_nodes_on_rotated_surface.push_back( static_cast<unsigned int>(global_node[jj]) );
    dir_nodes.push_back( static_cast<unsigned int>(global_node[jj]) );
  }

  const auto num_dir_nodes_before_dedup = dir_nodes_on_rotated_surface.size();
  VEC_T::sort_unique_resize( dir_nodes_on_rotated_surface );
  VEC_T::sort_unique_resize( dir_nodes );

  SYS_T::print_fatal_if( num_dir_nodes_before_dedup != dir_nodes_on_rotated_surface.size(),
      "Error: NodalBC_3D_rotating_wall repeated nodes detected in %s.\n",
      wall_file.c_str() );

  Create_ID( nFunc );

  std::cout<<"===> NodalBC_3D_rotating_wall specified by\n";
  std::cout<<"     "<<wall_file<<" :\n";
  std::cout<<"          num_node: "<<num_node<<", num_cell: "<<num_cell<<'\n';
}

