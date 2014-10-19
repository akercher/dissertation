/* data_io.h */

/*****************************************************/
/* Output states                                     */
/*---------------------------------------------------*/
/*****************************************************/
#ifdef MHD
void output_vtk_legacy(std::ofstream& output, Mesh mesh, Real gamma, 
		       StateArray state_array, RealArray current)
#else
void output_vtk_legacy(std::ofstream& output, Mesh mesh, Real gamma, StateArray state_array)
#endif
{

  Index i,j,k;
  State state;
  State prim_state;

  output << "# vtk DataFile Version 4.2" <<"\n";
#ifdef MHD
  output << "Ideal MHD" <<"\n";
#else
  output << "Compressible Euler" <<"\n";
#endif
  output << "ASCII" <<"\n";
  output << " " <<"\n";
  output << "DATASET UNSTRUCTURED_GRID" <<"\n";
  output << "POINTS " << mesh.npoin() << " DOUBLE"<<"\n";
  /* output << ncell_x << " " << ncell_y <<"\n"; */


  for(j = 0; j < mesh.ny ; j++){
    for(i = 0; i < mesh.nx; i++){

      k = mesh.nx*j + i;
      output << mesh.dx*Real(i) << " " << mesh.dy*Real(j) << " " << Real(0.0) <<"\n";

    }
  }
  output << " " <<"\n";

  output << "CELLS " << Index(2)*mesh.ncell() << " " << Index(4)*Index(2)*mesh.ncell() <<"\n";
  for(j = 0; j < mesh.ncell_y ; j++){
    for(i = 0; i < mesh.ncell_x; i++){

      k = mesh.nx*j + i;
      output << Index(3) << " " << k << " " << k + Index(1) << " " << k + mesh.nx + Index(1) <<"\n";
      output << Index(3) << " " << k << " " << k + mesh.nx + Index(1) << " " << k + mesh.nx <<"\n";

    }
  }
  output << " " <<"\n";
  output << "CELL_TYPES " << Index(2)*mesh.ncell() <<"\n";
  for(i = 0; i < Index(2)*mesh.ncell_x*mesh.ncell_y; i++){
    output << Index(5) <<"\n";
  }
  output << " " <<"\n";

  output << "POINT_DATA " << mesh.npoin() <<"\n";
  output << "SCALARS density DOUBLE" <<"\n";
  output << "LOOKUP_TABLE default" <<"\n";
  for(j = 0; j < mesh.ny ; j++){
    for(i = 0; i < mesh.nx; i++){

  	k = mesh.nx*j + i;

  	state = state_array[k];
  	output << density << "\n";

    }
  }

  output << "SCALARS pressure DOUBLE" <<"\n";
  output << "LOOKUP_TABLE default" <<"\n";
  for(j = 0; j < mesh.ny ; j++){
    for(i = 0; i < mesh.nx; i++){

  	k = mesh.nx*j + i;

  	state = state_array[k];
  	state = prim2cons_func(gamma,State(state_array[k]));
  	output << pressure << "\n";

    }
  }

  output << "VECTORS velocity DOUBLE" <<"\n";
  for(j = 0; j < mesh.ny ; j++){
    for(i = 0; i < mesh.nx; i++){

  	k = mesh.nx*j + i;

  	/* state = state_array[k]; */
	state = prim2cons_func(gamma,State(state_array[k]));

  	output << get_x(velocity) << " " << get_y(velocity)
  	       << " " << get_z(velocity) << "\n";

    }
  }

#ifdef MHD
  output << "VECTORS magnetic DOUBLE" <<"\n";
  for(j = 0; j < mesh.ny ; j++){
    for(i = 0; i < mesh.nx; i++){

  	k = mesh.nx*j + i;

  	/* state = state_array[k]; */
	state = prim2cons_func(gamma,State(state_array[k]));

  	output << get_x(bfield) << " " << get_y(bfield)
  	       << " " << get_z(bfield) << "\n";

    }
  }

  output << "SCALARS current DOUBLE" <<"\n";
  output << "LOOKUP_TABLE default" <<"\n";
  for(j = 0; j < mesh.ny ; j++){
    for(i = 0; i < mesh.nx; i++){

  	k = mesh.nx*j + i;
  	output << current[k] << "\n";

    }
  }

#endif

  output << " " <<"\n";
  output.close();

};
