/* data_io.h */

/*****************************************************/
/* Output states                                     */
/*---------------------------------------------------*/
/*****************************************************/
void output_vtk_legacy(std::ofstream& output, Mesh mesh, Real gamma, StateArray state_array)
{

  Index i,j,k;
  State state;
  State prim_state;
  /* char fname[240]; */
  /* FILE * fid; */

  /*  sprintf(fname,"/home/akercher/git/cfd/fvedge/dat/kh_instability_1.vtk"); */
  /* fid = fopen(fname,"w"); */

  /* fprintf(fid,"# vtk DataFile Version 4.2\n"); */
  /* fprintf(fid,"Incompressible Euler\n"); */
  /* fprintf(fid,"ASCII\n"); */
  /* fprintf(fid,"DATASET UNSTRUCTURED_GRID\n"); */
  /* fprintf(fid,"POINTS %d DOUBLE\n",mesh.npoin()); */
  /* for(j = 0; j < mesh.ny ; j++){ */
  /*   for(i = 0; i < mesh.nx; i++){ */
  /*     k = mesh.nx*j + i; */
  /*     fprintf(fid,"%e %e %e\n",mesh.dx*Real(i),mesh.dy*Real(j),Real(0.0)); */
  /*   } */
  /* } */
  /* fprintf(fid,"CELLS %d %d\n",Index(2)*mesh.ncell(),Index(4)*Index(2)*mesh.ncell()); */
  /* for(j = 0; j < mesh.ny ; j++){ */
  /*   for(i = 0; i < mesh.nx; i++){ */
  /*     k = mesh.nx*j + i; */
  /*     fprintf(fid,"%d %d %d %d\n",Index(3),k,k+Index(1),k + mesh.nx + Index(1)); */
  /*     fprintf(fid,"%d %d %d %d\n",Index(3),k,k + mesh.nx + Index(1),k+mesh.nx); */
  /*   } */
  /* } */
  /* fprintf(fid,"CELL_TYPES %d\n",Index(2)*mesh.ncell()); */
  /* for(i = 0; i < Index(2)*mesh.ncell_x*ncell_y; i++){ */
  /*   fprintf(fid,"%d\n",5); */
  /* } */
  /* fprintf(fid,"POINT_DATA %d\n",mesh.npoin()); */
  /* fprintf(fid,"SCALARS pressure DOUBLE\n"); */
  /* fprintf(fid,"LOOKUP_TABLE default\n"); */
  /* for(j = 0; j < mesh.ny ; j++){ */
  /*   for(i = 0; i < mesh.nx; i++){ */
  /* 	k = mesh.nx*j + i; */
  /* 	state = state_array[k]; */
  /* 	fprintf(fid,"%e\n",density); */
  /*   } */
  /* } */

  /* /\* fprintf(fid,"\n"); *\/ */
  /* /\* fprintf(fid,"\n"); *\/ */

  /* fclose(fid); */

  output << "# vtk DataFile Version 4.2" <<"\n";
  output << "Incompressible Euler" <<"\n";
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
  output << "SCALARS pressure DOUBLE" <<"\n";
  output << "LOOKUP_TABLE default" <<"\n";
  for(j = 0; j < mesh.ny ; j++){
    for(i = 0; i < mesh.nx; i++){

  	k = mesh.nx*j + i;

  	state = state_array[k];
  	output << density << "\n";

	/* state = prim2cons_func(gamma,State(state_array[k])); */
  	/* output << pressure << "\n"; */

    }
  }

  output << "VECTORS velocity DOUBLE" <<"\n";
  for(j = 0; j < mesh.ny ; j++){
    for(i = 0; i < mesh.nx; i++){

  	k = mesh.nx*j + i;

  	/* state = state_array[k]; */
	state = prim2cons_func(gamma,State(state_array[k]));

  	output << get_x(velocity) << " " << get_y(velocity)
  	       << " " << get_z(momentum)/density << "\n";

    }
  }

  output << " " <<"\n";
  output.close();

};
