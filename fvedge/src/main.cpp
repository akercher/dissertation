/*************************************************************************/
/* File   : main.c                                                       */
/* Author : A. Kercher                                                   */
/*-----------------------------------------------------------------------*/
/* References:                                                           */
/*   [1] J. Stone, & T. Gardiner, "A simple unsplit Godunov              */
/*       method for multidimensional MHD", New Astronomy 14,             */
/*       (2009), 139-148.                                                */
/*   [2] J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon         */
/*      "Athena: A new code for astrophysical MHD", ApJS, (2008)         */
/*   [3] Thrust - code at the speed of light.                            */
/*                url: https://github.com/thrust/thrust                  */
/*-----------------------------------------------------------------------*/
/* Copyright:                                                            */
/*                                                                       */
/*   This file is part of fvsmp.                                         */
/*                                                                       */
/*     fvedge is free software: you can redistribute it and/or modify    */
/*     it under the terms of the GNU General Public License version 3    */
/*     as published by the Free Software Foundation.                     */
/*                                                                       */
/*     fvedge is distributed in the hope that it will be useful,         */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     */
/*     GNU General Public License for more details.                      */
/*                                                                       */
/*     You should have received a copy of the GNU General Public License */
/*     along with fvsmp.  If not, see <http://www.gnu.org/licenses/>.    */
/*                                                                       */
/*************************************************************************/

#include "constrained_transport.h"
// #include "edge_solver.h"
#include "mesh_generate.h"
#include "coloring.h"
#include "data_io.h"

int main(int argc, char* argv[]){


  std::ofstream outputr;
  Real gamma;// = Real(5.0)/Real(3.0);
  Real da,da_inv;
  Real dt_fct;
  Index edge_bound_index; // index of first boundary edge
  Index fct; // = 0;
  Index mhd_source; // = 0;
  Real tf; // = 0.2;
  Real dt;
  Real cwm;
  Real cdiss;
  Index max_steps; 
  Index rk_stages; 
  Mesh mesh;
  Field field;
  Real disc_x;
  State state_l, state_r;
  RealArray dual_vol;
  RealArray wave_speed;
  Vector4Array lsq_inv; // least squares inverse
  StateArray state; // conservative state variables
  StateArray state_n; // conservative state variables at time = n
  StateGradiantArray lsq_grad; // gradiant for least squares approximation
  InterpStateArray interp_states; // interpolated primitive state variables at interface
  StateArray residual; //residuals
  StateArray flux; //flux at interface
  EdgeArray edge;
  BoundaryNodeArray bnode;
  BoundaryFaceArray bface;
  InterfaceArray face;
  StateArray antidiffusion; //antidiffusion

  Index nsteps_out = 100;
  Index output_count;
  char base_name[240];
  char file_name[240];
  std::ofstream output;
  std::ifstream input;
  std::ifstream mesh_data;

  input.open(argv[1]); // <-- opens input file
  
  max_steps = Index(1.0e4);

  /* initialize field */
  field.init();

  Index prob = Index(0);
  Index ct = Index(1);
  read_configuration(input, fct, ct, prob, field.Cour, cdiss, cwm, max_steps, 
		     rk_stages, mesh, disc_x, tf, gamma, state_l, state_r);

  /* initialize mesh */
  mesh.init();
  // mesh.iface_d = iface_; //1:tri, 0:quad
  if (mesh.ndim < 2) ct = Index(0);

  mesh.generate();

  thr::transform_n(make_device_counting_iterator(),mesh.ncell(),mesh.cells.begin(),connectivity(mesh.ncell_x));

  /* initialize conservative state variables */
  state.resize(mesh.npoin());

  /* initialize conservative state variables at time = n*/
  state_n.resize(mesh.npoin());

  /* Initialize dual vol */
  dual_vol.resize(mesh.npoin());

  // initialize edges
  edge.resize(mesh.nface());

  /* initialize residuals */
  residual.resize(mesh.npoin());

  /* initialize antidiffusive flux */
  antidiffusion.resize(mesh.nface());

  /* initialize least squares inverse*/
  lsq_inv.resize(mesh.npoin());

  /* initialize least squares gradiant*/
  lsq_grad.resize(mesh.npoin());

  // /* initialize difference of primitive state variables */
  // state_diff.resize(mesh.ncell());

  /* initialize difference of primitive state variables */
  interp_states.resize(mesh.nface());

  // initialize boundary nodess
  bnode.resize(mesh.nbnode());
  thr::fill_n(bnode.begin(),mesh.nbnode(),BoundaryNode(Coordinate(Real(0.0),Real(0.0)),
						       Index(-100),
						       Index(-100)));

  // initialize boundary faces
  bface.resize(mesh.nboun_x() + mesh.nboun_y());

  // initialize wave speeds
  wave_speed.resize(mesh.npoin());

  // define iterators
  StateIterator state_iter(state.begin());
  StateIterator antidiffusion_iter(antidiffusion.begin());
  EdgeIterator edge_iter(edge.begin());
  BoundaryNodeIterator bnode_iter(bnode.begin());
  BoundaryFaceIterator bface_iter(bface.begin());
  InterpStateIterator interp_states_iter(interp_states.begin());
  StateIterator residual_iter(residual.begin());
  Vector4Iterator lsq_inv_iter(lsq_inv.begin());
  StateGradiantIterator lsq_grad_iter(lsq_grad.begin());
  RealIterator wave_speed_iter(wave_speed.begin());

  RealArray emf_z; // electromotive force for CT 
  RealArray emf_z_poin; // electromotive force at points
  RealArray bn_edge; // normal B at interface
  RealArray bn_edge_n; // normal B at interface at time n

  // if (ct > Index(0)){

    emf_z.resize(mesh.ncell());
    thr::fill_n(emf_z.begin(),emf_z.size(),Real(0.0));
    
    emf_z_poin.resize(mesh.npoin());
    thr::fill_n(emf_z_poin.begin(),emf_z_poin.size(),Real(0.0));
    
    bn_edge.resize(mesh.nface());
    thr::fill_n(bn_edge.begin(),bn_edge.size(),Real(0.0));
    
    bn_edge_n.resize(mesh.nface());
    thr::fill_n(bn_edge_n.begin(),bn_edge_n.size(),Real(0.0));
    
  // }

  RealIterator emf_z_iter(emf_z.begin());
  RealIterator emf_z_poin_iter(emf_z_poin.begin());
  RealIterator bn_edge_iter(bn_edge.begin());    

  // create offset for coloring algorithm
  Offset offset;  
  offset.iface_d = mesh.iface_d; //1:tri, 0:quad
  offset.ncolors_per_dim = Index(2)*mesh.ndim;
  offset.init(mesh.ndim, mesh.ncell_x, mesh.ncell_y, mesh.nboun_x(), mesh.nboun_y(),
  	      mesh.btype_x, mesh.btype_y);

  edge_iter = edge.begin();
  edge_bound_index = Index(0);
  // Index ncolors = offset.ncolors - offset.nboun_colors;
  Index interior_ncolors = offset.ncolors - offset.ncolors_per_dim;
  for(Index color_index = 0; color_index < interior_ncolors; color_index++)
    {
      thr::transform_n(make_device_counting_iterator(),
  		       offset.faces_per_color[color_index],
  		       edge_iter,
  		       edges_init_2d(color_index,
  				     offset.ncolors_per_dim,
  				     mesh.ncell_x,
  				     mesh.ncell_y,
  				     mesh.btype_x,
  				     mesh.btype_y,
				     offset.iface_d,
  				     mesh.dx,
  				     mesh.dy));
      
      edge_iter += offset.faces_per_color[color_index];
      edge_bound_index += offset.faces_per_color[color_index];
    }

  for(Index color_index = interior_ncolors; color_index < offset.ncolors; color_index++)
    {
      // printf("faces_per_color[%d] = %d\n",color_index,Index(offset.faces_per_color[color_index]));
      thr::transform_n(make_device_counting_iterator(),
		       offset.faces_per_color[color_index],
		       thr::make_zip_iterator(thr::make_tuple(edge_iter,
							      bface_iter)),
		       edge_bounds_init_2d<thr::tuple<Edge,BoundaryFace> >(color_index,
									   offset.ncolors_per_dim,
									   (offset.iface_d*mesh.ndim),
									   mesh.ncell_x,
									   mesh.ncell_y,
									   mesh.btype_x,
									   mesh.btype_y,
									   offset.iface_d,
									   mesh.dx,
									   mesh.dy,
									   bnode_iter));
      
      edge_iter += offset.faces_per_color[color_index];
      bface_iter += offset.faces_per_color[color_index];
    }
  // reset iterator
  edge_iter = edge.begin();
  bface_iter = bface.begin();

  // for (Index i=0; i<mesh.nface();i++){
  // // for (Index i=0; i<12;i++){
  //   print_edges_host(i,edge[i]);
  // }

  /*-----------------------------------------------------------------*/
  /* Initialize node centered consevative state variables            */
  /*-----------------------------------------------------------------*/  

  if (prob == Index(0)){
    sprintf(base_name,"bin/constant_output");
    thr::fill_n(state.begin(),state.size(),State(Real(1.0),
						 Vector(Real(1.0),Real(1.0),Real(1.0)),
						 Real(1.0),
						 Vector(Real(0.0),Real(0.0),Real(0.0))));  
  }
  else if (prob == Index(1)){
    sprintf(base_name,"bin/constant_output");
    thr::fill_n(state.begin(),state.size(),State(Real(1.0),
						 Vector(Real(1.0),Real(1.0),Real(1.0)),
						 Real(1.0),
						 Vector(Real(1.0),Real(1.0),Real(1.0))));  
  }
  else if (prob == Index(2)){
    gamma = Real(5.0)/Real(3.0);
    nsteps_out = 100;
    sprintf(base_name,"bin/orzag_tang");
    mesh.btype_x = Index(1);
    mesh.btype_y = Index(1);
    thr::transform_n(make_device_counting_iterator(),
		     state.size(),
		     state.begin(),
		     orszag_tang_init(mesh.ncell_x,
				      mesh.dx,
				      mesh.dy));
    thr::transform_n(state.begin(),state.size(),state.begin(),prim2cons(gamma));
  }

  else if(prob == Index(3)){  
    gamma  = 1.4;
    nsteps_out = 500;
    sprintf(base_name,"bin/kh_instability");
    mesh.btype_x = Index(1);
    mesh.btype_y = Index(1);
    thr::transform_n(make_device_counting_iterator(),
		     state.size(),
		     state.begin(),
		     kh_instability_init(mesh.nx,
					 mesh.dx,
					 mesh.dy,
					 mesh.Lx,
					 mesh.Ly));
    thr::transform_n(state.begin(),state.size(),state.begin(),prim2cons(gamma));
  }
  
  else if(prob == Index(4)){  
    gamma = Real(5.0)/Real(3.0);
    nsteps_out = 50;
    sprintf(base_name,"bin/blast_wave");
    mesh.btype_x = Index(1);
    mesh.btype_y = Index(1);
    thr::transform_n(make_device_counting_iterator(),
		     state.size(),
		     state.begin(),
		     blast_wave_init(mesh.nx,
				     mesh.dx,
				     mesh.dy,
				     mesh.Lx,
				     mesh.Ly));
    thr::transform_n(state.begin(),state.size(),state.begin(),prim2cons(gamma));
  }
  else{
    mesh.btype_x = Index(0);
    mesh.btype_y = Index(1);
    sprintf(base_name,"bin/shock_tube");
    nsteps_out = 100;

    thr::transform_n(make_transform_iterator(
					     make_transform_iterator(make_device_counting_iterator(),
								     cells_init(mesh.ndim,
										mesh.nx,
										mesh.ny,
										mesh.dx,
										mesh.dy)),
					     shock_tube_init(disc_x,
							     state_l,
							     state_r)),
		     mesh.npoin(),
		     state.begin(),
		     prim2cons(gamma));

    // Real angle = Real(0.087266);
    // Real angle = Real(0.0);
    // thr::transform_n(state.begin(),
    // 		     state.size(),
    // 		     state.begin(),
    // 		     rotate_field(angle));
  }

  /*-----------------------------------------------------------------*/
  /* Initialize CT variables                                         */
  /*-----------------------------------------------------------------*/  
  // if(ct > Index(0)){
    /* initialize emf */

    for(Index i=0; i < offset.ncolors; i++){
      thr::transform_n(edge_iter,
		       offset.faces_per_color[i],	    
		       bn_edge_iter,
		       init_interface_bfield(state_iter));

      edge_iter += offset.faces_per_color[i];
      bn_edge_iter += offset.faces_per_color[i];
    }
    // reset iterator
    edge_iter = edge.begin();
    bn_edge_iter = bn_edge.begin();
  // }
  // for(Index i = 0; i < mesh.nface(); i++){
  //   printf("interface_bn[%d] = %f\n",i,Real(interface_bn[i]));
  // }
  
  
  // calculate dual volume at nodes
  thr::transform_n(make_device_counting_iterator(),
		   dual_vol.size(),
		   dual_vol.begin(),
		   calc_dual_vol(offset.iface_d,mesh.nx,mesh.ny,mesh.dx,mesh.dy));
  
  // calculate inverse least squares matrix
  thr::fill_n(lsq_inv.begin(),lsq_inv.size(),Vector4(Real(0.0),Real(0.0),Real(0.0),Real(0.0)));
  
  // loop over edges
  for(Index i=0; i < offset.ncolors; i++){
    thr::for_each_n(edge_iter,
  		    offset.faces_per_color[i],	    
  		    least_sq_inv_matrix(mesh.dx,
  					mesh.dy,
  					lsq_inv_iter));
    
    edge_iter += offset.faces_per_color[i];
  }
  // reset iterator
  edge_iter = edge.begin();

  // for(Index i = 0; i < mesh.npoin(); i++){
  //   printf("[%d] %f %f %f %f\n",i,
  // 	   thr::get<0>(Vector4(lsq_inv[i])),
  // 	   thr::get<1>(Vector4(lsq_inv[i])),
  // 	   thr::get<2>(Vector4(lsq_inv[i])),
  // 	   thr::get<3>(Vector4(lsq_inv[i])));
  //   // print_states_host(i,state_iter[i]);
  // }

  thr::transform_n(lsq_inv.begin(),lsq_inv.size(),lsq_inv.begin(),least_sq_inv_2x2());

  // for(Index i = 0; i < mesh.npoin(); i++){
  //   printf("[%d] %f %f %f %f\n",i,
  // 	   thr::get<0>(Vector4(lsq_inv[i])),
  // 	   thr::get<1>(Vector4(lsq_inv[i])),
  // 	   thr::get<2>(Vector4(lsq_inv[i])),
  // 	   thr::get<3>(Vector4(lsq_inv[i])));
  //   // print_states_host(i,state_iter[i]);
  // }

  output_count = Index(0);
  Index ksteps = Index(0);
  Real time = Real(0.0);

  /* timing variables */
  Timer program_timer, face_timer, cell_timer;
  Real faces_per_cpu_sec, cells_per_cpu_sec;
  Real faces_per_wall_sec, cells_per_wall_sec;
  Real face_cycles_per_cpu_sec, cell_cycles_per_cpu_sec;
  Real face_cycles_per_wall_sec, cell_cycles_per_wall_sec;

  cells_per_wall_sec = Real(0.0);
  faces_per_wall_sec = Real(0.0);
  cell_cycles_per_cpu_sec = Real(0.0);
  face_cycles_per_cpu_sec = Real(0.0);
  cell_cycles_per_wall_sec = Real(0.0);
  face_cycles_per_wall_sec = Real(0.0);

  program_timer.start();

  /*-----------------------------------------------------------------*/
  /* Main loop                                                       */
  /*-----------------------------------------------------------------*/

  while((time < tf) & (ksteps < max_steps)){
    
    // store state at begining of time step for two-step Runge-Kutta 
    thr::copy(state.begin(),state.end(),state_n.begin());

    // store edge states at begining of time-step
    thr::copy(bn_edge.begin(),bn_edge.end(),bn_edge_n.begin());

    
    thr::fill_n(wave_speed.begin(),wave_speed.size(),Real(0.0));

    /*-----------------------------------------------------------------*/
    /* Compute time-step                                               */
    /*-----------------------------------------------------------------*/  
    
    ksteps += 1;
    
    Real rk_coeff;
    
    face_timer.start();
    cell_timer.start();
    
    for(Index istage = 0; istage < rk_stages; istage++){ ////////////////////////////////////////////////////
      
      /*-----------------------------------------------------------------*/
      /* Compute Residuals at step n                                     */
      /*-----------------------------------------------------------------*/
      
      rk_coeff = Real(1.0)/(Real(istage) + Real(1.0));
      set_state_const(Real(0.0),residual);
      
      // set gradiants to zero
      thr::fill_n(lsq_grad.begin(),lsq_grad.size(),
		  StateGradiant(State(Real(0.0),
				      Vector(Real(0.0),Real(0.0),Real(0.0)),
				      Real(0.0),
				      Vector(Real(0.0),Real(0.0),Real(0.0))),
				State(Real(0.0),
				      Vector(Real(0.0),Real(0.0),Real(0.0)),
				      Real(0.0),
				      Vector(Real(0.0),Real(0.0),Real(0.0)))));
      
      // compute gradiants
      for(Index i=0; i < offset.ncolors; i++){
	thr::for_each_n(edge_iter,
			offset.faces_per_color[i],	    
			least_sq_gradiants(mesh.dx,
					   mesh.dy,
					   gamma,
					   state_iter,
					   lsq_grad_iter));
	edge_iter += offset.faces_per_color[i];
      }
      // reset iterator
      edge_iter = edge.begin();

      //  for(Index i = 0; i < mesh.npoin(); i++){
      // 	 printf("[%d] %f %f\n",i,get_x(thr::get<1>(get_x(StateGradiant(lsq_grad_iter[i])))),
      // 	 	get_y(thr::get<1>(get_y(StateGradiant(lsq_grad_iter[i])))));
      // }
      
      // multiple by precomputed inverse
      thr::transform_n(lsq_inv.begin(),lsq_inv.size(),lsq_grad.begin(),lsq_grad.begin(),
		       least_sq_grad_inv());

      //  for(Index i = 0; i < mesh.npoin(); i++){
      // 	 printf("[%d] %f %f\n",i,get_x(thr::get<1>(get_x(StateGradiant(lsq_grad_iter[i])))),
      // 	 	get_y(thr::get<1>(get_y(StateGradiant(lsq_grad_iter[i])))));
      // }


  if(1 == 1){///////////////////////////////////////////////////////////////////////////////////////////// 

      // interpolate to interface
      interp_states_iter = interp_states.begin();
      for(Index i=0; i < offset.ncolors; i++){
	thr::transform_n(edge_iter,
			 offset.faces_per_color[i],
			 interp_states_iter,
			 gradiant_reconstruction(gamma,
						 state_iter,
						 lsq_grad_iter));
	offset.update(i);
	edge_iter += offset.faces_per_color[i];
	interp_states_iter += offset.faces_per_color[i];
      }
      // reset iterators
      edge_iter = edge.begin();
      interp_states_iter = interp_states.begin();

      //  for(Index i = 0; i < mesh.nface(); i++){
      // 	 printf("[%d][%d] %f %f\n",get_x(thr::get<2>(Edge(edge[i]))),get_y(thr::get<2>(Edge(edge[i]))),
      // 		get_x(thr::get<1>(get_x(InterpState(interp_states_iter[i])))),
      // 	 	get_y(thr::get<1>(get_y(InterpState(interp_states_iter[i])))));
      // }
      
      // build residual
      for(Index i=0; i < offset.ncolors; i++){
	thr::transform_n(edge_iter,
			 offset.faces_per_color[i],
			 interp_states_iter,
			 antidiffusion_iter,
			 residual_op(gamma,
				     wave_speed_iter,
				     emf_z_iter,
				     residual_iter));
	
	edge_iter += offset.faces_per_color[i];
	interp_states_iter += offset.faces_per_color[i];
	antidiffusion_iter += offset.faces_per_color[i];
      }
      // reset iterators
      edge_iter = edge.begin();
      interp_states_iter = interp_states.begin();
      antidiffusion_iter = antidiffusion.begin();      

      // apply boundary conditions
      if (mesh.btype_x == Index(1)){
	// periodic left/right
	thr::for_each_n(make_device_counting_iterator(),
			mesh.ny,
			periodic_bcs(mesh.nx,
				     mesh.ny,
				     mesh.nx,
				     wave_speed_iter,
				     emf_z_iter,				   
				     residual_iter));
      }
      else{
	// outflow left/right
	for(Index i = interior_ncolors; i < (offset.ncolors - Index(2)); i++){
	  thr::for_each_n(bface_iter,
	  		  offset.faces_per_color[i],
	  		  outflow_bcs(mesh.iface_d,
	  			      gamma,
	  			      wave_speed_iter,
	  			      bnode_iter,
	  			      emf_z_iter,
	  			      state_iter,
	  			      residual_iter));
	  
	  bface_iter += offset.faces_per_color[i];
	}
	// reset iterators
	bface_iter = bface.begin();

      }

      // apply boundary conditions
      if (mesh.btype_y == Index(1)){
	// periodic top/bottom
	thr::for_each_n(make_device_counting_iterator(),
			mesh.nx,
			periodic_bcs(mesh.nx,
				     mesh.ny,
				     Index(1),
				     wave_speed_iter,
				     emf_z_iter,
				     residual_iter));
	
      }
      else{
	bface_iter += Index(2);
	// outflow top/bottom
	for(Index i = (interior_ncolors + Index(2)); i < (offset.ncolors); i++){
	  thr::for_each_n(bface_iter,
	  		  offset.faces_per_color[i],
	  		  outflow_bcs(mesh.iface_d,
	  			      gamma,
	  			      wave_speed_iter,
	  			      bnode_iter,
	  			      emf_z_iter,
	  			      state_iter,
	  			      residual_iter));
	  
	  bface_iter += offset.faces_per_color[i];
	}
	// reset iterators
	bface_iter = bface.begin();
      }

      // printf("\n");
      // for(Index i = 0; i < mesh.npoin(); i++){
      // 	print_states_host(i,State(residual_iter[i]));
      // }


      /*-----------------------------------------------------------------*/
      /* compute emfs                                                    */
      /*-----------------------------------------------------------------*/  	  	        
      for(Index i=0; i < offset.ncolors; i++){
	thr::for_each_n(thr::make_zip_iterator(thr::make_tuple(edge_iter,
							       antidiffusion_iter)),
			offset.faces_per_color[i],
			emf_upwind_calc<thr::tuple<Edge,State> >(mesh.nx,
								 mesh.ny,
								 emf_z_iter,
								 state_iter));
	
	edge_iter += offset.faces_per_color[i];
	antidiffusion_iter += offset.faces_per_color[i];
      }
      // reset iterators
      edge_iter = edge.begin();
      antidiffusion_iter = antidiffusion.begin();      

      // for(Index i = 0; i < mesh.ncell(); i++){
      // 	printf("[%d] %f\n",i,Real(emf_z_iter[i]));
      // }      


      /*-----------------------------------------------------------------*/
      /* Update solution: Two-stage Runge-Kutta                          */
      /*    1. state = state_n - (dt/vol)*residual(state_n)              */
      /*    2. state^{n+1} = 0.5*(state^n + state)                       */
      /*                     - 0.5*(dt/vol)*residual(state)              */
      /*-----------------------------------------------------------------*/

      /*-----------------------------------------------------------------*/
      /* compute time step if stage one                                  */
      /*-----------------------------------------------------------------*/  	  	        
      if(istage < 1){

	Real init(MAXIMUM_NUM); // initial value for reduce operation

	dt  = thr::reduce(make_transform_iterator(thr::make_zip_iterator(thr::make_tuple(dual_vol.begin(),
											 wave_speed.begin())),
						  time_step<thr::tuple<Real,Real> >()),
			  make_transform_iterator(thr::make_zip_iterator(thr::make_tuple(dual_vol.end(),
											 wave_speed.end())),
						  time_step<thr::tuple<Real,Real> >()),
			  init,array_min());
	
	if((time + dt) > tf) dt = tf - time;
	time += dt;
	
	if((ksteps-Index(1)) % 100 == 0){

	  std::cout << "ksteps = " << ksteps << " , " 
		    << "time = " << time << " , "
	    // << "max_wave_speed = " << field.max_wave_speed << " , " 
		    << "dt = " << dt << " , " 
		    << "cells_per_sec = " << cells_per_wall_sec << " , "  
		    << "face_per_sec = " << faces_per_wall_sec // << " , "  
		    << std::endl;
	}

      }
      /*-----------------------------------------------------------------*/
      /* Arverage states if stage two                                    */
      /*-----------------------------------------------------------------*/  	  	  
      else{
	thr::transform(state_n.begin(),
		       state_n.end(),
		       state.begin(),
		       state.begin(),
		       sum_and_scale_states(half));
      }
      /*-----------------------------------------------------------------*/
      /* Time integration                                                */
      /*-----------------------------------------------------------------*/  
      
      thr::transform_n(thr::make_zip_iterator(thr::make_tuple(state.begin(),residual.begin())),
		       state.size(),
		       dual_vol.begin(),
		       state.begin(),
		       integrate_time<thr::tuple<State,State> >(field.Cour*rk_coeff*dt));
      
      for(Index i=0; i < offset.ncolors; i++)
	{
	  thr::transform_n(edge_iter,
			   offset.faces_per_color[i],	    
			   bn_edge_iter,
			   bn_edge_iter,
			   integrate_ct(mesh.nx,
					mesh.ny,
					Index(1),
					field.Cour*rk_coeff*dt,
					emf_z_iter,
					state_iter));
	  edge_iter += offset.faces_per_color[i];
	  bn_edge_iter += offset.faces_per_color[i];
	}
      // reset iterators
      edge_iter = edge.begin();
      bn_edge_iter = bn_edge.begin();
      
      // for(Index i = 0; i < mesh.nface(); i++){
      // 	printf("[%d] %f\n",i,Real(bn_edge_iter[i]));
      // }
      // for(Index i = 0; i < mesh.nface(); i++){
      // 	print_states_host(i,State(residual_iter[i]));
      // 	// print_states_host(i,get_x(InterpState(interp_states_iter[i])));
      // }
    }////////////////////////////////////////////////////////////////////////////////////
      
    } ////////////// end loop over rk-stages
    face_timer.stop();
    cell_timer.stop();
    
    faces_per_wall_sec = mesh.nface()/face_timer.elapsed_wall_time();
    cells_per_wall_sec = mesh.ncell()/cell_timer.elapsed_wall_time();

    if((ksteps-Index(1)) % nsteps_out == 0){

      sprintf(file_name,"%s_%05d.vtk",base_name,output_count);
      output.open(file_name);
      output_vtk_legacy(output, mesh, gamma, state);
      output_count += Index(1);
    }
    
  }  // end time integration
  
  program_timer.stop();

  cell_cycles_per_cpu_sec = Real(ksteps)*Real(mesh.ncell())/program_timer.elapsed_cpu_time();
  cell_cycles_per_wall_sec = Real(ksteps)*Real(mesh.ncell())/program_timer.elapsed_wall_time();

  face_cycles_per_cpu_sec = Real(ksteps)*Real(mesh.nface())/program_timer.elapsed_cpu_time();
  face_cycles_per_wall_sec = Real(ksteps)*Real(mesh.nface())/program_timer.elapsed_wall_time();

  // std::cout << mesh.dx << " , " << mesh.dy*mesh.dx << std::endl;

  // -----------------------------------------------------------------
  std::cout << " " << std::endl;
  Index istart = Index(0.0);//(mesh.ndim - Index(1))*mesh.ncell_x;
  Index iend = mesh.ndim*mesh.ncell_x;
  // Index iend = mesh.ncell_x*mesh.ncell_y;

  // for(Index i = istart; i < iend; i++){
  //   print_states_host(i,State(state_iter[i]));
  // }

  // for(Index i = 0; i < mesh.nface(); i++){
  //   printf("interface_bn[%d] = %f\n",i,Real(interface_bn[i]));
  // }


  std::cout << "ksteps = " << ksteps << " , " << "time = " << time << " , "
	    // << "max_wave_speed = " << field.max_wave_speed << " , " 
	    // << "dt = " << dt << " , " 
	    // << "cells_per_sec = " << cells_per_sec << " , " 
	    // << "face_per_sec = " << faces_per_sec //<< " , " 
	    << std::endl;

  std::cout << " " << std::endl;
  std::cout << "CPU time elapsed : " << program_timer.elapsed_cpu_time() << "\n"
	    << "cell-cycles/cpu-sec : " << cell_cycles_per_cpu_sec << "\n"
	    << "face-cycles/cpu-sec : " << face_cycles_per_cpu_sec
	    << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Wall time elapsed : " << program_timer.elapsed_wall_time() << "\n"
	    << "cell-cycles/wall-sec : " << cell_cycles_per_wall_sec << "\n"
	    << "face-cycles/wall-sec : " << face_cycles_per_wall_sec //<< " , "
	    // << "Avg. cells per sec.: " << avg_cells_per_sec/Real(ksteps) << " , "
	    // << "Avg. faces per sec.: " << avg_faces_per_sec/Real(ksteps) //<< " , "
	    << std::endl;
  std::cout << " " << std::endl;

  // for(Index i = 0; i < mesh.npoin(); i++){
  //   print_states_host(i,State(state_iter[i]));
  //   // print_states_host(i,get_x(InterpState(interp_states_iter[i])));
  // }

  // output.open("bin/output.dat");
  // output_states(output, mesh.ncell_x, mesh.ncell_y, state);

  // output.open("dat/kh_instability.vtk");
  sprintf(file_name,"%s_%05d.vtk",base_name,output_count);
  output.open(file_name);
  output_vtk_legacy(output, mesh, gamma, state);

}
