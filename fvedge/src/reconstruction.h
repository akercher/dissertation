/*******************************************************************/
/* File   : reconstruction.h                                       */
/* Author : A. Kercher                                             */
/*-----------------------------------------------------------------*/
/*******************************************************************/

/* #include "thrust_wrapper.h" */
/* #include "defs.h" */
/* #include <iostream> */
/* #include <fstream> */
/* #include <sstream> */
/* #include <string> */
/* #include <cmath> */
/* #include <float.h> */
/* #include <sys/time.h> */
/* #include <ctime> */

/*-------------------------------------------------*/
/* Prototypes                                      */
/*-------------------------------------------------*/
__host__ __device__ State cons2prim_func (Real gamma, State state);
__host__ __device__ State prim2cons_func (Real gamma, State state);

/*****************************************************/
/*                                                   */  
/* Compute gradiants with least squares              */  
/*                                                   */  
/*                                                   */  
/*---------------------------------------------------*/
/*  Input :                                          */
/*  Output :                                         */
/*****************************************************/
struct gradiant_reconstruction : public thr::unary_function<Edge,InterpState>
{
  Real _gamma;
  StateIterator _state_iter;
  StateGradiantIterator _lsq_grad_iter;

 gradiant_reconstruction(Real gamma,
			 StateIterator state_iter,
			 StateGradiantIterator lsq_grad_iter)
   : _gamma(gamma) 
    ,_state_iter(state_iter)
    ,_lsq_grad_iter(lsq_grad_iter) {}

  __host__ __device__
    InterpState operator()(const Edge& edge) const
  {

    Index i;
    Real dx,dy;
    Real ddi,dvxi,dvyi,dvzi,dpgi,dbxi,dbyi,dbzi;
    Real ddj,dvxj,dvyj,dvzj,dpgj,dbxj,dbyj,dbzj;
    Real sgn;
    State grad_x,grad_y;
    State interp_state_i, interp_state_j;
    Real df[8], dp[8], dm[8],dwi[8],dwj[8],da[8];

    Coordinate edge_vec = thr::get<1>(Edge(edge));
    Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec) + get_y(edge_vec)*get_y(edge_vec));
    Real edge_vec_mag_inv = Real(1.0)/edge_vec_mag;
    Coordinate normal = Coordinate(get_x(edge_vec)*edge_vec_mag_inv,get_y(edge_vec)*edge_vec_mag_inv);

    Index index_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index index_j = thr::get<1>(thr::get<2>(Edge(edge)));

    State cons_state_i = State(this->_state_iter[index_i]);
    State cons_state_j = State(this->_state_iter[index_j]);

    // convert primitive state variables
    State state_i;
    State state_j;

    state_i = cons2prim_func(this->_gamma,cons_state_i);
    state_j = cons2prim_func(this->_gamma,cons_state_j);

    // gradiants at state i
    grad_x = get_x(StateGradiant(this->_lsq_grad_iter[index_i]));
    grad_y = get_y(StateGradiant(this->_lsq_grad_iter[index_i]));

    /* printf("[%d][%d] %f %f\n",index_i,index_j,get_x(thr::get<1>(grad_x)),get_x(thr::get<1>(grad_y))); */


    ddi = half*edge_vec_mag*(thr::get<0>(grad_x)*get_x(normal) + thr::get<0>(grad_y)*get_y(normal));
    dvxi = half*edge_vec_mag*(get_x(thr::get<1>(grad_x))*get_x(normal) 
			      + get_x(thr::get<1>(grad_y))*get_y(normal));
    dvyi = half*edge_vec_mag*(get_y(thr::get<1>(grad_x))*get_x(normal) 
			      + get_y(thr::get<1>(grad_y))*get_y(normal));
    dvzi = half*edge_vec_mag*(get_z(thr::get<1>(grad_x))*get_x(normal) 
			      + get_z(thr::get<1>(grad_y))*get_y(normal));
    dpgi = half*edge_vec_mag*(thr::get<2>(grad_x)*get_x(normal) 
			      + thr::get<2>(grad_y)*get_y(normal));
    dbxi = half*edge_vec_mag*(get_x(thr::get<3>(grad_x))*get_x(normal) 
			      + get_x(thr::get<3>(grad_y))*get_y(normal));
    dbxi = half*edge_vec_mag*(get_y(thr::get<3>(grad_x))*get_x(normal) 
			      + get_y(thr::get<3>(grad_y))*get_y(normal));
    dbxi = half*edge_vec_mag*(get_z(thr::get<3>(grad_x))*get_x(normal) 
			      + get_z(thr::get<3>(grad_y))*get_y(normal));

    // gradiants at state j
    grad_x = get_x(StateGradiant(this->_lsq_grad_iter[index_j]));
    grad_y = get_y(StateGradiant(this->_lsq_grad_iter[index_j]));

    ddj = half*edge_vec_mag*(thr::get<0>(grad_x)*get_x(normal) 
			     + thr::get<0>(grad_y)*get_y(normal));
    dvxj = half*edge_vec_mag*(get_x(thr::get<1>(grad_x))*get_x(normal) 
			      + get_x(thr::get<1>(grad_y))*get_y(normal));
    dvyj = half*edge_vec_mag*(get_y(thr::get<1>(grad_x))*get_x(normal) 
			      + get_y(thr::get<1>(grad_y))*get_y(normal));
    dvzj = half*edge_vec_mag*(get_z(thr::get<1>(grad_x))*get_x(normal) 
			      + get_z(thr::get<1>(grad_y))*get_y(normal));
    dpgj = half*edge_vec_mag*(thr::get<2>(grad_x)*get_x(normal) 
			      + thr::get<2>(grad_y)*get_y(normal));
    dbxj = half*edge_vec_mag*(get_x(thr::get<3>(grad_x))*get_x(normal) 
			      + get_x(thr::get<3>(grad_y))*get_y(normal));
    dbxj = half*edge_vec_mag*(get_y(thr::get<3>(grad_x))*get_x(normal) 
			      + get_y(thr::get<3>(grad_y))*get_y(normal));
    dbxj = half*edge_vec_mag*(get_z(thr::get<3>(grad_x))*get_x(normal) 
			      + get_z(thr::get<3>(grad_y))*get_y(normal));

    // simple interpolation
    thr::get<0>(interp_state_i) = density_i + ddi;
    get_x(thr::get<1>(interp_state_i)) = get_x(velocity_i) + dvxi;
    get_y(thr::get<1>(interp_state_i)) = get_y(velocity_i) + dvyi;
    get_z(thr::get<1>(interp_state_i)) = get_z(velocity_i) + dvzi;
    thr::get<2>(interp_state_i) = pressure_i + dpgi;
    get_x(thr::get<3>(interp_state_i)) = get_x(bfield_i) + dbxi;
    get_y(thr::get<3>(interp_state_i)) = get_y(bfield_i) + dbyi;
    get_z(thr::get<3>(interp_state_i)) = get_z(bfield_i) + dbzi;

    thr::get<0>(interp_state_j) = density_j - ddj;
    get_x(thr::get<1>(interp_state_j)) = get_x(velocity_j) - dvxj;
    get_y(thr::get<1>(interp_state_j)) = get_y(velocity_j) - dvyj;
    get_z(thr::get<1>(interp_state_j)) = get_z(velocity_j) - dvzj;
    thr::get<2>(interp_state_j) = pressure_j - dpgj;
    get_x(thr::get<3>(interp_state_j)) = get_x(bfield_j) - dbxj;
    get_y(thr::get<3>(interp_state_j)) = get_y(bfield_j) - dbyj;
    get_z(thr::get<3>(interp_state_j)) = get_z(bfield_j) - dbzj;

    // Van Albada limiter
    /* face derivatives */
    df[0] = half*(density_j - density_i);
    df[1] = half*(get_x(velocity_j) - get_x(velocity_i));
    df[2] = half*(get_y(velocity_j) - get_y(velocity_i));
    df[3] = half*(get_z(velocity_j) - get_z(velocity_i));
    df[4] = half*(pressure_j - pressure_i);
    df[5] = half*(get_x(bfield_j) - get_x(bfield_i));
    df[6] = half*(get_y(bfield_j) - get_y(bfield_i));
    df[7] = half*(get_z(bfield_j) - get_z(bfield_i));

    dwi[0] = ddi;
    dwi[1] = dvxi;
    dwi[2] = dvyi;
    dwi[3] = dvzi;
    dwi[4] = dpgi;
    dwi[5] = dbxi;
    dwi[6] = dbyi;
    dwi[7] = dbzi;

    dwj[0] = ddj;
    dwj[1] = dvxj;
    dwj[2] = dvyj;
    dwj[3] = dvzj;
    dwj[4] = dpgj;
    dwj[5] = dbxj;
    dwj[6] = dbyj;
    dwj[7] = dbzj;

    // slope at i
    Real eps2 = (Real(0.3)*edge_vec_mag)*(Real(0.3)*edge_vec_mag)*(Real(0.3)*edge_vec_mag);
    for(i=0;i<Index(8);i++){
      dm[i] = Real(2.0)*dwi[i] - df[i];
      dp[i] = df[i];
      sgn = Real(1.0);
      if (std::fabs(dm[i]*dp[i]) > Real(0.0)){
	sgn = (dm[i]*dp[i])/(std::fabs(dm[i]*dp[i]));
      }
      da[i] = half*(sgn + Real(1.0))*(((dp[i]*dp[i] + eps2)*dm[i] + (dm[i]*dm[i] + eps2)*dp[i])
    				      /((dm[i]*dm[i] + dp[i]*dp[i] + Real(2.0)*eps2)));
    }

    /* printf("[%d][%d] dwi[1] = % f da[1] = %f dm[1] = %f dp[1] = %f df[1] = %f eps2 = %f\n",index_i,index_j,dwi[1],da[1],dm[1],dp[1],df[1],eps2); */

    // if negative pressure use constant interp.
    if((pressure_i + da[4]) <= Real(0.0)){
      for(i=0;i<Index(8);i++){
	da[i] = Real(0.0);
      }
    }

    // value at state i
    thr::get<0>(interp_state_i) = density_i + da[0];
    get_x(thr::get<1>(interp_state_i)) = get_x(velocity_i) + da[1];
    get_y(thr::get<1>(interp_state_i)) = get_y(velocity_i) + da[2];
    get_z(thr::get<1>(interp_state_i)) = get_z(velocity_i) + da[3];
    thr::get<2>(interp_state_i) = pressure_i + da[4];
    get_x(thr::get<3>(interp_state_i)) = get_x(bfield_i) + da[5];
    get_y(thr::get<3>(interp_state_i)) = get_y(bfield_i) + da[6];
    get_z(thr::get<3>(interp_state_i)) = get_z(bfield_i) + da[7];
    

    // slope at j
    for(i=0;i<Index(8);i++){
      dm[i] = -(Real(2.0)*dwj[i] - df[i]);
      dp[i] = -df[i];
      sgn = (dm[i]*dp[i] + Real(1.0))/(std::abs(dm[i]*dp[i]) + Real(1.0));
      da[i] = half*(sgn + Real(1.0))*(((dp[i]*dp[i] + eps2)*dm[i] + (dm[i]*dm[i] + eps2)*dp[i])
    				      /((dm[i]*dm[i] + dp[i]*dp[i] + Real(2.0)*eps2)));
    }

    // if negative pressure use constant interp.
    if((pressure_j + da[4]) <= Real(0.0)){
      for(i=0;i<Index(8);i++){
	da[i] = Real(0.0);
      }
    }
    // value at state j
    thr::get<0>(interp_state_j) = density_j + da[0];
    get_x(thr::get<1>(interp_state_j)) = get_x(velocity_j) + da[1];
    get_y(thr::get<1>(interp_state_j)) = get_y(velocity_j) + da[2];
    get_z(thr::get<1>(interp_state_j)) = get_z(velocity_j) + da[3];
    thr::get<2>(interp_state_j) = pressure_j + da[4];
    get_x(thr::get<3>(interp_state_j)) = get_x(bfield_j) + da[5];
    get_y(thr::get<3>(interp_state_j)) = get_y(bfield_j) + da[6];
    get_z(thr::get<3>(interp_state_j)) = get_z(bfield_j) + da[7];

    

    /* printf("[%d][%d] intp_di = %f intp_dj = %f\n",index_i,index_j, */
    /* 	   thr::get<0>(interp_state_i),thr::get<0>(interp_state_j)); */

    return InterpState(State(interp_state_i),State(interp_state_j));

  }
};

/*****************************************************/
/*                                                   */  
/* Compute gradiants with least squares              */  
/*                                                   */  
/*                                                   */  
/*---------------------------------------------------*/
/*  Input :                                          */
/*  Output :                                         */
/*****************************************************/
struct least_sq_gradiants : public thr::unary_function<Edge,void>
{
  Real _dx;
  Real _dy;
  Real _gamma;
  StateIterator _state_iter;
  StateGradiantIterator _lsq_grad_iter;

 least_sq_gradiants(Real dx,
		    Real dy,
		    Real gamma,
		    StateIterator state_iter,
		    StateGradiantIterator lsq_grad_iter)
   : _dx(dx) 
    ,_dy(dy)
    ,_gamma(gamma)
    ,_state_iter(state_iter)
    ,_lsq_grad_iter(lsq_grad_iter) {}

  __host__ __device__
    void operator()(const Edge& edge) const
  {
    
    Real dxi,dyi,dxj,dyj;
    Real dd,dvx,dvy,dvz,dpg,dbx,dby,dbz;
    State grad_x,grad_y;

    Coordinate edge_vec = thr::get<1>(Edge(edge));
    Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec) + get_y(edge_vec)*get_y(edge_vec));

    Index index_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index index_j = thr::get<1>(thr::get<2>(Edge(edge)));

    State cons_state_i = State(this->_state_iter[index_i]);
    State cons_state_j = State(this->_state_iter[index_j]);

    // convert primitive state variables
    State state_i;
    State state_j;

    state_i = cons2prim_func(this->_gamma,cons_state_i);
    state_j = cons2prim_func(this->_gamma,cons_state_j);

    dxi = std::fabs(get_x(edge_vec));
    dyi = std::fabs(get_y(edge_vec));

    if(index_i > index_j){
      dxi = -dxi;
      dyi = -dyi;
    }

    dxj = -dxi;
    dyj = -dyi;

    // calculate difference in primitive variables
    dd = density_j - density_i;
    dvx = get_x(velocity_j) - get_x(velocity_i);
    dvy = get_y(velocity_j) - get_y(velocity_i);
    dvz = get_z(velocity_j) - get_z(velocity_i);
    dpg = pressure_j - pressure_i;
    dbx = get_x(bfield_j) - get_x(bfield_i);
    dbx = get_y(bfield_j) - get_y(bfield_i);
    dbx = get_z(bfield_j) - get_z(bfield_i);

    // gradiants for cell i
    grad_x = get_x(StateGradiant(this->_lsq_grad_iter[index_i]));
    grad_y = get_y(StateGradiant(this->_lsq_grad_iter[index_i]));

    thr::get<0>(grad_x) += dd*dxi;
    get_x(thr::get<1>(grad_x)) += dvx*dxi;
    get_y(thr::get<1>(grad_x)) += dvy*dxi;
    get_z(thr::get<1>(grad_x)) += dvz*dxi;
    thr::get<2>(grad_x) += dpg*dxi;
    get_x(thr::get<3>(grad_x)) += dbx*dxi;
    get_y(thr::get<3>(grad_x)) += dby*dxi;
    get_z(thr::get<3>(grad_x)) += dbz*dxi;

    thr::get<0>(grad_y) += dd*dyi;
    get_x(thr::get<1>(grad_y)) += dvx*dyi;
    get_y(thr::get<1>(grad_y)) += dvy*dyi;
    get_z(thr::get<1>(grad_y)) += dvz*dyi;
    thr::get<2>(grad_y) += dpg*dyi;
    get_x(thr::get<3>(grad_y)) += dbx*dyi;
    get_y(thr::get<3>(grad_y)) += dby*dyi;
    get_z(thr::get<3>(grad_y)) += dbz*dyi;

    this->_lsq_grad_iter[index_i] = StateGradiant(State(grad_x),State(grad_y));

    // gradiants for cell j
    grad_x = get_x(StateGradiant(this->_lsq_grad_iter[index_j]));
    grad_y = get_y(StateGradiant(this->_lsq_grad_iter[index_j]));

    thr::get<0>(grad_x) -= dd*dxj;
    get_x(thr::get<1>(grad_x)) -= dvx*dxj;
    get_y(thr::get<1>(grad_x)) -= dvy*dxj;
    get_z(thr::get<1>(grad_x)) -= dvz*dxj;
    thr::get<2>(grad_x) -= dpg*dxj;
    get_x(thr::get<3>(grad_x)) -= dbx*dxj;
    get_y(thr::get<3>(grad_x)) -= dby*dxj;
    get_z(thr::get<3>(grad_x)) -= dbz*dxj;

    thr::get<0>(grad_y) -= dd*dyj;
    get_x(thr::get<1>(grad_y)) -= dvx*dyj;
    get_y(thr::get<1>(grad_y)) -= dvy*dyj;
    get_z(thr::get<1>(grad_y)) -= dvz*dyj;
    thr::get<2>(grad_y) -= dpg*dyj;
    get_x(thr::get<3>(grad_y)) -= dbx*dyj;
    get_y(thr::get<3>(grad_y)) -= dby*dyj;
    get_z(thr::get<3>(grad_y)) -= dbz*dyj;

    this->_lsq_grad_iter[index_j] = StateGradiant(State(grad_x),State(grad_y));

  }
};

/*****************************************************/
/*                                                   */  
/* Multiply gradiants by pre-computed inverse.       */  
/*                                                   */  
/*                                                   */  
/*---------------------------------------------------*/
/*  Input :                                          */
/*  Output :                                         */
/*****************************************************/
struct least_sq_grad_inv : public thr::binary_function<Vector4,StateGradiant,StateGradiant>
{

  __host__ __device__
    StateGradiant operator()(const Vector4& lsq_inv, const StateGradiant& lsq_grad) const
  {
    
    Real dx,dy;
    Real A[4];
    State grad_x,grad_y;

    A[0] = thr::get<0>(Vector4(lsq_inv));
    A[1] = thr::get<1>(Vector4(lsq_inv));
    A[2] = thr::get<2>(Vector4(lsq_inv));
    A[3] = thr::get<3>(Vector4(lsq_inv));

    grad_x = get_x(StateGradiant(lsq_grad));
    grad_y = get_y(StateGradiant(lsq_grad));

    dx = thr::get<0>(grad_x);
    dy = thr::get<0>(grad_y);
    thr::get<0>(grad_x) = A[0]*dx + A[1]*dy;
    thr::get<0>(grad_y) = A[2]*dx + A[3]*dy;

    dx = get_x(thr::get<1>(grad_x));
    dy = get_x(thr::get<1>(grad_y));
    get_x(thr::get<1>(grad_x)) = A[0]*dx + A[1]*dy;
    get_x(thr::get<1>(grad_y)) = A[2]*dx + A[3]*dy;

    dx = get_y(thr::get<1>(grad_x));
    dy = get_y(thr::get<1>(grad_y));
    get_y(thr::get<1>(grad_x)) = A[0]*dx + A[1]*dy;
    get_y(thr::get<1>(grad_y)) = A[2]*dx + A[3]*dy;

    dx = get_z(thr::get<1>(grad_x));
    dy = get_z(thr::get<1>(grad_y));
    get_z(thr::get<1>(grad_x)) = A[0]*dx + A[1]*dy;
    get_z(thr::get<1>(grad_y)) = A[2]*dx + A[3]*dy;

    dx = thr::get<2>(grad_x);
    dy = thr::get<2>(grad_y);
    thr::get<2>(grad_x) = A[0]*dx + A[1]*dy;
    thr::get<2>(grad_y) = A[2]*dx + A[3]*dy;

    dx = get_x(thr::get<3>(grad_x));
    dy = get_x(thr::get<3>(grad_y));
    get_x(thr::get<3>(grad_x)) = A[0]*dx + A[1]*dy;
    get_x(thr::get<3>(grad_y)) = A[2]*dx + A[3]*dy;

    dx = get_y(thr::get<3>(grad_x));
    dy = get_y(thr::get<3>(grad_y));
    get_y(thr::get<3>(grad_x)) = A[0]*dx + A[1]*dy;
    get_z(thr::get<3>(grad_y)) = A[2]*dx + A[3]*dy;

    dx = get_z(thr::get<3>(grad_x));
    dy = get_z(thr::get<3>(grad_y));
    get_z(thr::get<3>(grad_x)) = A[0]*dx + A[1]*dy;
    get_z(thr::get<3>(grad_y)) = A[2]*dx + A[3]*dy;

    return StateGradiant(State(grad_x),State(grad_y));
 
  }
};

/* ------------------------------------------------------------------------------------*/
/* *************************************************************************************/
/*  --- Inverse Matrix for 2x2 Least-Squares Gradient Reconstruction ---               */
/*                                                                                     */
/*  Construct a matrix for the linear least-squares(LSQ) gradient reconstruction.      */
/*  (unweighted LSQ; more accurate than weighted ones to my knowledge.)                */
/*                                                                                     */
/*  Note: it requires at least 2 non-colinear neighbors.                               */
/*                                                                                     */
/*  Example: Consider constructing (ux,uy) at i with the following stencil.            */
/*                                                                                     */
/*       3 o     o 2                                                                   */
/*          \   /                                                                      */
/*           \ /                                                                       */
/*          i *-----o 1                                                                */
/*           /|                                                                        */
/*          / |                                                                        */
/*         /  o 5      *: node in interest (i)                                         */
/*        o 4          o: neighbors (k = 1,2,3,4,5)                                    */
/*                                                                                     */
/*   5 equations:                                                                      */
/*     (x1-xi)*ux + (y1-yi)*uy = (u1-ui)                                               */
/*     (x2-xi)*ux + (y2-yi)*uy = (u2-ui)                                               */
/*     (x3-xi)*ux + (y3-yi)*uy = (u3-ui)                                               */
/*     (x4-xi)*ux + (y4-yi)*uy = (u4-ui)                                               */
/*     (x5-xi)*ux + (y5-yi)*uy = (u5-ui)                                               */
/*                                                                                     */
/*   This system is written in the matrix form:                                        */
/*                                                                                     */
/*         A*x = b,  x=(ux,uy), A=5x2 matrix, b=5x1 matrix                             */
/*                                                                                     */
/*   The least-squares problem is                                                      */
/*                                                                                     */
/*       A^T*A*x = A^T*b, (T denotes the transpose: A^T=2x5 matrix)                    */
/*                                                                                     */  
/*   which is                                                                          */
/*                                                                                     */
/*   [sum_k (xk-xi)^2]*ux       + [sum_k (xk-xi)*(yk-yi)]*uy = [sum_k (uk-ui)*(xk-xi)] */
/*   [sum_k (xk-xi)*(yk-yi)]*ux + [sum_k (yk-yi)]*uy         = [sum_k (uk-ui)*(yk-yi)] */
/*                                                                                     */
/*  This subroutine computes the inverse of (A^T*A) at every node (which depends       */
/*  only on the grid geometry), so that the gradient at a node can be computed         */
/*  by a matrix-vector multiplication, i.e., (A^T*A)^{-1}*(A^T*b),                     */
/*  (only A^T*b needs to be re-computed).                                              */
/*                                                                                     */
/*  ---------------------------------------------------------------------------------- */
/*   Input:  inode = node number                                                       */
/*                                                                                     */
/*  Output:  nodes(inode)%lsq_inv_2x2 = inverse matrix for LSQ reconstruction          */
/*  ---------------------------------------------------------------------------------- */
/*                                                                                     */
/*  References:                                                                        */
/*    [1] Original version for node-center (edge-based) finite volume discretization   */
/*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com).                      */
/*                                                                                     */
/* *********************************************************************************** */
struct least_sq_inv_matrix : public thr::unary_function<Edge,void>
{
  Real _dx;
  Real _dy;
  Vector4Iterator _lsq_inv_iter;

 least_sq_inv_matrix(Real dx,
		     Real dy,
		     Vector4Iterator lsq_inv_iter)
   : _dx(dx) 
    ,_dy(dy)
    ,_lsq_inv_iter(lsq_inv_iter) {}

  __host__ __device__
    void operator()(const Edge& edge) const
  {
    
    Real dx,dy;
    Real A[4];

    Index index_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index index_j = thr::get<1>(thr::get<2>(Edge(edge)));

    Coordinate edge_vec = thr::get<1>(Edge(edge));
    Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec) + get_y(edge_vec)*get_y(edge_vec));
    
    dx = std::fabs(get_x(edge_vec));
    dy = std::fabs(get_y(edge_vec));

    if(index_i > index_j){
      dx = -dx;
      dy = -dy;
    }

    A[0] = thr::get<0>(Vector4(this->_lsq_inv_iter[index_i]));
    A[1] = thr::get<1>(Vector4(this->_lsq_inv_iter[index_i]));
    A[2] = thr::get<2>(Vector4(this->_lsq_inv_iter[index_i]));
    A[3] = thr::get<3>(Vector4(this->_lsq_inv_iter[index_i]));

    A[0] += dx*dx;
    A[1] += dx*dy;
    A[2] += dx*dy;
    A[3] += dy*dy;

    this->_lsq_inv_iter[index_i] = Vector4(A[0],A[1],A[2],A[3]);

    A[0] = thr::get<0>(Vector4(this->_lsq_inv_iter[index_j]));
    A[1] = thr::get<1>(Vector4(this->_lsq_inv_iter[index_j]));
    A[2] = thr::get<2>(Vector4(this->_lsq_inv_iter[index_j]));
    A[3] = thr::get<3>(Vector4(this->_lsq_inv_iter[index_j]));
    
    A[0] += dx*dx;
    A[1] += dx*dy;
    A[2] += dx*dy;
    A[3] += dy*dy;

    this->_lsq_inv_iter[index_j] = Vector4(A[0],A[1],A[2],A[3]);

  }
};


/*****************************************************/
/*                                                   */  
/* Invert least squares matrix                       */  
/*                                                   */  
/*                                                   */  
/*---------------------------------------------------*/
/*  Input : left/right states                        */
/*  Output : interpolated left/right states          */
/*****************************************************/
struct least_sq_inv_2x2 : public thr::unary_function<Vector4,Vector4>
{

  __host__ __device__
    Vector4 operator()(const Vector4& lsq_inv) const
  {
    
    Real A[4],B[4];

    A[0] = thr::get<0>(Vector4(lsq_inv));
    A[1] = thr::get<1>(Vector4(lsq_inv));
    A[2] = thr::get<2>(Vector4(lsq_inv));
    A[3] = thr::get<3>(Vector4(lsq_inv));

    Real det = A[0]*A[3] - A[1]*A[2];
        
    Real det_inv = Real(1.0)/det;

    B[0] = A[3]*det_inv;
    B[1] = -A[2]*det_inv;
    B[2] = -A[1]*det_inv;
    B[3] = A[0]*det_inv;

    return Vector4(B[0],B[1],B[2],B[3]);

  }
};



/*****************************************************/
/* Calculate left, right and center differences for  */
/*      structured grids.                            */
/*                                                   */
/*  References:                                      */
/*    [1] J. Stone, T. Gardiner, P. Teuben,          */
/*        J. Hawley, & J. Simon "Athena: A new code  */
/*        for astrophysical MHD", ApJS, (2008).      */
/*                                                   */  
/*---------------------------------------------------*/
/*  Input : left/right states                        */
/*  Output : interpolated left/right states          */
/*****************************************************/
struct state_differences : public thr::unary_function<Interface,void>
{
  Real _gamma;
  StateIterator _state_iter;
  StateIterator _state_diff_iter;

 state_differences(Real gamma,
		   StateIterator state_iter,
		   StateIterator state_diff_iter)
   : _gamma(gamma) 
    ,_state_iter(state_iter)
    ,_state_diff_iter(state_diff_iter) {}

  __host__ __device__
    void operator()(const Interface& interface) const
  {
    Index index_i = thr::get<0>(thr::get<0>(Interface(interface)));
    Index index_j = thr::get<1>(thr::get<0>(Interface(interface)));

    // modify for outflow conditions
    if (index_i < 0) index_i = index_j;
    if (index_j < 0) index_j = index_i;    

    /* State cons_state_i = State(this->_state_iter[index_i]); */
    /* State cons_state_j = State(this->_state_iter[index_j]); */

    State state_i;
    State state_j;

    state_i = cons2prim_func(this->_gamma,State(this->_state_iter[index_i]));
    state_j = cons2prim_func(this->_gamma,State(this->_state_iter[index_j]));

    Real di = density_i;
    Real vxi = get_x(velocity_i);
    Real vyi = get_y(velocity_i);
    Real vzi = get_z(velocity_i);
    Real pgi = pressure_i;
    Real bxi = get_x(bfield_i);
    Real byi = get_y(bfield_i);
    Real bzi = get_z(bfield_i);

    Real dj = density_j;
    Real vxj = get_x(velocity_i);
    Real vyj = get_y(velocity_i);
    Real vzj = get_z(velocity_i);
    Real pgj = pressure_j;
    Real bxj = get_x(bfield_j);
    Real byj = get_y(bfield_j);
    Real bzj = get_z(bfield_j);

    Real dd,dvx,dvy,dvz,dpg,dbx,dby,dbz;

    dd = dj - di;
    dvx = vxj - vxi;
    dvy = vyj - vyi;
    dvz = vzj - vzi;
    dpg = pgj - pgi;
    dbx = bxj - bxi;
    dby = byj - byi;
    dbz = bzj - bzi;

    this->_state_diff_iter[Index(index_j)] = State(Real(dd),
						   Vector(dvx,dvy,dvz),
						   Real(dpg),
						   Vector(dbx,dby,dbz));
    
    return;
  }
};

/*****************************************************/
/* Piece-wise linear reconstruction for              */
/*      structured grids.                            */
/*                                                   */
/*  References:                                      */
/*    [1] J. Stone, T. Gardiner, P. Teuben,          */
/*        J. Hawley, & J. Simon "Athena: A new code  */
/*        for astrophysical MHD", ApJS, (2008).      */
/*                                                   */  
/*---------------------------------------------------*/
/*  Input : left/right states                        */
/*  Output : interpolated left/right states          */
/*****************************************************/
/* struct interp_plm : public thr::unary_function<Interface,void> */
/* { */
/*   Index _neigen; */
/*   Real _gamma; */
/*   StateIterator _state_iter; */
/*   StateIterator _state_diff_iter; */
/*   InterpStateIterator _interp_states_iter; */

/*  interp_plm(Index neigen, */
/* 	    Real gamma, */
/* 	    StateIterator state_iter, */
/* 	    StateIterator state_diff_iter, */
/* 	    InterpStateIterator interp_states_iter) */

/*    : _neigen(neigen)  */
/*     ,_gamma(gamma) */
/*     ,_state_iter(state_iter) */
/*     ,_state_diff_iter(state_diff_iter) */
/*     ,_interp_states_iter(interp_states_iter) {} */

/*   __host__ __device__ */
/*     void operator()(const Interface& interface) const */
/*   { */

/*     Index i,j; */
/*     Real a, a_sq; */
/*     Real ev[this->_neigen],rem[this->_neigen][this->_neigen],lem[this->_neigen][this->_neigen]; */
/*     Real dwc[this->_neigen],dwl[this->_neigen]; */
/*     Real dwr[this->_neigen],dwg[this->_neigen]; */
/*     Real dac[this->_neigen],dal[this->_neigen];   */
/*     Real dar[this->_neigen],dag[this->_neigen],da[this->_neigen]; */
/*     Real wlv[this->_neigen],wrv[this->_neigen]; */
/*     Real dw[this->_neigen],dwm[this->_neigen]; */

/*     // calculate left, right, and center differences */
/*     Index index_i = thr::get<0>(thr::get<0>(Interface(interface))); */
/*     Index index_j = thr::get<1>(thr::get<0>(Interface(interface))); */

/*     // modify for outflow conditions */
/*     if (index_i < 0) index_i = index_j; */
/*     if (index_j < 0) index_j = index_i;     */

/*     State state_i = State(this->_state_diff_iter[index_i]); */
/*     State state_j = State(this->_state_diff_iter[index_j]); */

/*     // left difference for state_i is stored at index_i */
/*     dwl[0] = density_i; */
/*     dwl[1] = get_x(velocity_i); */
/*     dwl[2] = get_y(velocity_i); */
/*     dwl[3] = get_z(velocity_i); */
/*     dwl[4] = pressure_i; */
/*     /\* dwl[5] = get_x(bfield_i); *\/ */
/*     /\* dwl[6] = get_y(bfield_i); *\/ */
/*     /\* dwl[7] = get_z(bfield_i); *\/ */

/*     // right difference for state_i is stored at index_j */
/*     dwr[0] = density_j; */
/*     dwr[1] = get_x(velocity_j); */
/*     dwr[2] = get_y(velocity_j); */
/*     dwr[3] = get_z(velocity_j); */
/*     dwr[4] = pressure_j; */
/*     /\* dwr[5] = get_x(bfield_j); *\/ */
/*     /\* dwr[6] = get_y(bfield_j); *\/ */
/*     /\* dwr[7] = get_z(bfield_j); *\/ */


/*     state_i = State(this->_state_iter[index_i]); */
/*     state_j = State(this->_state_iter[index_j]); */

/*     Real di = density_i; */
/*     Real vxi = get_x(velocity_i); */
/*     Real vyi = get_y(velocity_i); */
/*     Real vzi = get_z(velocity_i); */
/*     Real pgi = pressure_i; */
/*     Real bxi = get_x(bfield_i); */
/*     Real byi = get_y(bfield_i); */
/*     Real bzi = get_z(bfield_i); */

/*     Real dj = density_j; */
/*     Real vxj = get_x(velocity_j); */
/*     Real vyj = get_y(velocity_j); */
/*     Real vzj = get_z(velocity_j); */
/*     Real pgj = pressure_j; */
/*     Real bxj = get_x(bfield_j); */
/*     Real byj = get_y(bfield_j); */
/*     Real bzj = get_z(bfield_j); */

/*     // center difference calculated from left and right differences */
/*     dwc[0] = dwr[0] - dwl[0] + Real(2.0)*di; */
/*     dwc[1] = dwr[1] - dwl[1] + Real(2.0)*vxi; */
/*     dwc[2] = dwr[2] - dwl[2] + Real(2.0)*vyi; */
/*     dwc[3] = dwr[3] - dwl[3] + Real(2.0)*vzi; */
/*     dwc[4] = dwr[4] - dwl[4] + Real(2.0)*pgi; */
/*     /\* dwc[5] = dwr[5] - dwl[5] + Real(2.0)*bxi; *\/ */
/*     /\* dwc[6] = dwr[6] - dwl[6] + Real(2.0)*byi; *\/ */
/*     /\* dwc[7] = dwr[7] - dwl[7] + Real(2.0)*bzi; *\/ */

/*     for(i=0; i<this->_neigen; i++){ */
/*       dwg[i] = 0.0; */
/*       if (dwl[i]*dwr[i] > Real(0.0)) { */
/*         dwg[i] = Real(2.0)*dwl[i]*dwr[i]/(dwl[i]+dwr[i]); */
/*       }  */
/*     } */

/*     // set left and right eigenmatrix to zero. */
/*     for (i=0; i<this->_neigen; i++) { */
/*       for (j=0; j<this->_neigen; j++) { */
/* 	rem[i][j] = Real(0.0); */
/* 	lem[i][j] = Real(0.0); */
/*       } */
/*     } */

/*     a_sq = this->_gamma*pgi/di; */
/*     a = std::sqrt(a); */

/*     /\* Compute eigenvalues (eq. A2) *\/ */
/*     ev[0] = vxi - a; */
/*     ev[1] = vxi; */
/*     ev[2] = vxi; */
/*     ev[3] = vxi; */
/*     ev[4] = vxi + a; */
    
/*     /\* Right-eigenvectors, stored as COLUMNS (eq. A3) *\/     */
/*     rem[0][0] = Real(1.0); */
/*     rem[1][0] = -a/di; */
/*     /\*rem[2][0] = 0.0; *\/ */
/*     /\*rem[3][0] = 0.0; *\/ */
/*     rem[4][0] = a_sq; */
    
/*     rem[0][1] = Real(1.0); */
/*     /\*rem[1][1] = 0.0; *\/ */
/*     /\*rem[2][1] = 0.0; *\/ */
/*     /\*rem[3][1] = 0.0; *\/ */
/*     /\*rem[4][1] = 0.0; *\/ */
    
/*     /\*rem[0][2] = 0.0; *\/ */
/*     /\*rem[1][2] = 0.0; *\/ */
/*     rem[2][2] = Real(1.0); */
/*     /\*rem[3][2] = 0.0; *\/ */
/*     /\*rem[4][2] = 0.0; *\/ */
    
/*     /\*rem[0][3] = 0.0; *\/ */
/*     /\*rem[1][3] = 0.0; *\/ */
/*     /\*rem[2][3] = 0.0; *\/ */
/*     rem[3][3] = Real(1.0); */
/*     /\*rem[4][3] = 0.0; *\/ */
    
/*     rem[0][4] = Real(1.0); */
/*     rem[1][4] = -rem[1][0]; */
/*     /\*rem[2][4] = 0.0; *\/ */
/*     /\*rem[3][4] = 0.0; *\/ */
/*     rem[4][4] = a_sq; */
    
/*     /\* Left-eigenvectors, stored as ROWS (eq. A4) *\/     */
/*     /\*lem[0][0] = 0.0; *\/ */
/*     lem[0][1] = -half*di/a; */
/*     /\*lem[0][2] = 0.0; *\/ */
/*     /\*lem[0][3] = 0.0; *\/ */
/*     lem[0][4] = half/a_sq; */
    
/*     lem[1][0] = Real(1.0); */
/*     /\*lem[1][1] = 0.0; *\/ */
/*     /\*lem[1][2] = 0.0; *\/ */
/*     /\*lem[1][3] = 0.0; *\/ */
/*     lem[1][4] = -Real(1.0)/a_sq; */
    
/*     /\*lem[2][0] = 0.0; *\/ */
/*     /\*lem[2][1] = 0.0; *\/ */
/*     lem[2][2] = Real(1.0); */
/*     /\*lem[2][3] = 0.0; *\/ */
/*     /\*lem[2][4] = 0.0; *\/ */
    
/*     /\*lem[3][0] = 0.0; *\/ */
/*     /\*lem[3][1] = 0.0; *\/ */
/*     /\*lem[3][2] = 0.0; *\/ */
/*     lem[3][3] = Real(1.0); */
/*     /\*lem[3][4] = 0.0; *\/ */
    
/*     /\*lem[4][0] = 0.0; *\/ */
/*     lem[4][1] = -lem[0][1]; */
/*     /\*lem[4][2] = 0.0; *\/ */
/*     /\*lem[4][3] = 0.0; *\/ */
/*     lem[4][4] = lem[0][4]; */
    
/*     // differences in terms of characteristic variables */
/*     for (i=0; i<this->_neigen; i++) { */
/*       dac[i] = lem[i][0]*dwc[0]; */
/*       dal[i] = lem[i][0]*dwl[0]; */
/*       dar[i] = lem[i][0]*dwr[0]; */
/*       dag[i] = lem[i][0]*dwg[0]; */
/*       for (j=1; j<this->_neigen; j++) { */
/* 	dac[i] += lem[i][j]*dwc[j]; */
/* 	dal[i] += lem[i][j]*dwl[j]; */
/* 	dar[i] += lem[i][j]*dwr[j]; */
/* 	dag[i] += lem[i][j]*dwg[j]; */
/*       } */
/*     } */
    
/*     /\* Apply monotonicity constraints to characteristic differences *\/ */
/*     Real lim_slope1, lim_slope2, sgn_dac; */
/*     for (i=0; i<this->_neigen; i++) { */
/*       da[i] = Real(0.0); */
/*       if (dal[i]*dar[i] > Real(0.0)){ */
/* 	sgn_dac = dac[i]/std::abs(dac[i]); */
/*         lim_slope1 = fmin(    std::fabs(dal[i]),std::fabs(dar[i])); */
/*         lim_slope2 = fmin(half*std::fabs(dac[i]),std::fabs(dag[i])); */
/*         da[i] = sgn_dac*fmin(Real(2.0)*lim_slope1,lim_slope2); */
/*       } */
/*     } */

/*     // Project monotonic slopes in characteristic back to primitive variables  */
/*     for (i=0; i<this->_neigen; i++){ */
/*       dwm[i] = da[0]*rem[i][0]; */
/*       for (j=1; j<this->_neigen; j++){ */
/*         dwm[i] += da[j]*rem[i][j]; */
/*       } */
/*     } */

/*     // ensure interpolated values are bounded by neighboring cells */
/*     wlv[0] = di - half*dwm[0]; */
/*     wlv[1] = vxi - half*dwm[1]; */
/*     wlv[2] = vyi - half*dwm[2]; */
/*     wlv[3] = vzi - half*dwm[3]; */
/*     wlv[4] = pgi - half*dwm[4]; */
/*     /\* wlv[5] = bxi - half*dwm[5]; *\/ */
/*     /\* wlv[6] = byi - half*dwm[6]; *\/ */
/*     /\* wlv[7] = bzi - half*dwm[7]; *\/ */

/*     wrv[0] = di + half*dwm[0]; */
/*     wrv[1] = vxi + half*dwm[1]; */
/*     wrv[2] = vyi + half*dwm[2]; */
/*     wrv[3] = vzi + half*dwm[3]; */
/*     wrv[4] = pgi + half*dwm[4]; */
/*     /\* wrv[5] = bxi + half*dwm[5]; *\/ */
/*     /\* wrv[6] = byi + half*dwm[6]; *\/ */
/*     /\* wrv[7] = bzi + half*dwm[7]; *\/ */

/*     Real C; */

/*     C = wrv[0] + wlv[0]; */
/*     wrv[0] = fmax(fmin(di,dj),wrv[0]); */
/*     wrv[0] = fmin(fmax(di,dj),wrv[0]); */
/*     wlv[0] = C - wrv[0]; */

/*     C = wrv[1] + wlv[1]; */
/*     wrv[1] = fmax(fmin(vxi,vxj),wrv[1]); */
/*     wrv[1] = fmin(fmax(vxi,vxj),wrv[1]); */
/*     wlv[1] = C - wrv[1]; */

/*     C = wrv[2] + wlv[2]; */
/*     wrv[2] = fmax(fmin(vyi,vyj),wrv[2]); */
/*     wrv[2] = fmin(fmax(vyi,vyj),wrv[2]); */
/*     wlv[2] = C - wrv[2]; */

/*     C = wrv[3] + wlv[3]; */
/*     wrv[3] = fmax(fmin(vzi,vzj),wrv[3]); */
/*     wrv[3] = fmin(fmax(vzi,vzj),wrv[3]); */
/*     wlv[3] = C - wrv[3]; */

/*     C = wrv[4] + wlv[4]; */
/*     wrv[4] = fmax(fmin(pgi,pgj),wrv[4]); */
/*     wrv[4] = fmin(fmax(pgi,pgj),wrv[4]); */
/*     wlv[4] = C - wrv[4]; */

/*     C = wrv[5] + wlv[5]; */
/*     wrv[5] = fmax(fmin(bxi,bxj),wrv[5]); */
/*     wrv[5] = fmin(fmax(bxi,bxj),wrv[5]); */
/*     wlv[5] = C - wrv[5]; */

/*     C = wrv[6] + wlv[6]; */
/*     wrv[6] = fmax(fmin(byi,byj),wrv[6]); */
/*     wrv[6] = fmin(fmax(byi,byj),wrv[6]); */
/*     wlv[6] = C - wrv[6]; */

/*     C = wrv[7] + wlv[7]; */
/*     wrv[7] = fmax(fmin(bzi,bzj),wrv[7]); */
/*     wrv[7] = fmin(fmax(bzi,bzj),wrv[7]); */
/*     wlv[7] = C - wrv[7]; */



/*     wlv[0] = fmax(fmin(di,di),wlv[0]); */

/*     /\* for (i=0; i<this->_neigen; i++){ *\/ */
/*     /\*   C = wrv[n] + wlv[n]; *\/ */
/*     /\*   wlv[n] = fmax(fmin(pW[i][n],pW[i-1][n]),wlv[i]); *\/ */
/*     /\*   wlv[n] = fmin(fmax(pW[i][n],pW[i-1][n]),wlv[i]); *\/ */
/*     /\*   wrv[n] = C - wlv[n]; *\/ */

/*     /\*   wrv[n] = fmax(fmin(pW[i][n],pW[i+1][n]),wrv[i]); *\/ */
/*     /\*   wrv[n] = fmin(fmax(pW[i][n],pW[i+1][n]),wrv[i]); *\/ */
/*     /\*   wlv[n] = C - wrv[i]; *\/ */
/*     /\* } *\/ */

/*     /\* for (i=0; i<this->_neigen; i++){ *\/ */
/*     /\*   dw[i] = wrv[i] - wlv[i]; *\/ */
/*     /\* }     *\/ */

/*     printf("wlv[0] = %f wrv[0] = %f\n",wlv[0],wrv[0]); */

/*     State interp_state_i = State(wlv[0], */
/* 				 Vector(wlv[1],wlv[2],wlv[3]), */
/* 				 wlv[4], */
/* 				 Vector(Real(0.0),Real(0.0),Real(0.0))); */
/* 				 /\* Vector(wlv[5],wlv[6],wlv[7])); *\/ */

/*     State interp_state_j = State(wrv[0], */
/* 				 Vector(wrv[1],wrv[2],wrv[3]), */
/* 				 wrv[4], */
/* 				 Vector(Real(0.0),Real(0.0),Real(0.0))); */
/* 				 /\* Vector(wrv[5],wrv[6],wrv[7])); *\/ */

/*     this->_interp_states_iter[index_i] = InterpState(State(interp_state_j),State(interp_state_i)); */
    
/*   } */
/* }; */

/*****************************************************/
/* Convert to primitive variables                    */
/*---------------------------------------------------*/
/* Input : State ucons : tuple of cons. variables    */
/* Output : State wprim : tuple of prim. variables   */
/*****************************************************/
__host__ __device__
State cons2prim_func (Real gamma, State state)
{
  
  Real d,di,vx,vy,vz,pg,bx,by,bz,en,ke,pb;
  
  d = density;
  di = Real(1.0)/d;
  vx = get_x(momentum)*di;
  vy = get_y(momentum)*di;
  vz = get_z(momentum)*di;
  en = energy;
  bx = get_x(bfield);
  by = get_y(bfield);
  bz = get_z(bfield);
  ke = Real(0.5)*(vx*vx + vy*vy + vz*vz);
  pb = Real(0.5)*(bx*bx + by*by + bz*bz);
  
  pg = thr::max(PRESSURE_MIN,
		((gamma - Real(1.0))*(en - d*ke - pb)));
  
  return State(d,
	       Vector(vx, vy, vz),
	       Real(pg),
	       Vector(bx, by, bz));

};

/*****************************************************/
/* Convert to conservative variables                 */
/*---------------------------------------------------*/
/* Input : State wprim : tuple of prim. variables    */
/* Output : State ucons : tuple of cons. variables   */
/*****************************************************/
__host__ __device__
State prim2cons_func (Real gamma, State state)
{    
  Real d,mx,my,mz,pg,bx,by,bz,en,ke,pb;
  
  d = density;//get_x(State(state));
  mx = get_x(velocity)*d;
  my = get_y(velocity)*d;
  mz = get_z(velocity)*d;
  pg = pressure;
  bx = get_x(bfield);
  by = get_y(bfield);
  bz = get_z(bfield);
  ke = Real(0.5)*(mx*mx + my*my + mz*mz);
  pb = Real(0.5)*(bx*bx + by*by + bz*bz);
  
  en = pg/(gamma - Real(1.0)) + ke/d + pb;
  
  return State(d,
	       Vector(mx, my, mz),
	       en,
	       Vector(bx, by, bz));  
};
