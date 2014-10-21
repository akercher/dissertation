/* mesh_generate.h */

void Mesh::generate()
{
  
  Index i,j,k;
  Index ndof;
  Index nelem;
  Index idof;
  Real x[nx*ny], y[nx*ny]; // structured grid data

  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      k = i + j*nx;      // lexcographic ordering
      x[k] = dx*Real(i);
      y[k] = dx*Real(j);
    }
  }
  
  cells.resize(ncell());

};

/*****************************************************/
/* Connectivity                                      */
/*---------------------------------------------------*/
/*****************************************************/
struct connectivity : public thr::unary_function<Index,Quad>
{

  Index _ncell_x;

 connectivity(Index ncell_x)
   : _ncell_x(ncell_x) {}

  __host__ __device__
    Quad operator()(const Index& index) const
  {

    Index i1,i2,i3,i4;
    
    i1 = index;
    i2 = index + Index(1);
    i3 = index + this->_ncell_x + Index(1);
    i4 = index + this->_ncell_x;

    /* printf("[%d] %d %d %d %d\n",index,i1,i2,i3,i4); */
    /* return Quad(Index(i1),Index(i2),Index(i3),Index(i4)); */
    return Quad(i1,i2,i3,i4);

  }
};

/*******************************************************************/
/* Dual volume calculation                                         */
/*-----------------------------------------------------------------*/
/*******************************************************************/
/* template<typename Tuple> */
struct calc_dual_vol : public thr::unary_function<Index,Real>
{

  Index _iedge_d;
  Index _nx;
  Index _ny;
  Index _btype_x;
  Index _btype_y;
  Real _dx;
  Real _dy;

 calc_dual_vol(Index iedge_d,
	       Index nx,
	       Index ny,
	       Index btype_x,
	       Index btype_y,
	       Real dx,
	       Real dy)
   : _iedge_d(iedge_d)
    ,_nx(nx)
    ,_ny(ny)
    ,_btype_x(btype_x)
    ,_btype_y(btype_y)
    ,_dx(dx)
    ,_dy(dy) {}

  __host__ __device__
    Real operator()(const Index& index) const
  {

    Real vol;

    vol = this->_dx*this->_dy;

    if(this->_iedge_d > Index(0)){
      if (index == Index(0)) vol *= third; // bottom left corner
      else if (index == (this->_nx*(this->_ny - Index(1)))) vol *= (half*third); // top left corner
      else if (index == (this->_nx - Index(1))) vol *= (half*third); // bottom right corner
      else if (index == (this->_nx*this->_ny - Index(1))) vol *= third; // top right corner
      else if (index % (this->_nx) == 0) vol *= half; // left boundary
      else if ((index + Index(1)) % (this->_nx) == 0) vol *= half; // right boundary
      else if (index < (this->_nx)) vol *= half; // bottom boundary
      else if (index > (this->_nx*(this->_ny - Index(1)))) vol *= half; // top boundary
      /* else if (index == (this->_nx - Index(1))) vol *= half; */
      /* else if (index == ((this->_ny - Index(1))*this->_nx)) vol *= half; */
      /* else if (index == (this->_ny*this->_nx - Index(1))) vol *= half; */
    }
    else{
      if (this->_btype_x == Index(0)){
	if (index == Index(0)) vol *= (half); // bottom left corner
	else if (index == (this->_nx*(this->_ny - Index(1)))) vol *= (half); // top left corner
	else if (index == (this->_nx - Index(1))) vol *= (half); // bottom right corner
	else if (index == (this->_nx*this->_ny - Index(1))) vol *= half; // top right corner
	else if (index % (this->_nx) == 0) vol *= half; // left boundary
	else if ((index + Index(1)) % (this->_nx) == 0) vol *= half; // right boundary
	/* else if (index < (this->_nx)) vol *= half; // bottom boundary */
	/* else if (index > (this->_nx*(this->_ny - Index(1)))) vol *= half; // top boundary */
      }

      if (this->_btype_y == Index(0)){
	if (index == Index(0)) vol *= (half); // bottom left corner
	else if (index == (this->_nx*(this->_ny - Index(1)))) vol *= (half); // top left corner
	else if (index == (this->_nx - Index(1))) vol *= (half); // bottom right corner
	else if (index == (this->_nx*this->_ny - Index(1))) vol *= half; // top right corner
	/* else if (index % (this->_nx) == 0) vol *= half; // left boundary */
	/* else if ((index + Index(1)) % (this->_nx) == 0) vol *= half; // right boundary */
	else if (index < (this->_nx)) vol *= half; // bottom boundary
	else if (index > (this->_nx*(this->_ny - Index(1)))) vol *= half; // top boundary
      }
    }

    return vol;

    
  }
};

/*******************************************************************/
/* Directed area sum                                               */
/*-----------------------------------------------------------------*/
/*******************************************************************/
/* template<typename Tuple> */
struct directed_area_sum : public thr::unary_function<Edge,void>
{

  RealIterator _sum_iter;

 directed_area_sum(RealIterator sum_iter)
   
   : _sum_iter(sum_iter) {}
  
  __host__ __device__
    void operator()(const Edge& edge) const
  {
   
    Index point_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index point_j = thr::get<1>(thr::get<2>(Edge(edge)));
    
    // directed area vector
    Coordinate area_vec = thr::get<0>(Edge(edge));
    Real area_vec_mag = std::sqrt(get_x(area_vec)*get_x(area_vec)
				  + get_y(area_vec)*get_y(area_vec));
    Real area_vec_mag_inv = Real(1.0)/area_vec_mag;

    Real area_sum_i = this->_sum_iter[point_i];
    Real area_sum_j = this->_sum_iter[point_j];

    area_sum_i += area_vec_mag;
    area_sum_j -= area_vec_mag;

    this->_sum_iter[point_i] = area_sum_i;
    this->_sum_iter[point_j] = area_sum_j;
    
  }
};
