/* coloring.h */

/*****************************************************/
/* Offsets for coloring                              */
/*---------------------------------------------------*/
/*****************************************************/
class Offset
{
 public:

  /* variables */

  /* Index offset; */
  Index ncolors;
  Index nboun_colors;
  Index ncolors_per_dim;
  Index iface_d;
  IndexArray faces_per_color;
  IndexArray bounds_per_color;
  InterfaceIterator face_start;
  RealIterator bn_start;
  StateIterator state_start;
  StateIterator flux_start;
  CoordinateIterator coor_start;

  Offset(); // constructor
  ~Offset(); // deconstructor

  /* functions */
  void init(Index ndim, Index ncell_x, Index ncell_y, Index nboun_x, Index nboun_y,
	    Index btype_x, Index btype_y);
  /* void reset(InterfaceIterator interface_iter, StateIterator state_iter); */
  void reset(InterfaceIterator interface_iter);
  void update(Index i);

};

Offset::Offset() {};
Offset::~Offset() {};

void Offset::init(Index ndim, Index ncell_x, Index ncell_y, Index nboun_x, Index nboun_y, 
		  Index btype_x, Index btype_y)
{

  Index nface_x = (ncell_x - Index(1))*ncell_y;
  Index nface_y = (ncell_y - Index(1))*ncell_x;
  Index nface_d;// = ncell_x*ncell_y;
  /* Index iface_d; */

  /* iface_d = 0; */
  /* if (nface_d > Index(0)) iface_d = 1; */

  nface_d = 0;
  if (iface_d > Index(0)) nface_d = ncell_x*ncell_y;

  /* ncolors_per_dim = Index(4); */
  ncolors = ncolors_per_dim*ndim;
  nboun_colors = Index(2);
  
  if(ndim < 2)
    {
      ncolors_per_dim = Index(2);
      ncolors = ncolors_per_dim*ndim;

      faces_per_color.resize(ncolors);
      bounds_per_color.resize(nboun_colors);

      Index nface = nface_x + nface_y + nboun_x + nboun_y;

      if((nface  % Index(2)) == 0) 
	{
	  faces_per_color[Index(0)] = nface / ncolors;
	  faces_per_color[Index(1)] = nface / ncolors;
	}
      else
	{ 
	  faces_per_color[Index(0)] = (nface - Index(1)) / ncolors + Index(1);
	  faces_per_color[Index(1)] = (nface - Index(1)) / ncolors;
	}
      bounds_per_color[Index(0)] = Index(2)*ncell_y;
      bounds_per_color[Index(1)] = Index(0);
    }
  else ////////////////////////// 2D
    {

      ncolors += (ncolors_per_dim + iface_d*ndim);
      faces_per_color.resize(ncolors);
      bounds_per_color.resize(nboun_colors);

      Index nface = nface_x + nface_y;

      if((ncell_x % Index(2)) == 0) // ncell_x is even 
	{
 	  faces_per_color[Index(0)] = half*ncell_x*ncell_y;
	  faces_per_color[Index(1)] = (half*ncell_x - Index(1))*ncell_y;

	  if (ncolors_per_dim > Index(2))
	    {
	      faces_per_color[Index(2)] = faces_per_color[Index(1)];
	      faces_per_color[Index(3)] = faces_per_color[Index(2)];
	      faces_per_color[Index(1)] = faces_per_color[Index(0)];
	      if((ncell_y % Index(2)) == 0) 
		{
		  faces_per_color[Index(0)] = half*Index(faces_per_color[Index(0)]);
		  faces_per_color[Index(1)] = half*Index(faces_per_color[Index(1)]);
		  faces_per_color[Index(2)] = half*Index(faces_per_color[Index(2)]);
		  faces_per_color[Index(3)] = half*Index(faces_per_color[Index(3)]);

		}
	      else
		{
		  faces_per_color[Index(0)] = half*half*(ncell_y + Index(1))*ncell_x;
		  faces_per_color[Index(1)] = half*half*(ncell_y - Index(1))*ncell_x;
		  faces_per_color[Index(2)] = half*half*(ncell_y + Index(1))*(ncell_x - Index(2.0));
		  faces_per_color[Index(3)] = half*half*(ncell_y - Index(1))*(ncell_x - Index(2.0));
		}
	    }
	}
      else // ncell_x is odd
      	{
 	  faces_per_color[Index(0)] = half*(ncell_x - Index(1))*ncell_y;
	  faces_per_color[Index(1)] = faces_per_color[0];

	  if (ncolors_per_dim > Index(2))
	    {
	      faces_per_color[Index(2)] = faces_per_color[Index(1)];
	      faces_per_color[Index(3)] = faces_per_color[Index(2)];
	      faces_per_color[Index(1)] = faces_per_color[Index(0)];
	      if((ncell_y % Index(2)) == 0) 
		{
		  faces_per_color[Index(0)] = half*Index(faces_per_color[Index(0)]);
		  faces_per_color[Index(1)] = half*Index(faces_per_color[Index(1)]);
		  faces_per_color[Index(2)] = half*Index(faces_per_color[Index(2)]);
		  faces_per_color[Index(3)] = half*Index(faces_per_color[Index(3)]);;
		}
	      else
		{
		  faces_per_color[Index(0)] = half*half*(ncell_y + Index(1))*(ncell_x - Index(1));
		  faces_per_color[Index(1)] = half*half*(ncell_y - Index(1))*(ncell_x - Index(1));
		  faces_per_color[Index(2)] = half*half*(ncell_y + Index(1))*(ncell_x - Index(1));
		  faces_per_color[Index(3)] = half*half*(ncell_y - Index(1))*(ncell_x - Index(1));
		}
	    }

      	}

      if((ncell_y % Index(2)) == 0) // ncell_y is even 
	{
 	  faces_per_color[ncolors_per_dim] = half*ncell_y*ncell_x;
	  faces_per_color[ncolors_per_dim + Index(1)] = (half*ncell_y - Index(1))*ncell_x;

	  if (ncolors_per_dim > Index(2))
	    {
	      faces_per_color[Index(6)] = faces_per_color[Index(5)];
	      faces_per_color[Index(7)] = faces_per_color[Index(6)];
	      faces_per_color[Index(5)] = faces_per_color[Index(4)];
	      if((ncell_x % Index(2)) == 0)
	      	{
		  faces_per_color[Index(4)] = half*Index(faces_per_color[Index(4)]);
		  faces_per_color[Index(5)] = half*Index(faces_per_color[Index(5)]);
		  faces_per_color[Index(6)] = half*Index(faces_per_color[Index(6)]);
		  faces_per_color[Index(7)] = half*Index(faces_per_color[Index(7)]);
	      	}
	      else
	      	{
	      	  faces_per_color[Index(4)] = half*half*(ncell_x + Index(1))*ncell_y;
	      	  faces_per_color[Index(5)] = half*half*(ncell_x - Index(1))*ncell_y;
	      	  faces_per_color[Index(6)] = half*half*(ncell_x + Index(1))*(ncell_y - Index(2));
		  faces_per_color[Index(7)] = half*half*(ncell_x - Index(1))*(ncell_y - Index(2));
	      	}
	    }


	}
      else // ncell_y is odd
      	{
 	  faces_per_color[ncolors_per_dim] = half*(ncell_y - Index(1))*ncell_x;
	  faces_per_color[ncolors_per_dim + Index(1)] = faces_per_color[ncolors_per_dim];

	  if (ncolors_per_dim > Index(2))
	    {
	      faces_per_color[Index(6)] = faces_per_color[Index(5)];
	      faces_per_color[Index(7)] = faces_per_color[Index(6)];
	      faces_per_color[Index(5)] = faces_per_color[Index(4)];
	      if((ncell_x % Index(2)) == 0) 
		{
		  faces_per_color[Index(4)] = half*Index(faces_per_color[Index(4)]);
		  faces_per_color[Index(5)] = half*Index(faces_per_color[Index(5)]);
		  faces_per_color[Index(6)] = half*Index(faces_per_color[Index(6)]);
		  faces_per_color[Index(7)] = half*Index(faces_per_color[Index(7)]);
		}
	      else
		{
		  faces_per_color[Index(4)] = half*half*(ncell_x + Index(1))*(ncell_y - Index(1));
		  faces_per_color[Index(5)] = half*half*(ncell_x - Index(1))*(ncell_y - Index(1));
		  faces_per_color[Index(6)] = half*half*(ncell_x + Index(1))*(ncell_y - Index(1));
		  faces_per_color[Index(7)] = half*half*(ncell_x - Index(1))*(ncell_y - Index(1));
		}
	    }

      	}
      

      // diagonal faces
      if (nface_d > Index(0)){
	if((ncell_x % Index(2)) == 0){ // ncell_x is even
	  faces_per_color[Index(8)] = half*ncell_x*ncell_y;
	  faces_per_color[Index(9)] = half*ncell_x*ncell_y;
	}
	else{ // ncell_x is odd
	  faces_per_color[Index(8)] = half*(ncell_x + Index(1))*ncell_y;	  
	  faces_per_color[Index(9)] = half*(ncell_x - Index(1))*ncell_y;	  
	}
      }

      if (ncolors_per_dim < Index(3)){
	faces_per_color[Index(2)*ncolors_per_dim] = Index(2)*ncell_y;
	faces_per_color[Index(2)*ncolors_per_dim + Index(1)] = Index(2)*ncell_x;
	if(btype_x == 1) faces_per_color[Index(2)*ncolors_per_dim] = Index(1)*ncell_y;
	if(btype_y == 1) faces_per_color[Index(2)*ncolors_per_dim + Index(1)] = Index(1)*ncell_x;
      }
      else if (ncolors_per_dim > Index(2))
	{
	  if((ncell_y % Index(2)) == 0)
	    {
	      faces_per_color[ncolors - Index(4)] = ncell_y;
	      faces_per_color[ncolors - Index(3)] = ncell_y;	      
	      /* if(btype_x == 1) */
	      /* 	{ */
	      /* 	  faces_per_color[Index(2)*ncolors_per_dim] = half */
	      /* 	    *faces_per_color[Index(2)*ncolors_per_dim]; */
	      /* 	  faces_per_color[Index(2)*ncolors_per_dim + Index(1)] = half */
	      /* 	    *faces_per_color[Index(2)*ncolors_per_dim + Index(1)]; */
	      /* 	} */
	    }
	  else
	    {
	      faces_per_color[ncolors - Index(4)] = ncell_y; //+ Index(1);
	      faces_per_color[ncolors - Index(3)] = ncell_y; //- Index(1);
	      /* if(btype_x == 1) */
	      /* 	{ */
	      /* 	  faces_per_color[Index(2)*ncolors_per_dim] = ncell_y + Index(1); */
	      /* 	  faces_per_color[Index(2)*ncolors_per_dim + Index(1)] = ncell_y - Index(1); */

	      /* 	  faces_per_color[Index(2)*ncolors_per_dim] = half */
	      /* 	    *faces_per_color[Index(2)*ncolors_per_dim]; */
	      /* 	  faces_per_color[Index(2)*ncolors_per_dim + Index(1)] = half */
	      /* 	    *faces_per_color[Index(2)*ncolors_per_dim + Index(1)]; */
	      /* 	} */
	    }

	  if((ncell_x % Index(2)) == 0)
	    {
	      faces_per_color[ncolors - Index(2)] = ncell_x;
	      faces_per_color[ncolors - Index(1)] = ncell_x;	      
	      /* if(btype_y == 1) */
	      /* 	{ */
	      /* 	  faces_per_color[ncolors - Index(2)] = half*Index(faces_per_color[ncolors - Index(2)]); */
	      /* 	  faces_per_color[ncolors - Index(1)] = half*Index(faces_per_color[ncolors - Index(1)]); */
	      /* 	} */
	    }
	  else
	    {
	      faces_per_color[ncolors - Index(2)] = ncell_x; // + Index(1);	      
	      faces_per_color[ncolors - Index(1)] = ncell_x; // - Index(1);	      

	      /* if(btype_y == 1) */
	      /* 	{ */
	      /* 	  faces_per_color[ncolors - Index(2)] = ncell_x + Index(1);	       */
	      /* 	  faces_per_color[ncolors - Index(1)] = ncell_x - Index(1);	       */

	      /* 	  faces_per_color[ncolors - Index(2)] = half*Index(faces_per_color[ncolors - Index(2)]); */
	      /* 	  faces_per_color[ncolors - Index(1)] = half*Index(faces_per_color[ncolors - Index(1)]); */
	      /* 	} */
	      
	    }

	}
      
      /* bounds_per_color[0] = Index(2)*ncell_y; */
      /* bounds_per_color[1] = Index(2)*ncell_x; */
      
    }
};

void Offset::reset(InterfaceIterator interface_iter)
{
  face_start = interface_iter;
  /* state_start = state_iter; */
  /* flux_start = flux_iter; */
  /* bn_start = bn_iter; */
  /* coor_start = coor_iter; */
};

void Offset::update(Index i)
{
  face_start += faces_per_color[i];
  /* state_start += faces_per_color[i]; */
  /* flux_start += faces_per_color[i]; */
  /* bn_start += faces_per_color[i]; */
  /* coor_start += faces_per_color[i]; */
};

/*****************************************************/
/* initialize edges                                  */
/*---------------------------------------------------*/
/*                                                   */
/*  nx = number of points along x-axis               */
/*  ny = number of point along y-axis                */
/*  k =  nx*ny is the number of points               */
/*                                                   */
/*  Interior edges:                                  */
/*  nedge_x = (ncell_x -1)*ncell_y, normal = (1,0)   */
/*  nedge_y = (ncell_y-1)*ncell_x, normal = (0,1)    */
/*                                                   */
/*  Total edges including boundary:                  */
/*  nedge = nedge_x + nedge_y                        */
/*          + 2*ncell_x + 2*ncell_y                  */
/*                                                   */
/*    2D edges:                                      */
/*                                                           */
/*                                                           */
/*    -----5---------6---------7---------8---------9-----    */
/*    |         |         |         |         |         |    */
/*    4   k=20  4   k=21 10   k=22  5  k=23  11   k=24  5    */
/*    |  (0,1)  |  (1,1)  |  (2,1)  |  (3,1)  |  (4,1)  |    */
/*    |         |         |         |         |         |    */
/*    -----5---------6---------7---------8---------9-----    */
/*    |         |         |         |         |         |    */
/*    4   k=15  4   k=16 10   k=17  5  k=18  11   k=19  5    */
/*    |  (0,1)  |  (1,1)  |  (2,1)  |  (3,1)  |  (4,1)  |    */
/*    |         |         |         |         |         |    */
/*    -----5---------6---------7---------8---------9-----    */
/*    |         |         |         |         |         |    */
/*    4   k=10  4   k=11 10   k=12  5  k=13  11   k=14  5    */
/*    |  (0,1)  |  (1,1)  |  (2,1)  |  (3,1)  |  (4,1)  |    */
/*    |         |         |         |         |         |    */
/*    -----17--------18--------19--------20--------21----    */
/*    |         |         |         |         |         |    */
/*    2   k=5   2   k=6   8   k=7   3  k=8    9   k=9   3    */
/*    |  (0,0)  |  (1,0)  |  (2,0)  |  (3,0)  |  (4,0)  |    */
/*    |         |         |         |         |         |    */
/*    -----12--------13--------14--------15--------16----    */
/*    |         |         |         |         |         |    */
/*    0   k=0   0   k=1   6   k=2   1  k=3    7   k=4   1    */
/*    |  (0,0)  |  (1,0)  |  (2,0)  |  (3,0)  |  (4,0)  |    */ 
/*    |         |         |         |         |         |    */
/*    -----0---------1---------2---------3---------4-----    */
/*                                                           */
/*                                                           */
/*************************************************************/
struct edges_init_2d : public thr::unary_function<Index,Edge>
{
  
  Index _color_index;
  Index _ncolors_per_dim;
  Index _ncell_x,_ncell_y;
  Index _btype_x, _btype_y;
  Index _iedge_d;
  Real _dx,_dy;
  

 edges_init_2d(Index color_index,
	       Index ncolors_per_dim,
	       Index ncell_x, Index ncell_y,
	       Index btype_x, Index btype_y, 
	       Index iedge_d,
	       Real dx, Real dy)
   : _color_index(color_index)
    ,_ncolors_per_dim(ncolors_per_dim)
    ,_ncell_x(ncell_x), _ncell_y(ncell_y)
    ,_btype_x(btype_x), _btype_y(btype_y)
    ,_iedge_d(iedge_d)
    ,_dx(dx), _dy(dy) {}

  __host__ __device__
    Edge operator()(const Index& index) const
  {

    Edge edge;

    Index nface_x, nface_y, nface, nfaces;
    Index nx = this->_ncell_x + Index(1);
    Index ny = this->_ncell_y + Index(1);
    Index edges_per_color;
    Index remainder;
    Index dim_index;
    Index color_dim_index;
    Index index_x, index_y;
    Index index_i, index_j;
    Index point_i, point_j;
    Index i;
    Index j;
    Real anx,any; // directed area vector
    Real enx,eny; // edge vector
    Real bface_i, bface_j;
    Real area_1,area_2;

    area_1 = Real(1.0);
    area_2 = Real(0.0);
    if(this->_iedge_d > Index(0)){
      area_1 = third;
      area_2 = Real(2.0)*third;
    }

    bface_i = Real(0.0);
    bface_j = Real(0.0);

    

    if(this->_ncolors_per_dim < 3){

      Index colors_per_dim = Index(2);

      dim_index = this->_color_index/colors_per_dim; // 0,1
      color_dim_index = this->_color_index % colors_per_dim; // 0,1,0,1
      if(this->_color_index < 2){ //x-faces
	  
	Index nedge = this->_ncell_x - Index(1);
	remainder = nedge % colors_per_dim;
	nedge = Index(nedge/colors_per_dim);
	if (color_dim_index < remainder) nedge += Index(1);

	index_x = colors_per_dim*(index % nedge) + color_dim_index;
	index_y = Index(index/nedge);
	
	index_i = this->_ncell_x*index_y + index_x;
	index_j = index_i + Index(1);
	
	i  = index_j % this->_ncell_x;
	j = (index_j - i)/this->_ncell_x;
	
	if(index_j < this->_ncell_x){
	      
	  point_i = i;
	}
	else{
	  point_i = index_j + Index(index_j / this->_ncell_x);
	}
	point_j = point_i + nx;		
	
	anx = Real(0.0);
	any = this->_dx;

	enx = Real(0.0);
	eny = this->_dy;


      }
      else{ //y-faces
	Index nedge = this->_ncell_x;
	
	index_x = index % nedge;
	index_y = colors_per_dim*(Index(index/nedge)) + color_dim_index;
	
	index_i = this->_ncell_x*index_y + index_x;
	index_j = index_i + this->_ncell_x;
	
	i  = index_j % this->_ncell_x;
	j = (index_j - i)/this->_ncell_x;
	
	point_i = nx*j + i;
	point_j = point_i + Index(1);		
	
	anx = this->_dy;
	any = Real(0.0);

	enx = -this->_dx;	
	eny = Real(0.0);
      }
    }

    else{ /////////////////// 4 colors

    Index colors_per_pass;
    
    colors_per_pass = Index(1);
    if (this->_ncolors_per_dim > Index(2)) colors_per_pass = Index(2);
    
    dim_index = this->_color_index/this->_ncolors_per_dim; // 0,1
    color_dim_index = this->_color_index % this->_ncolors_per_dim; // 0,1,2,3,0,1,2,3
    
    if(this->_color_index < this->_ncolors_per_dim){ ////////// x-faces
      
      Index nedge = this->_ncell_x - Index(1);
      Index remainder = nedge % (this->_ncolors_per_dim/colors_per_pass);
      nedge = Index(colors_per_pass*nedge/this->_ncolors_per_dim);
      
      if (color_dim_index < remainder*colors_per_pass) nedge += Index(1);
      
      Index offset_x;
      Index offset_y;
      offset_x = color_dim_index - (color_dim_index % colors_per_pass);
      offset_y = (color_dim_index % colors_per_pass);
      if (this->_ncolors_per_dim > Index(2)){
	
	if (color_dim_index >= this->_ncolors_per_dim/colors_per_pass){
	  offset_x = Index(1);
	}
      }
      index_x = (index % nedge)*this->_ncolors_per_dim/colors_per_pass + offset_x;	    
      index_y = Index(index/nedge)*colors_per_pass + offset_y;
      
      index_i = this->_ncell_x*index_y + index_x;
      index_j = index_i + Index(1);
	
      /* printf("index = %d color_dim_index = %d nface_x = %d index_x = %d index_y = %d index_i = %d index_j = %d\n",index,color_dim_index,nedge,index_x,index_y,index_i,index_j); */

      i  = index_j % this->_ncell_x;
      j = (index_j - i)/this->_ncell_x;
	
      if(index_j < this->_ncell_x){
	point_i = i;
      }
      else{
	point_i = index_j + Index(index_j / this->_ncell_x);
      }
      point_j = point_i + nx;		
      
      anx = -this->_dx*area_1;//third;
      any = this->_dy*area_2;//third*Real(2.0);
      
      enx = Real(0.0);
      eny = this->_dy;
      
    }
	
    else if(this->_color_index < Index(2.0)*this->_ncolors_per_dim){ ////////// y-faces
      
      Index nedge = this->_ncell_x;
      if (color_dim_index % colors_per_pass == Index(0)) nedge += nedge % colors_per_pass;
      else nedge -= nedge % colors_per_pass;
      
      Index offset_y;
      offset_y = Index(0);
      if (color_dim_index >= this->_ncolors_per_dim/colors_per_pass){
	offset_y =  Index(1);//(color_dim_index % colors_per_pass);
      }
      
      index_x = (index * colors_per_pass) % (nedge) + color_dim_index % colors_per_pass;
      index_y = Index(colors_per_pass*index/nedge)*this->_ncolors_per_dim/colors_per_pass + offset_y;
      
      index_i = this->_ncell_x*index_y + index_x;
      index_j = index_i + this->_ncell_x;
      
      /* printf("index = %d color_dim_index = %d nface_y = %d index_x = %d index_y = %d index_i = %d index_j = %d\n",index,color_dim_index,nedge,index_x,index_y,index_i,index_j); */
      
      i  = index_j % this->_ncell_x;
      j = (index_j - i)/this->_ncell_x;
	
      point_i = nx*j + i;
      point_j = point_i + Index(1);		

      // swap index
      i = point_i;
      point_i = point_j;
      point_j = i;

      anx = -this->_dy*area_2;//third*Real(2.0);
      any = this->_dx*area_1;//third;

      enx = -this->_dx;
      eny = Real(0.0);
      
    }

    else{ // diagonal faces

      colors_per_pass = Index(1);

      color_dim_index = this->_color_index % Index(2); // 0,1

      Index nedge = this->_ncell_x - Index(1);
      Index remainder = nedge % (this->_ncolors_per_dim/colors_per_pass);
      nedge = Index(colors_per_pass*nedge/this->_ncolors_per_dim);

      if (color_dim_index < remainder*colors_per_pass) nedge += Index(1);
      
      Index offset_x;
      Index offset_y;
      offset_x = color_dim_index - (color_dim_index % colors_per_pass);
      offset_y = (color_dim_index % colors_per_pass);
      if (this->_ncolors_per_dim > Index(2)){
	
      	if (color_dim_index >= this->_ncolors_per_dim/colors_per_pass){
      	  offset_x = Index(1);
      	}
      }
      index_x = (index % nedge)*this->_ncolors_per_dim/colors_per_pass + offset_x;
      index_y = Index(index/nedge)*colors_per_pass + offset_y;
      
      index_x += Index(2)*(index_y/this->_ncell_y);
      index_y = index_y % this->_ncell_y;

      index_i = this->_ncell_x*index_y + index_x;
      index_j = index_i + Index(1);
	
      index_j = index_i;

      i  = index_j % this->_ncell_x;
      j = (index_j - i)/this->_ncell_x;

      /* printf("index = %d color_dim_index = %d nface_x = %d index_x = %d index_y = %d index_i = %d index_j = %d i = %d j = %d\n",index,color_dim_index,nedge,index_x,index_y,index_i,index_j,i,j); */
	
      if(index_j < this->_ncell_x){
      	point_i = i;
      }
      else{
      	point_i = index_j + Index(index_j / this->_ncell_x);
      }
      point_j = point_i + nx + Index(1);
      
      /* point_i -= (Index(1)); */

      i = point_i;
      point_i = point_j;
      point_j = i;

      Real d = std::sqrt(this->_dx*this->_dx + this->_dy*this->_dy);

      anx = -third*this->_dy;
      any = -third*this->_dx;

      enx = -d*half*std::sqrt(Real(2.0));
      eny = -d*half*std::sqrt(Real(2.0));
      
      /* if (point_i < nx) bface_i = Real(1.0); */
      /* if (point_j >= (ny - Index(1))*nx) bface_j = Real(1.0); */
      
    }
    }// 4 colors

    return Edge(Coordinate(anx,any),Coordinate(enx,eny),IndexPair(point_i,point_j),IndexPair(index_i,index_j));

  }
};


/*****************************************************/
/*                                                   */
/*---------------------------------------------------*/
/*                                                   */
/*                                                   */
/*    -----4---------5---------6---------7-----      */
/*    |         |         x         |         |      */
/*    6   k=12  6   k=13  3  k=14   7  k=15   7      */
/*    |  (0,3)  |  (1,3)  x  (2,3)  |  (3,3)  |      */
/*    |         |         x         |         |      */
/*    -----4---------5---------6---------7-----      */
/*    |         |         x         |         |      */
/*    4   k=8   4   k=9   2  k=10   5  k=11   5      */
/*    |  (0,2)  |  (1,2)  x  (2,2)  |  (3,2)  |      */
/*    |         |         x         |         |      */
/*    xxxx 0 xxxxxxx 1 xxxxxxx 2 xxxxxxx 3 xxxx      */
/*    |         |         x         |         |      */
/*    2   k=4   2   k=5   1   k=6   3  k= 7   3      */
/*    |  (0,1)  |  (1,1)  x  (2,1)  |  (3,1)  |      */
/*    |         |         x         |         |      */
/*    -----0---------1---------2---------3-----      */
/*    |         |         x         |         |      */
/*    0   k=0   0   k=1   0   k=2   1  k= 3   1      */
/*    |  (0,0)  |  (1,0)  x  (2,0)  |  (3,0)  |      */
/*    |         |         x         |         |      */
/*    ---- 0-------- 1-------- 2-------- 3-----      */
/*                                                   */
/*                                                   */
/*****************************************************/
template<typename Tuple>
struct edge_bounds_init_2d : public thr::unary_function<Index,Tuple>
{
  
  Index _color_index;
  Index _ncolors_per_dim;
  Index _iface_d;
  Index _ncell_x,_ncell_y;
  Index _btype_x, _btype_y;
  Index _iedge_d;
  Real _dx,_dy;
  BoundaryNodeIterator _bnode_iter;
  
 edge_bounds_init_2d(Index color_index,
		     Index ncolors_per_dim,
		     Index iface_d,
		     Index ncell_x, Index ncell_y,
		     Index btype_x, Index btype_y, 
		     Index iedge_d,
		     Real dx, Real dy,
		     BoundaryNodeIterator bnode_iter)
   : _color_index(color_index)
    ,_ncolors_per_dim(ncolors_per_dim)
    ,_iface_d(iface_d)
    ,_ncell_x(ncell_x), _ncell_y(ncell_y)
    ,_btype_x(btype_x), _btype_y(btype_y)
    ,_dx(dx), _dy(dy) 
    ,_iedge_d(iedge_d)
    ,_bnode_iter(bnode_iter){}

  __host__ __device__
    Tuple operator()(const Index& index) const
  {
    Index cnt, nface_x, nface_y, nface;
    Index colors_per_pass;
    Index dim_index;
    Index colors_dim_index;
    Index index_x, index_y;
    Index index_i, index_j;
    Index faces_per_color;
    Index remainder;
    Index point_i;
    Index point_j;
    Index periodic_point_i;
    Index periodic_point_j;
    Index belem; // boundary element of adjacent face
    Index i;
    Index j;
    Index nx = this->_ncell_x + Index(1);
    Index ny = this->_ncell_y + Index(1);
    Real anx,any; // directed area vector
    Real enx,eny; // edge vector
    Real snx,sny; // scaled face outward normal    
    Real bface_i, bface_j;
    Real area_1, area_2;
    Index btype;
    Edge edge;
    BoundaryFace bface;

    colors_per_pass = Index(1);
    if (this->_ncolors_per_dim > Index(2)) colors_per_pass = Index(2);

    area_1 = half;
    area_2 = half;
    if(this->_iedge_d > Index(0)){
      area_1 = third;
      area_2 = Real(2.0)*third;
    }

    index_i = Index(-2);
    index_j = Index(-2);

    bface_i = Real(0.0);
    bface_j = Real(0.0);

    if(((this->_color_index - this->_iface_d) % Index(4)) < 1){

      btype = this->_btype_x;

      /* if(this->_btype_x == 0){ // outflow x-faces first pass */
      if(index % Index(2) == 0){ // outflow left boundary first pass

	anx = this->_dy*area_1;//third;
	any = -this->_dx*area_2;//*third*Real(2.0);
	
	enx = Real(0.0);
	eny = -this->_dy;
	
	snx = -this->_dy;
	sny = Real(0.0);
	
	index_i = Index(-1);
	index_j = colors_per_pass*half*index*this->_ncell_x;
	
	belem = index_j;
	
	i  = index_j % this->_ncell_x;
	j = (index_j - i)/this->_ncell_x;
	
	point_j = j*nx;
	point_i = point_j + nx;

	periodic_point_i = point_i + Index(this->_ncell_x);
	periodic_point_j = point_j + Index(this->_ncell_x);

	if(btype == Index(1)){
	  index_i = index_j + this->_ncell_x - Index(1);
	}

      }
      else{ // outflow right boundary first pass
	
	index_i = colors_per_pass*half*(index + Index(1))*this->_ncell_x - Index(1);
	index_j = Real(-1);
	
	belem = index_i;
	
	i  = index_i % this->_ncell_x;
	j = (index_i - i)/this->_ncell_x;
	
	point_i = index_i + j + Index(1);
	point_j = point_i + nx;
	
	anx = -this->_dy*area_1;//*third;
	any = this->_dx*area_2;//*third*Real(2.0);
	
	enx = Real(0.0);
	eny = this->_dy;
	
	snx = this->_dy;
	sny = Real(0.0);

	periodic_point_i = -Index(1);
	periodic_point_j = -Index(1);
	
	if(btype == Index(1)){
	  index_j = index_i - this->_ncell_x + Index(1);
	}

      }	
    }
    
    else if(((this->_color_index - this->_iface_d) % Index(4)) < 2 
	    && (this->_color_index- this->_iface_d) > Index(5)){

      btype = this->_btype_x;

      /* if(this->_btype_x == 0){ // outflow x-faces second pass */
	
      if(index % Index(2) == 0){ // outflow RIGHT boundary second pass
	
	index_i = colors_per_pass*half*(index + Index(1))*this->_ncell_x - Real(1.0);
	index_j = Real(-1);
	
	belem = index_i;
	
	i  = index_i % this->_ncell_x;
	j = (index_i - i)/this->_ncell_x;
	
	point_i = index_i + j + Index(1);
	point_j = point_i + nx;
	
	anx = -this->_dy*area_1;//*third;
	any = this->_dx*area_2;//*third*Real(2.0);

	enx = Real(0.0);
	eny = this->_dy;
	
	snx = this->_dy;
	sny = Real(0.0);

	periodic_point_i = -Index(1);
	periodic_point_j = -Index(1);

	if(btype == Index(1)){
	  index_j = index_i - this->_ncell_x + Index(1);
	}
	
      }
      else{ // outflow LEFT boundary second pass
	
	index_i = Index(-1);
	index_j = colors_per_pass*half*index*this->_ncell_x;
	
	belem = index_j;
	
	i  = index_j % this->_ncell_x;
	j = (index_j - i)/this->_ncell_x;
	
	point_j = j*nx;
	point_i = point_j + nx;
	
	anx = Real(0.0);
	any = -half*this->_dy;
	
	anx = this->_dy*area_1;//*third;
	any = -this->_dx*area_2;//*third*Real(2.0);

	enx = Real(0.0);
	eny = -this->_dy;
	
	snx = -this->_dy;
	sny = Real(0.0);

	periodic_point_i = point_i + Index(this->_ncell_x);
	periodic_point_j = point_j + Index(this->_ncell_x);

	if(btype == Index(1)){
	  index_i = index_j + this->_ncell_x - Index(1);
	}
	
      }
	
    }
	
    else if(((this->_color_index - this->_iface_d)% Index(4)) < 3 
	    && (this->_color_index - this->_iface_d) > Index(5)){
      
      btype = this->_btype_y;

      /* if(this->_btype_y == 0){ // outflow */
	
      if (index % Index(2) == Index(0)){ // outflow bottom boundary first pass
	
	index_i = Index(-1);
	index_j = index;
	
	belem = index_j;
	
	i  = index_j % this->_ncell_x;
	j = (index_j - i)/this->_ncell_x;
	
	point_i = i;
	point_j = point_i + Index(1);
	
	anx = this->_dy*area_2;//*third*Real(2.0);
	any = -this->_dx*area_1;//*third;
	
	enx = this->_dx;
	eny = Real(0.0);
	
	snx = Real(0.0);
	sny = -this->_dx;

	periodic_point_i = point_i + Index(this->_ncell_y)*nx;
	periodic_point_j = point_j + Index(this->_ncell_y)*nx;
	
	if(btype == Index(1)){
	  index_i = index_j + (this->_ncell_y - Index(1))*this->_ncell_x;
	}

      }
      else{ // outflow top boundary first pass
	
	index_i = this->_ncell_x*(this->_ncell_y - Index(1)) + (index % this->_ncell_x);
	index_j = Real(-1);
	
	belem = index_i;
	
	i  = index_i % this->_ncell_x;
	j = (index_i - i)/this->_ncell_x;
	
	point_j = (j+Index(1))*nx + i;
	point_i = point_j + Index(1);
	
	anx = -this->_dy*area_2;//*third*Real(2.0);
	any = this->_dx*area_1;//*third;
	
	enx = -this->_dx;
	eny = Real(0.0);
	
	snx = Real(0.0);
	sny = this->_dx;
	
	periodic_point_i = -Index(1);
	periodic_point_j = -Index(1);

	if(btype == Index(1)){
	  index_j = index_i - (this->_ncell_y - Index(1))*this->_ncell_x;
	}

      }
    }
    else{ //if((this->_color_index % Index(4)) > 2)
	  
      btype = this->_btype_y;

      /* if(this->_btype_y == 0){ // outflow */
      
      if (index % Index(2) > Index(0)) // outflow bottom boundary second pass
	{
	  index_i = Index(-1);
	  index_j = index;
	  
	  belem = index_j;
	  
	  i  = index_j % this->_ncell_x;
	  j = (index_j - i)/this->_ncell_x;
	  
	  point_i = i;
	  point_j = point_i + Index(1);
	  
	  anx = this->_dy*area_2;//*third*Real(2.0);
	  any = -this->_dx*area_1;//*third;
	  
	  enx = this->_dx;
	  eny = Real(0.0);
	  
	  snx = Real(0.0);
	  sny = -this->_dx;

	  periodic_point_i = point_i + Index(this->_ncell_y)*nx;
	  periodic_point_j = point_j + Index(this->_ncell_y)*nx;
	  
	if(btype == Index(1)){
	  index_i = index_j + (this->_ncell_y - Index(1))*this->_ncell_x;
	}

	}
      else{ // outflow top boundary second pass
	
	index_i = this->_ncell_x*(this->_ncell_y - Index(1)) + (index % this->_ncell_x);
	index_j = Real(-1);
	
	belem = index_i;
	
	i = index_i % this->_ncell_x;
	j = (index_i - i)/this->_ncell_x;
	
	point_j = (j+Index(1))*nx + i;
	point_i = point_j + Index(1);
	
	anx = -this->_dy*area_2;//*third*Real(2.0);
	any = this->_dx*area_1;//*third;
	
	enx = -this->_dx;
	eny = Real(0.0);
	
	snx = Real(0.0);
	sny = this->_dx;

	periodic_point_i = -Index(1);
	periodic_point_j = -Index(1);

	if(btype == Index(1)){
	  index_j = index_i - (this->_ncell_y - Index(1))*this->_ncell_x;
	}

	
      }
    }
    
    Index bi,bj;
    Index ndiv;

    if(point_i < nx){
      bi = point_i;
    }
    else if(point_i >= nx*this->_ncell_y){
      bi = point_i - (nx - Index(2))*(ny - Index(2)); 
    }
    else if((point_i % nx) < Index(1)){
      ndiv = point_i/nx;
      bi = point_i - (nx - Index(2))*(ndiv - Index(1.0));
    }
    else{
      ndiv = (point_i + Index(1))/nx;
      bi = point_i - (nx - Index(2))*(ndiv - Index(1));

    }

    if(point_j < nx){
      bj = point_j;
    }
    else if(point_j >= nx*this->_ncell_y){
      bj = point_j - (nx - Index(2))*(ny - Index(2)); 
    }
    else if((point_j % nx) < Index(1)){
      ndiv = point_j/nx;
      bj = point_j - (nx - Index(2))*(ndiv - Index(1.0));

    }
    else{
      ndiv = (point_j + Index(1))/nx;
      bj = point_j - (nx - Index(2))*(ndiv - Index(1));

    }

    Real px,py;

    px = get_x(thr::get<0>(BoundaryNode(this->_bnode_iter[bi])));
    py = get_y(thr::get<0>(BoundaryNode(this->_bnode_iter[bi])));

    px += half*snx;
    py += half*sny;

    this->_bnode_iter[bi] = BoundaryNode(Coordinate(px,py),point_i,btype);

    px = get_x(thr::get<0>(BoundaryNode(this->_bnode_iter[bj])));
    py = get_y(thr::get<0>(BoundaryNode(this->_bnode_iter[bj])));

    px += half*snx;
    py += half*sny;

    this->_bnode_iter[bj] = BoundaryNode(Coordinate(px,py),point_j,btype);

    /* printf("[%d][%d] bnode[%d] = %d bnode[%d] = %d\n",point_i,point_j,bi,thr::get<1>(BoundaryNode(this->_bnode_iter[bi])),bj,thr::get<1>(BoundaryNode(this->_bnode_iter[bj]))); */


    bface = BoundaryFace(Coordinate(snx,sny),IndexPair(bi,bj),
    			 IndexPair(periodic_point_i,periodic_point_j),
    			 Index(belem));
    edge = Edge(Coordinate(half*anx,half*any),Coordinate(enx,eny),
		IndexPair(point_i,point_j),IndexPair(index_i,index_j));
    
    /* printf("[%d][%d] bnode[%d] = %d bnode[%d] = %d\n",point_i,point_j,bi, */
    /* 	   thr::get<1>(BoundaryNode(this->_bnode_iter[Index(bi)])), */
    /* 	   bj,thr::get<1>(BoundaryNode(this->_bnode_iter[bj]))); */

    return Tuple(Edge(edge),BoundaryFace(bface));

  }
};
