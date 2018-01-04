#ifndef __VCGLIB__TEXTCOOORD_OPTIMIZATION
#define __VCGLIB__TEXTCOOORD_OPTIMIZATION

#include <cmath>

#include <vcg/container/simple_temporary_data.h>
#include <vcg/complex/complex.h>

#include "vertex_position.h"

/*

SINGLE PATCH TEXTURE OPTIMIZATIONS

A set of classes to perform optimizations of disk->disk parametrizations.

Requires texture coords to be defined per vertex (do replicate seams!).

*/


/// Quadratic equation solver
static float QuadraticRoots(float a, float b, float c, float *x1, float *x2)
{
    float delta = b*b - 4.0f*a*c;

    if (delta >= 0) {
        *x1 = (-b + std::sqrt(delta)) / (2*a);
        *x2 = (-b - std::sqrt(delta)) / (2*a);
    }

/*
    if (delta > 0) {
        if (b >= 0) {
            *x1 = (-b - std::sqrt(delta)) / 2.0f*a;
            *x2 = 2.0f*c / (-b - std::sqrt(delta));
        } else {
            *x1 = 2.0f*c / (-b + std::sqrt(delta));
            *x2 = (-b + std::sqrt(delta)) / 2.0f*a;
        }
    }
    */
    return delta;
}

// this functions takes the 2d base coordinates and the offsets for each vertex, and returns the smallest positive
// step that prevents a triangle flip, or numeric_limits::max() if the step is unbounded
/// Negli appunti pij = (u_ij, v_ij), uguale per pki
static float ComputeStepSizeNoFlip(const Point2d& pi, const Point2d& pj, const Point2d& pk, const Point2f& di, const Point2f& dj, const Point2f& dk)
{
    Point2f dji = dj - di; Point2f dki = dk - di;
    Point2d pji = pj - pi; Point2d pki = pk - pi;

    float a = dji[0]*dki[1] - dji[1]*dki[0];
    float b = pji[0]*dki[1] + pki[1]*dji[0] - pji[1]*dki[0] - pki[0]*dji[1];
    float c = pji[0]*pki[1] - pji[1]*pki[0];

    float x1, x2;
    float delta = QuadraticRoots(a, b, c, &x1, &x2);
    if (delta < 0) return std::numeric_limits<float>::max();
    //assert(x1 != 0 && x2 != 0);
    if (x2 < x1) std::swap(x1, x2);
    if (x1 > 0) return x1;
    else if (x2 > 0) return x2;
    else return std::numeric_limits<float>::max();
}




namespace vcg
{
namespace tri
{

// helper function (checks that coords are inside -1..+1)
template <class ScalarType>
bool testParamCoordsPoint(const vcg::Point2<ScalarType> &p)
{
    ScalarType eps=(ScalarType)0.00001;
    if (!((p.X()>=-1-eps) && (p.X()<=1+eps)  && (p.Y()>=-1-eps) && (p.Y()<=1+eps)))
        return false;
    return true;
}

/* Base class for all Texture Optimizers*/
template<class MESH_TYPE>
class TexCoordOptimization{
protected:
    MESH_TYPE &m;
    SimpleTempData<typename MESH_TYPE::VertContainer, int > isFixed;
public:

    /* Types */
    typedef MESH_TYPE MeshType;
    typedef typename MESH_TYPE::VertexIterator VertexIterator;
    typedef typename MESH_TYPE::FaceIterator FaceIterator;
    typedef typename MESH_TYPE::VertexType VertexType;
    typedef typename MESH_TYPE::FaceType FaceType;
    typedef typename MESH_TYPE::ScalarType ScalarType;


    /* Access functions */
    const MeshType & Mesh() const { return m; }
    MeshType & Mesh() { return m; }

    /* Constructior */
    TexCoordOptimization(MeshType &_m) : m(_m), isFixed(_m.vert, 0)
    {
    }

    // initializes on current geometry
    virtual void TargetCurrentGeometry() = 0;

    // performs an interation. Returns largest movement.
    virtual ScalarType Iterate() = 0;

    // performs an iteration (faster, but it does not tell how close it is to stopping)
    virtual void IterateBlind() = 0;

    // performs <steps> iteration
    virtual ScalarType IterateN(int step)
    {
        for (int i=0; i<step-1; i++) this->IterateBlind();
        if (step>1) return this->Iterate();
        else return 0;
    }

    // performs iterations until convergence.
    virtual int IterateUntilConvergence(ScalarType threshold=0.0001, int maxite=5000)
    {
        int i = 0;
        while (Iterate()>threshold) {
            if (i++>maxite) return i;
        }
        return i;
    }

    // desctuctor: free temporary field
    virtual ~TexCoordOptimization()
    {
    }

    // set the current border as fixed (forced to stay in position during text optimization)
    void SetBorderAsFixed()
    {
        for (VertexIterator v=m.vert.begin(); v!=m.vert.end(); v++) {
            isFixed[v]=(v->IsB())?1:0;
        }
    }

    // everything moves, no vertex must fixed during texture optimization)
    void SetNothingAsFixed()
    {
        for (VertexIterator v=m.vert.begin(); v!=m.vert.end(); v++) {
            isFixed[v]=0;
        }
    }

    // fix a given vertex
    void FixVertex(const VertexType *v, bool fix=true)
    {
        isFixed[v]=(fix)?1:0;
    }

    bool IsFixed(const VertexType *v)
    {
        return (isFixed[v]);
    }

    bool Fixed(const FaceType* f)
    {
        return ((isFixed[f->V(0)])&&(isFixed[f->V(1)])&&(isFixed[f->V(2)]));
    }

    virtual void SetSpeed(ScalarType)
    {
    }

    virtual ScalarType GetSpeed()
    {
        return 0;
    }

    virtual void SetTheta(int)
    {
        assert(0);
    }

    virtual int GetTheta()
    {
        assert(0);
        return 0;
    }

};

/*
AREA PRESERVING TEXTURE OPTIMIZATION

as in: Degener, P., Meseth, J., Klein, R. 
       "An adaptable surface parameterization method."
       Proc. of the 12th International Meshing Roundtable, 201-213 [2003].

Features:
  
:) - Balances angle and area distortions (best results!).
:) - Can choose how to balance area and angle preservation (see SetTheta)
       theta=0 -> pure conformal (use MIPS instead!)
       theta=3 -> good balance between area and angle preservation
       theta>3 -> care more about area than about angles
:( - Slowest method.
:( - Requires a fixed boundary, else expands forever in texture space (unless theta=0).
:( - Diverges in presence of flipped faces (unless theta=0).
:( - Requires a speed parameter to be set. 
       Speed too large => when close, bounces back and forth around minimum, w/o getting any closer.
       Lower speed => longer convercence times
*/

template<class MESH_TYPE, class VertexPosFct = DefaultVertexPosition<MESH_TYPE>>
class AreaPreservingTexCoordOptimization:public TexCoordOptimization<MESH_TYPE>{
public:
  /* Types */
  typedef MESH_TYPE MeshType;
  typedef typename MESH_TYPE::VertexIterator VertexIterator;
  typedef typename MESH_TYPE::FaceIterator FaceIterator;
  typedef typename MESH_TYPE::VertexType VertexType;
  typedef typename MESH_TYPE::FaceType FaceType;
  typedef typename MESH_TYPE::ScalarType ScalarType;
  typedef typename MESH_TYPE::CoordType CoordType;
  

private:
  typedef TexCoordOptimization<MESH_TYPE> Super; // superclass (commodity)
  
  // extra data per face: [0..3] -> cotangents. [4] -> area*2
  SimpleTempData<typename MESH_TYPE::FaceContainer, Point4<ScalarType> > data;

  // vertex gradients
  SimpleTempData<typename MESH_TYPE::VertContainer, Point2<ScalarType> > sum;

  // per face gradient components (sumX is a vector of 3d points where sumX[k][i] contains
  // the gradient term of vertex i wrt to face k, same goes for sumY)
  std::vector<CoordType> sumX;
  std::vector<CoordType> sumY;

  SimpleTempData<typename MESH_TYPE::VertContainer, Point2<ScalarType> > lastDir;
  /*SimpleTempData<typename MESH_TYPE::VertContainer,omp_lock_t> lck;*/
  SimpleTempData<typename MESH_TYPE::VertContainer, ScalarType > vSpeed;
  
  ScalarType totArea;
  ScalarType speed;
  
  int theta;
  
  VertexPosFct VertexPos;

public:
  
   
  // constructor and destructor
  AreaPreservingTexCoordOptimization(MeshType &_m, VertexPosFct vpf = VertexPosFct{})
      :Super(_m),data(_m.face),sum(_m.vert),lastDir(_m.vert),vSpeed(_m.vert,1), VertexPos{vpf} {
    speed=(ScalarType)0.00005;
    theta=3;
	
	/*for (int i=0;i<m.vert.size();i++)
		omp_init_lock(&lck[i]);*/

  }
  
  ~AreaPreservingTexCoordOptimization(){
   /* data.Stop();
    sum.Stop();
    Super::isFixed.Stop();*/
  }
  
  void SetSpeed(ScalarType _speed){
    speed=_speed;
  }

  ScalarType GetSpeed(){
    return speed;
  }
  
  // sets the parameter theta:
  // good parameters are in 1..3
  //  0 = converge to pure conformal, ignore area preservation
  //  3 = good balance between area and conformal
  // >3 = area more important, angle preservation less important
  void SetTheta(int _theta){
    theta=_theta;
  }

  int GetTheta(){
    return theta;
  }
  
  void IterateBlind(){
    /* todo: do as iterate, but without */ 
    Iterate();
  }

  ScalarType Area(int i)
  {
    FaceType *f = &(Super::m.face[i]);
    double val = 0;
    if (!(Super::isFixed[f->V(0)] && Super::isFixed[f->V(1)] && Super::isFixed[f->V(2)]))
      val = (f->V(1)->T().P()-f->V(0)->T().P())^(f->V(2)->T().P()-f->V(0)->T().P());

	/* bool b0=testParamCoords(f->V(0));
	 bool b1=testParamCoords(f->V(1));
	 bool b2=testParamCoords(f->V(2));*/
			
	 if(!((fabs(val)<3.14)&&(fabs(val)>=0.0)))
	 {
		 printf("v0 %lf,%lf \n",f->V(0)->T().U(),f->V(0)->T().V());
		 printf("v1 %lf,%lf \n",f->V(1)->T().U(),f->V(1)->T().V());
		 printf("v2 %lf,%lf \n",f->V(2)->T().U(),f->V(2)->T().V());
		 printf("Area Value %lf \n",val);
		 //system("pause");
	 }

	 return fabs(val);
  }

  void InitSum()
  {
    int k;
    auto n = Super::m.vert.size();
    auto n1 = Super::m.face.size();
    for (k = 0; k < n; k++) {
      sum[k]=Point2<ScalarType>(0,0);
    }
    for (k = 0; k < n1; k++) {
      sumX[k].X()=0;
      sumX[k].Y()=0;
      sumX[k].Z()=0;
      sumY[k].X()=0;
      sumY[k].Y()=0;
      sumY[k].Z()=0;
    }
  }

  ScalarType getProjArea()
  {
    int k;
    int n = Super::m.face.size();
    ScalarType tot_proj_area=0;
   //# pragma omp parallel for
    for (k = 0; k < n; k++) {
      tot_proj_area+=Area(k);
    }
    return (tot_proj_area);
  }

vcg::Point2<ScalarType> VertValue(const int &face,const int &vert,const double &scale)
{
         FaceType *f=&Super::m.face[face];
	 /*int i=0;*/
	 vcg::Point2<ScalarType> t0=(f->V0(vert)->T().P());
     vcg::Point2<ScalarType> t1=(f->V1(vert)->T().P());
     vcg::Point2<ScalarType> t2=(f->V2(vert)->T().P());
	 ScalarType area2 = fabs((t1-t0) ^ (t2-t0));
	 ScalarType  a = (t1-t0).Norm(),
				 b =  ((t1-t0) * (t2-t0))/a,
			     c = area2 / a,
			    
				 m0= data[face][vert] / area2,
				 m1= data[face][(vert+1)%3] / area2,
				 m2= data[face][(vert+2)%3] / area2,
				  
				 mx= (b-a)/area2,
				 my= c/area2, // 1.0/a
				 mA= data[face][3]/area2* scale,
				 e = m0*((b-a)*(b-a)+c*c) + m1*(b*b+c*c) + m2*a*a, // as obvious
				 M1= mA + 1.0/mA,
				 M2= mA - 1.0/mA,
				 px= e*my,
				 py=-e*mx,
				 qx= m1*b+ m2*a,
				 qy= m1*c,

				 
				 dx=pow(M1,theta-1)
				  	 *(px*(M1+ theta*M2) - 2.0*qx*M1), 

				 dy=pow(M1,theta-1)
					   *(py*(M1+ theta*M2) - 2.0*qy*M1), 

				 gy= dy/c,
				 gx= (dx - gy*b) / a;

				  // 3d gradient
				Point2<ScalarType> val=( (t1-t0) * gx + (t2-t0) * gy ) * data[face][3]; 

				return val;
}

void UpdateSum(const double &scale)
{
         int n=Super::m.face.size();
	 int k;
	 FaceType *f;
	 ScalarType myscale=scale;
	 vcg::Point2<ScalarType> val0,val1,val2;
	  for (k=0;k<n; k++) {
                          f=&Super::m.face[k];
			  val0=VertValue(k,0,myscale);
			  val1=VertValue(k,1,myscale);
			  val2=VertValue(k,2,myscale);
			  sumX[k].V(0)=val0.X();		     
			  sumX[k].V(1)=val1.X();
			  sumX[k].V(2)=val2.X();
			  sumY[k].V(0)=val0.Y();		     
			  sumY[k].V(1)=val1.Y();
			  sumY[k].V(2)=val2.Y();
	  }
}


void SumVertex()
{
	for (unsigned int j=0; j<Super::m.face.size(); j++) 
	{
		for (int i=0;i<3;i++)
		{
                        VertexType *v=Super::m.face[j].V(i);
			sum[v].X()+=sumX[j].V(i);
			sum[v].Y()+=sumY[j].V(i);
		}
	}
}

ScalarType Iterate()
{

    InitSum();
    ScalarType tot_proj_area=getProjArea();
    double scale = tot_proj_area / totArea ;
    UpdateSum(scale);

    ScalarType max=0; // max displacement

    SumVertex();

    for (unsigned int j=0; j<Super::m.vert.size(); j++)
    {
        VertexType *v=&Super::m.vert[j];

        if (  !Super::isFixed[v] ) //if (!v->IsB())
        {
            ScalarType n=sum[v].Norm();
            //printf("N %f \n",n);
            if ( n > 1 ) { sum[v]/=n; n=1.0;}

            if (lastDir[v]*sum[v]<0) vSpeed[v]*=(ScalarType)0.85;
            else vSpeed[v]/=(ScalarType)0.92;
            lastDir[v]= sum[v];

     /* if ( n*speed<=0.1 );
      {*/

            /*vcg::Point2f goal=v->T().P()-(sum[v] * (speed * vSpeed[v]) );
            bool isOK=testParamCoordsPoint<ScalarType>(goal);
            if (isOK)
                v->T().P()-=(sum[v] * (speed * vSpeed[v]) ); */
            v->T().P()-=(sum[v] * (speed * vSpeed[v]) );


            n=n*speed * vSpeed[v];
            max=std::max(max,n);
 /* }*/
        }
    }
    return max;
}

  
  void TargetCurrentGeometry()
  {
    
   /* Super::isFixed.Start();
    data.Start();
    sum.Start();*/
    for (auto& v : Super::m.vert) lastDir[&v] = Point2f(0, 0);
    sumX.resize(Super::m.face.size());
    sumY.resize(Super::m.face.size());
    totArea=0;
    for (FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++) {
      //double area2 = ((f->V(1)->P() - f->V(0)->P() )^(f->V(2)->P() - f->V(0)->P() )).Norm();
      double area2 = ((VertexPos(&*f, 1) - VertexPos(&*f, 0)) ^ (VertexPos(&*f, 2) - VertexPos(&*f, 1))).Norm();
      totArea+=area2;
		  //if (  Super::isFixed[f->V1(0)] )
      for (int i=0; i<3; i++) {
        //data[f][i]=((f->V1(i)->P() - f->V0(i)->P()) * (f->V2(i)->P() - f->V0(i)->P())) / area2;
        data[f][i]= ((VertexPos.V1(&*f, i) - VertexPos.V0(&*f, i)) * (VertexPos.V2(&*f, i) - VertexPos.V0(&*f, i))) / area2;
        data[f][3]=area2;
      }
    }
  }
  
};


/*
MIPS TEXTURE OPTIMIZATION

Features:
  
:( - Targets angle distortions only (not ideal for texture mapping).
:) - Quite fast.
:) - Does not require fixed boundary (will auto-find a boundary -- up to a scale).
:) - Tends to nicely heal flipped faces 
:( - Requires a speed parameter to be set 
    (SHOULD NOT BE LIKE THIS. BETTER IMPLEMENTATION NEEDED). 
       Speed too large => when close, bounces back and forth around minimum, w/o getting any closer.
       Lower speed => longer convercence times

*/

template<class MESH_TYPE, class VertexPosFct = DefaultVertexPosition<MESH_TYPE>>
class MIPSTexCoordOptimization:public TexCoordOptimization<MESH_TYPE>{
public:
    /* Types */
    typedef MESH_TYPE MeshType;
    typedef typename MESH_TYPE::VertexIterator VertexIterator;
    typedef typename MESH_TYPE::FaceIterator FaceIterator;
    typedef typename MESH_TYPE::VertexType VertexType;
    typedef typename MESH_TYPE::FaceType FaceType;
    typedef typename MESH_TYPE::ScalarType ScalarType;


protected:
    typedef TexCoordOptimization<MESH_TYPE> Super; // superclass (commodity)

    // extra data per face: [0..3] -> cotangents.
    SimpleTempData<typename MESH_TYPE::FaceContainer, Point4<ScalarType> > data;
    SimpleTempData<typename MESH_TYPE::FaceContainer, Point3<ScalarType> > cot;
    SimpleTempData<typename MESH_TYPE::VertContainer, Point2<double> > sum;
    SimpleTempData<typename MESH_TYPE::VertContainer, Point2<double> > lastDir;
  SimpleTempData<typename MESH_TYPE::VertContainer, ScalarType > vSpeed;

    ScalarType totArea;
    ScalarType speed;

    VertexPosFct VertexPos;

public:


    // constructor and destructor
    MIPSTexCoordOptimization(MeshType &_m, VertexPosFct vpf = VertexPosFct{})
        : Super(_m),data(_m.face), cot(_m.face), sum(_m.vert), lastDir(_m.vert), vSpeed(_m.vert, 0.1), VertexPos{vpf}
    {
        speed = (ScalarType)0.001;
    }

    ~MIPSTexCoordOptimization()
    {
        /* data.Stop();
        sum.Stop();
        Super::isFixed.Stop();*/
    }

    void SetSpeed(ScalarType _speed)
    {
        speed=_speed;
    }

    ScalarType GetSpeed()
    {
        return speed;
    }

    void IterateBlind()
    {
        /* todo: do as iterate, but without */
        Iterate();
    }

    ScalarType Iterate()
    {

#define v0 (f.V(0)->T().P())
#define v1 (f.V(1)->T().P())
#define v2 (f.V(2)->T().P())
#define vi (f.V(i)->T().P())
#define vj (f.V(j)->T().P())
#define vk (f.V(k)->T().P())

        SimpleTempData<typename MESH_TYPE::VertContainer, Point2<ScalarType>> grad(Super::m.vert);
        for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) {
            //sum[v].Zero();
            sum[v]=Point2<double>(0,0);
            grad[v] = Point2<ScalarType>(0,0);
        }

        int a1 = 0;
        int a2 = 0;
        ScalarType totProjArea=0;
        ScalarType E_curr = 0;



        for (auto& f : Super::m.face) {
            ScalarType area2 = ((v1-v0) ^ (v2-v0));

            area2 > 0 ? a1++ : a2++;

            totProjArea+=area2;
            ScalarType o[3] = { // (opposite edge)^2
                (   v1-v2).SquaredNorm(),
                (v0   -v2).SquaredNorm(),
                (v0-v1   ).SquaredNorm(),
            };
            ScalarType e = (data[f][0] * o[0] +
                            data[f][1] * o[1] +
                            data[f][2] * o[2] ) / (area2*area2);
            //E_curr += (e / area2);



            ScalarType cotg[3] = {
                 cot[f][0],
                 cot[f][1],
                 cot[f][2]
             };

             ScalarType area3D = data[f][3];
             ScalarType areaUV = 0.5f * ((v1 - v0) ^ (v2 - v0));

             // dirichlet energy at the current face
             ScalarType e_d = cotg[0] * o[0] + cotg[1] * o[1] + cotg[2] * o[2];

             E_curr += (area3D*0.5f*e_d) / areaUV;


            for (int i=0; i<3; i++) {
                int j=(i+1)%3, k=(i+2)%3;
                ScalarType p=(vj-vi)*(vk-vi);
                ScalarType gy= e*(o[k]-p) - 2*data[f][j];
                ScalarType gx= e*(o[j]-p) - 2*data[f][k];

                // speed free mode: (a try!)
                //sum[f->V(i)]+= ( (vj-vi) * gx + (vk-vi) * gy );// / area2;

                // speed mode:
                sum[f.V(i)]+= ( (vj-vi) * gx + (vk-vi) * gy ) / area2;


                 ScalarType gu_cotan_term = - cotg[k]*vj.X() - cotg[j]*vk.X() + (cotg[k]+cotg[j])*vi.X();
                 ScalarType gv_cotan_term = - cotg[k]*vj.Y() - cotg[j]*vk.Y() + (cotg[k]+cotg[j])*vi.Y();

                 ScalarType gu_area = ( vj.Y() - vk.Y());
                 ScalarType gv_area = (-vj.X() + vk.X());

                 ScalarType gu = area3D * (gu_cotan_term*areaUV - 0.5f*e_d*gu_area) / std::pow(areaUV, 2);
                 ScalarType gv = area3D * (gv_cotan_term*areaUV - 0.5f*e_d*gv_area) / std::pow(areaUV, 2);

                 grad[f.V(i)].X() += gu;
                 grad[f.V(i)].Y() += gv;


                //grad[f.V(i)] += Point2<ScalarType>(gx, gy);
            }
        }

        //if (a1 > 0 && a2 > 0) std::cout << "Parameterization contains folds " << a1 << " " << a2 << std::endl;
        //std::cout << "Energy at the current iteration is " << E_curr << std::endl;


/*
        // Compute descent step with inexact line search
        float c = 0.1f;
        float scalarTerm = 0;
        float gamma = 0.8f;
        float E_t = E_curr; // energy value for step t
        // store current positions and precompute decrease factor (p_k.dot(grad_Ek) == -grad_Ek.SquaredNorm() since
        // the descent direction p_k is -grad_Ek)

        /// TO TRY normalize each direction
        SimpleTempData<typename MESH_TYPE::VertContainer, Point2<ScalarType>> currentUV(Super::m.vert);
        for (auto& v : Super::m.vert) {
            currentUV[v] = v.T().P();
            scalarTerm += (-sum[v]) * grad[v];
        }
        scalarTerm = c * scalarTerm;

        // Compute max stepsize that guarantees no triangle inversion
        float t = std::numeric_limits<float>::max();
        for (auto& f : Super::m.face) {
            float tFace = ComputeStepSizeNoFlip(v0, v1, v2, -sum[f.V(0)], -sum[f.V(1)], -sum[f.V(2)]);
            if (tFace < t) t = tFace;
        }
        std::cout << "initial t = " << t << std::endl;
        t = std::min(1.0f, t*0.8f);

        int numIter = 0;
        float maxMovement;
        while (E_t > E_curr + t * scalarTerm) {
            numIter++;
            maxMovement = 0;
            E_t = 0.0f;
            for (auto& v : Super::m.vert) {
                v.T().P() = currentUV[v] - t*sum[v];
                std::cout << t << "    " << sum[v].Norm() << "       " << t*sum[v].Norm() << std::endl;
                maxMovement = std::max(maxMovement, (t*sum[v]).Norm());
            }

            for (auto& f : Super::m.face) {
                ScalarType area2 = ((v1-v0) ^ (v2-v0));

                ScalarType o[3] = { // (opposite edge)^2
                    (   v1-v2).SquaredNorm(),
                    (v0   -v2).SquaredNorm(),
                    (v0-v1   ).SquaredNorm(),
                };
                ScalarType e = (data[f][0] * o[0] +
                                data[f][1] * o[1] +
                                data[f][2] * o[2] ) / (area2*area2);

                //E_t += (e / area2);
                ScalarType cotg[3] = {
                    cot[f][0],
                    cot[f][1],
                    cot[f][2]
                };

                ScalarType area3D = data[f][3];
                ScalarType areaUV = 0.5f * ((v1 - v0) ^ (v2 - v0));

                // dirichlet energy at the current face
                ScalarType e_d = cotg[0] * o[0] + cotg[1] * o[1] + cotg[2] * o[2];

                E_t += (area3D*0.5f*e_d) / areaUV;

            }

            // update step size for next iteration
            t *= gamma;
        }

        //std::cout << "Energy went from " << totalEnergy << " to " << E_t << " ITERATIONS=" << numIter
        //          << " STEPSIZE=" << t << std::endl;

        return maxMovement;
*/

        ScalarType max=0; // max displacement

        for (unsigned int j=0; j<Super::m.vert.size(); j++)
        {
            VertexType *v=&Super::m.vert[j];

            if (  !Super::isFixed[v] ) //if (!v->IsB())
            {
                ScalarType n=sum[v].Norm();
                //printf("N %f \n",n);
                if ( n > 0.001f ) { sum[v]/= (n/0.001f); n=0.001f;}

                if (lastDir[v]*sum[v]<0) vSpeed[v]*=(ScalarType)0.85;
                else vSpeed[v]/=(ScalarType)0.92;
                lastDir[v]= sum[v];

                v->T().P()-=(sum[v] * (speed * vSpeed[v]) );


                n=n*speed * vSpeed[v];
                max=std::max(max,n);
            }
        }
        return E_curr;




/*
        for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) {
            if (!Super::isFixed[v]) {
                // speed free mode: (a try!)
                //v->T().P()-=speed * sum[v] *totProjArea/totArea;

                // speed mode:

                ScalarType n=sum[v].Norm();
                if ( n > 1 ) { sum[v]/=n; n=1.0; }
                v->T().P()-=(sum[v] ) * speed ;
                if (max<n) max=n;
            }
        }*/
        //return max;
#undef v0
#undef v1
#undef v2
#undef vi
#undef vj
#undef vk
        //printf("rejected %d\n",rejected);
    }
/*
    void TargetCurrentGeometry()
    {
        totArea = 0;
        for (FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++) {
            double area2 = ((f->V(1)->P() - f->V(0)->P()) ^ (f->V(2)->P() - f->V(0)->P())).Norm();
            totArea += area2;
            for (int i=0; i<3; i++) {
                data[f][i] = std::max(vcg::Angle(f->V1(i)->P() - f->V0(i)->P(),
                                                 f->V2(i)->P() - f->V0(i)->P()),
                                      float(1e-8));
                // / area2;
            }
        }
    }
*/
    void TargetCurrentGeometry()
    {
        /* Super::isFixed.Start();
        data.Start();
        sum.Start();*/
        for (auto& v : Super::m.vert) lastDir[v].SetZero();
        totArea = 0;
        for (FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++) {
            double area2 = ((f->V(1)->P() - f->V(0)->P()) ^ (f->V(2)->P() - f->V(0)->P())).Norm();
            data[f][3] = area2 / 2.0;
            totArea += area2;
            for (int i=0; i<3; i++) {
                //data[f][i]=((f->V1(i)->P() - f->V0(i)->P()) * (f->V2(i)->P() - f->V0(i)->P()));
                data[f][i] = (VertexPos.V1(&*f, i) - VertexPos.V0(&*f, i)) * (VertexPos.V2(&*f, i) - VertexPos.V0(&*f, i));
                // / area2;

                Point3f a = VertexPos.V1(&*f, i) - VertexPos.V0(&*f, i);
                                Point3f b = VertexPos.V2(&*f, i) - VertexPos.V1(&*f, i);
                                Point3f c = VertexPos.V2(&*f, i) - VertexPos.V0(&*f, i);
                                float cotg = (a.SquaredNorm() + c.SquaredNorm() - b.SquaredNorm())/(2.0f*data[f][3])/2.0f;
                                cot[f][i] = cotg;

            }
        }
    }

};

template<class MESH_TYPE, class VertexPosFct = DefaultVertexPosition<MESH_TYPE>>
class SymmetricDirichletTexCoordOptimization:public TexCoordOptimization<MESH_TYPE>{
public:
    /* Types */
    typedef MESH_TYPE MeshType;
    typedef typename MESH_TYPE::VertexIterator VertexIterator;
    typedef typename MESH_TYPE::FaceIterator FaceIterator;
    typedef typename MESH_TYPE::VertexType VertexType;
    typedef typename MESH_TYPE::FaceType FaceType;
    typedef typename MESH_TYPE::ScalarType ScalarType;


protected:
    typedef TexCoordOptimization<MESH_TYPE> Super; // superclass (commodity)

    // extra data per face: [0..3] -> cotangents, 3D area
    SimpleTempData<typename MESH_TYPE::FaceContainer, Point4<ScalarType> > data;
    SimpleTempData<typename MESH_TYPE::VertContainer, Point2<ScalarType> > sum;

    VertexPosFct VertexPos;

public:

    // constructor and destructor
    SymmetricDirichletTexCoordOptimization(MeshType &_m, VertexPosFct vpf = VertexPosFct{}):Super(_m),data(_m.face),sum(_m.vert),VertexPos{vpf}
    {
    }

    ~SymmetricDirichletTexCoordOptimization()
    {
    }

    void IterateBlind()
    {
        /* todo: do as iterate, but without */
        Iterate();
    }


    ScalarType Iterate()
    {

//#define p0 (VertexPos(&*f, 0))
//#define p1 (VertexPos(&*f, 1))
//#define p2 (VertexPos(&*f, 2))
#define p0 (f.V(0)->P())
#define p1 (f.V(1)->P())
#define p2 (f.V(2)->P())
#define q0 (f.V(0)->T().P())
#define q1 (f.V(1)->T().P())
#define q2 (f.V(2)->T().P())
#define qi (f.V(i)->T().P())
#define qj (f.V(j)->T().P())
#define qk (f.V(k)->T().P())

        for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) {
            sum[v]=Point2<ScalarType>(0,0);
        }

        ScalarType totalEnergy = 0;
        for (auto& f : Super::m.face) {
            ScalarType o[3] = { // (opposite edge)^2
                (   q1-q2).SquaredNorm(),
                (q0   -q2).SquaredNorm(),
                (q0-q1   ).SquaredNorm(),
            };

            ScalarType cotg[3] = {
                std::tan(float(M_PI_2) - data[f][0]),
                std::tan(float(M_PI_2) - data[f][1]),
                std::tan(float(M_PI_2) - data[f][2]),
            };

            ScalarType area3D = data[f][3];
            ScalarType areaUV = 0.5f * ((q1 - q0) ^ (q2 - q0));

            // dirichlet energy at the current face
            ScalarType e_d = 0.5f*(cotg[0] * o[0] + cotg[1] * o[1] + cotg[2] * o[2]);

            ScalarType areaTermRatio = - (area3D*area3D) / (areaUV*areaUV*areaUV);
            ScalarType e_a = (1 + (area3D*area3D) / (areaUV*areaUV));

            //totalEnergy += e_d / (areaUV); // multiplied by the area in parameter space
            totalEnergy += (1 + (area3D*area3D)/(areaUV*areaUV)) * (e_d);

            assert(areaUV > 0);

            for (int i = 0; i < 3; ++i) {
                int j = (i+1)%3;
                int k = (i+2)%3;

                ScalarType gu_area = areaTermRatio * ( qj.Y() - qk.Y());
                ScalarType gv_area = areaTermRatio * (-qj.X() + qk.X());

                ScalarType gu_angle = (cotg[k]*(qi.X() - qj.X()) + cotg[j]*(qi.X() - qk.X()));
                ScalarType gv_angle = (cotg[k]*(qi.Y() - qj.Y()) + cotg[j]*(qi.Y() - qk.Y()));

                ScalarType gu = gu_area * e_d + e_a * gu_angle;
                ScalarType gv = gv_area * e_d + e_a * gv_angle;

                sum[f.V(i)].X() += gu;
                sum[f.V(i)].Y() += gv;

            }
        }

        // Compute descent step with inexact line search
        float c = 0.1f;
        float scalarTerm = 0;
        float gamma = 0.8f;
        float E_t = totalEnergy; // energy value for step t
        // store current positions and precompute decrease factor (p_k.dot(grad_Ek) == -grad_Ek.SquaredNorm() since
        // the descent direction p_k is -grad_Ek)

        /// TO TRY normalize each direction
        SimpleTempData<typename MESH_TYPE::VertContainer, Point2d> currentUV(Super::m.vert);
        for (auto& v : Super::m.vert) {
            //std::cout << sum[v].X() << " , " << sum[v].Y() << std::endl;
            currentUV[v] = v.T().P();
            scalarTerm += sum[v].SquaredNorm();
        }
        scalarTerm = - c * scalarTerm;

        // Compute max stepsize that guarantees no triangle inversion
        float t = std::numeric_limits<float>::max();
        for (auto& f : Super::m.face) {
            float tFace = ComputeStepSizeNoFlip(q0, q1, q2, -sum[f.V(0)], -sum[f.V(1)], -sum[f.V(2)]);
            if (tFace < t) t = tFace;
        }
        std::cout << "initial t = " << t << std::endl;
        t = std::min(1.0f, t*0.8f);

        int numIter = 0;
        while (E_t > totalEnergy + t * scalarTerm) {
            numIter++;
            t *= gamma;
            E_t = 0.0f;
            for (auto& v : Super::m.vert) {
                Point2d sv = Point2d(sum[v].X(), sum[v].Y());
                v.T().P() = currentUV[v] - double(t)*sv;
            }

            // compute energy
            for (auto& f : Super::m.face) {
                ScalarType o[3] = { // (opposite edge)^2
                    (   q1-q2).SquaredNorm(),
                    (q0   -q2).SquaredNorm(),
                    (q0-q1   ).SquaredNorm(),
                };
                ScalarType cotg[3] = {
                    std::tan(float(M_PI_2) - data[f][0]),
                    std::tan(float(M_PI_2) - data[f][1]),
                    std::tan(float(M_PI_2) - data[f][2]),
                };
                ScalarType areaUV = 0.5f * ((q1 - q0) ^ (q2 - q0));
                ScalarType area3D = data[f][3];
                ScalarType e_d = 0.5f * (cotg[0] * o[0] + cotg[1] * o[1] + cotg[2] * o[2]);
                E_t += (1 + (area3D*area3D)/(areaUV*areaUV)) * (e_d);
            }

            std::cout << E_t << std::endl;

        }

        std::cout << "Energy went from " << totalEnergy << " to " << E_t << " ITERATIONS=" << numIter
                  << " STEPSIZE=" << t << std::endl;

        return std::sqrt((scalarTerm / (-c)));// returns the magnitude of the gradient
#undef p0
#undef p1
#undef p2
#undef q0
#undef q1
#undef q2
#undef qi
#undef qj
#undef qk
    }

    void TargetCurrentGeometry()
    {
        for (auto& f : Super::m.face) {
            double area3D = ((VertexPos(&f, 1) - VertexPos(&f, 0)) ^ (VertexPos(&f, 2) - VertexPos(&f, 0))).Norm() / 2.0f;
            data[f][3] = area3D;
            for (int i=0; i<3; i++) {
                data[f][i] = std::max(vcg::Angle(VertexPos.V1(&f, i) - VertexPos.V0(&f, i),
                                                 VertexPos.V2(&f, i) - VertexPos.V0(&f, i)),
                                      float(1e-8));
            }
        }
    }
};

#if 0  // Temporarly commented out. It still have to be thoroughly tested...

template<class MESH_TYPE>
class WachspressTexCoordOptimization:public TexCoordOptimization<MESH_TYPE>{
public:
  /* Types */
  typedef MESH_TYPE MeshType;
  typedef typename MESH_TYPE::VertexIterator VertexIterator;
  typedef typename MESH_TYPE::FaceIterator FaceIterator;
  typedef typename MESH_TYPE::VertexType VertexType;
  typedef typename MESH_TYPE::FaceType FaceType;
  typedef typename MESH_TYPE::VertexType::TexCoordType::PointType PointType;
  typedef typename MESH_TYPE::VertexType::TexCoordType::PointType::ScalarType ScalarType;
  
private:
  class Factors{
    public:
    ScalarType data[3][2];
  };
  
  typedef TexCoordOptimization<MESH_TYPE> Super; // superclass (commodity)
  
  // extra data per face: [0..3] -> cotangents. [4] -> area*2
  SimpleTempData<typename MESH_TYPE::FaceContainer, Factors > factors;
  
public:
  
   
  // constructor and destructor
  WachspressTexCoordOptimization(MeshType &_m):Super(_m),factors(_m.face){
  }
  
  ~WachspressTexCoordOptimization(){
   /* factors.Stop();
    Super::isFixed.Stop();*/
  }

  void IterateBlind(){
    /* todo: do as iterate, but without */ 
    Iterate();
  }
    
  ScalarType Iterate(){
    
    ScalarType max; // max displacement
  	
	  
    for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) {
		  v->div=0; v->sum.SetZero();
	  }

    #define vi0 (f->V(i0)->P())
    #define vi1 (f->V(i1)->P())
    #define vi2 (f->V(i2)->P())
    #define EPSILON 1e-4

	  for (FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++){
		  ScalarType A=(((f->V(0)->P()) - (f->V(1)->P()))^((f->V(0)->P()) - (f->V(2)->P()))).Norm(); 
		  if (A<EPSILON) continue;
		  for (int i0=0; i0<3; i0++) {
			  int i1=(i0+1)%3,i2=(i0+2)%3;


  		  ScalarType fact = (vi1-vi0)*(vi2-vi0)/A;
				//fact=1;
			  if ( (!f->V(i1)->IsB()) /*|| f->V(o)->IsB()*/);{
				  f->V(i1)->sum += f->V(i0)->T().P() * fact;
				  f->V(i1)->div += fact;
			  }
			  if ( (!f->V(i2)->IsB()) /*|| f->V(o)->IsB()*/);{
				  f->V(i2)->sum += f->V(i0)->T().P() * fact;
				  f->V(i2)->div += fact;
			  }
		  }
	  }

	  for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) {
      if (  !Super::isFixed[v] )
		  if (v->div>0.001) {
			  v->T().P() = v->sum/v->div;
		  }
	  }
    return max; 	
  }
  
  
  void TargetCurrentGeometry(){
  }
  
};

#endif

#if 0

template<class MESH_TYPE> 
class MeanValueTexCoordOptimization:public TexCoordOptimization<MESH_TYPE>{
public:
  /* Types */
  typedef MESH_TYPE MeshType;
  typedef typename MESH_TYPE::VertexIterator VertexIterator;
  typedef typename MESH_TYPE::FaceIterator FaceIterator;
  typedef typename MESH_TYPE::VertexType VertexType;
  typedef typename MESH_TYPE::FaceType FaceType;
  typedef typename MESH_TYPE::VertexType::TexCoordType::PointType PointType;
  typedef typename MESH_TYPE::VertexType::TexCoordType::PointType::ScalarType ScalarType;
  
private:
  class Factors{
    public:
    ScalarType data[3][2];
  };
  
  typedef TexCoordOptimization<MESH_TYPE> Super; // superclass (commodity)
  
  // extra data per face: factors
  SimpleTempData<typename MESH_TYPE::FaceContainer, Factors > factors;
  
  // extra data per vertex: sums and div
  SimpleTempData<typename MESH_TYPE::VertContainer, PointType > sum;
  SimpleTempData<typename MESH_TYPE::VertContainer, ScalarType > div;
  
public:
  
   
  // constructor and destructor
  MeanValueTexCoordOptimization(MeshType &_m):Super(_m),factors(_m.face),sum(_m.vert),div(_m.vert){
  }
  
  ~MeanValueTexCoordOptimization(){
   /* factors.Stop();
    sum.Stop();
    div.Stop();
    Super::isFixed.Stop();*/
  }
  
  void IterateBlind(){
    /* todo: do as iterate, but without */ 
    Iterate();
  }
  
  ScalarType Iterate(){
    
    ScalarType max=0; // max displacement
  		  
	  for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) {
		  sum[v]=PointType(0,0);
		  div[v]=0; 
	  }

	  for (FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++){
		  for (int i=0; i<3; i++) 
		  for (int j=1; j<3; j++) {
			  int d=i, o=(i+3-j)%3;
			  sum[f->V(d)] += f->V(o)->T().P() * factors[f].data[i][j-1];
			  div[f->V(d)] += factors[f].data[i][j-1];
		  }
	  }

	  for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) 
    if (  !Super::isFixed[v] )
	  if (		div[v]>0.000001 ) {
		  PointType swap=v->T().P();
		  PointType goal=sum[v]/div[v];
		
		  v->T().P() = goal*(ScalarType)0.1+swap*(ScalarType)0.9;

		  //v->T().P()=v->RestUV*(1-v->Damp)+(sum[v]/div[v])*(v->Damp);
		  ScalarType temp=(swap-v->T().P()).SquaredNorm();
		  if (max<temp)
			  max=temp;
	  }
    return max; 	
  }
  
  void TargetEquilateralGeometry(){
	  for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) {
		  div[v]=0;
	  }
	  const ScalarType fact= 1.0 / sqrt(3.0);
	  for (FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++){
		  for (int i=0; i<3; i++) 
			  for (int j=0; j<2; j++) {
				  factors[f].data[i][j] = fact;
				  div[f->V(i)] += fact ;
			  }
	  }
  }

  void TargetCurrentGeometry(){
    
    /*Super::isFixed.Start();
    factors.Start();
    sum.Start();
    div.Start();*/

    for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) {
		  div[v]=0; 
	  };
	  for (FaceIterator  f=Super::m.face.begin(); f!=Super::m.face.end(); f++){
		  for (int i=0; i<3; i++) 
		  for (int j=1; j<3; j++) factors[f].data[i][j-1]=0;
	  };

    #define vs (f->V(s)->P())
    #define vd (f->V(d)->P())
    #define vo (f->V(o)->P())
    #define EPSILON 1e-4

	  for (FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++){
		  int s=0,d=1,o=2;
		  ScalarType A=((vs - vd)^(vs - vo)).Norm(); 
		  if (A<EPSILON) break;
		  for (int i=0; i<3; i++) 
		  for (int j=1; j<3; j++) {
			  d=i; s=(i+j)%3; o=(i+3-j)%3;
			  {

				  ScalarType dd=((vd-vs).Norm());
				  if (dd<=0.0001) continue;
				  ScalarType fact= ( ( vd -vo ).Norm() - ((vd-vo)*(vd-vs))/dd) /A;

				  //if (fact<0) printf("AAAGH!");
				  factors[f].data[d][j-1] = fact;
				  //f->V(d)->sum += f->V(o)->projected * fact;
				  div[f->V(d)] += fact ;
			  } //else {
				  //printf(".");
			  //}
		  }
	  }
	  for (FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++){
      for (int i=0; i<3; i++) 
		  for (int j=1; j<3; j++)   
		  	if (div[f->V(i)]>0.0001) {
		  	  //fact[i][j-1]/=div[f->V(i)];
		  	} 
		  	//else f->fact[i][j-1]=0.0;      
    }
		
	  /*
    for (f=face.begin(); f!=face.end(); f++)  {
	  	for (int i=0; i<3; i++) 
		  for (int j=1; j<3; j++) 
		  	if (f->V(i)->div>0.01) {
		  	  f->fact[i][j-1]/=f->V(i)->div;
		  	} 
		  	else f->fact[i][j-1]=0.0;
	  } */
	  
	  #undef vs 
    #undef vd 
    #undef vo

  }
  
};
#endif

/* texture coords general utility functions */
/*++++++++++++++++++++++++++++++++++++++++++*/

// returns false if any fold is present (faster than MarkFolds)
template<class MESH_TYPE>
bool IsTexCoordFoldFree(MESH_TYPE &m){
  
  assert(HasPerVertexTexCoord(m));
  
  typedef typename MESH_TYPE::VertexType::TexCoordType::PointType::ScalarType ScalarType;
  
  ScalarType lastsign=0;
  for (typename MESH_TYPE::FaceIterator f=m.face.begin(); f!=m.face.end(); f++){
    ScalarType sign=((f->V(1)->T().P()-f->V(0)->T().P()) ^ (f->V(2)->T().P()-f->V(0)->T().P()));
    if (sign!=0) {
      if (sign*lastsign<0) return false;
      lastsign=sign;
    }
  }
  return true;
}

// detects and marks folded faces, by setting their quality to 0 (or 1 otherwise)
// returns number of folded faces
template<class MESH_TYPE>
int MarkTexCoordFolds(MESH_TYPE &m){
  
  assert(HasPerVertexTexCoord(m));
  assert(m.HasPerFaceQuality());
  
  typedef typename MESH_TYPE::VertexType::TexCoordType::PointType::ScalarType ScalarType;
  
  SimpleTempData<typename MESH_TYPE::FaceContainer, short> sign(m.face);
  //sign.Start(0);
  
  // first pass, determine predominant sign
  int npos=0, nneg=0;
  ScalarType lastsign=0;
  for (typename MESH_TYPE::FaceIterator f=m.face.begin(); f!=m.face.end(); f++){
    ScalarType fsign=((f->V(1)->T().P()-f->V(0)->T().P()) ^ (f->V(2)->T().P()-f->V(0)->T().P()));
    if (fsign<0) { sign[f]=-1;  nneg++; }
    if (fsign>0) { sign[f]=+1; npos++; }
  }
  
  // second pass, detect folded faces
  int res=0;
  short gsign= (nneg>npos)?-1:+1;
  for (typename MESH_TYPE::FaceIterator f=m.face.begin(); f!=m.face.end(); f++){
    if (sign[f]*gsign<0){
      res++;
      f->Q()=0;
    } else f->Q()=1;
  }
  
  //sign.Stop();
  
  return res;
}

// Smooths texture coords.
// (can be useful to remove folds, 
//  e.g. these created when obtaining tecture coordinates after projections)
template<class MESH_TYPE>
void SmoothTexCoords(MESH_TYPE &m){
  
  assert(HasPerVertexTexCoord(m));
  
  typedef typename MESH_TYPE::VertexType::TexCoordType::PointType PointType;
  typedef typename MESH_TYPE::VertexType::TexCoordType::ScalarType ScalarType;

  SimpleTempData<typename MESH_TYPE::VertContainer, int> div(m.vert);
  SimpleTempData<typename MESH_TYPE::VertContainer, PointType > sum(m.vert);
  
 /* div.Start();
  sum.Start();*/
  
	for (typename MESH_TYPE::VertexIterator v=m.vert.begin(); v!=m.vert.end(); v++) {
		sum[v].SetZero();
    div[v]=0;
	}

	for (typename MESH_TYPE::FaceIterator f=m.face.begin(); f!=m.face.end(); f++){
		div[f->V(0)] +=2; sum[f->V(0)] += f->V(2)->T().P(); sum[f->V(0)] += f->V(1)->T().P();
		div[f->V(1)] +=2; sum[f->V(1)] += f->V(0)->T().P(); sum[f->V(1)] += f->V(2)->T().P();
		div[f->V(2)] +=2; sum[f->V(2)] += f->V(1)->T().P(); sum[f->V(2)] += f->V(0)->T().P();
	}

	for (typename MESH_TYPE::VertexIterator v=m.vert.begin(); v!=m.vert.end(); v++) 
	if (!v->IsB()) 
  {
		//if (v->div>0) {
	    if (div[v]>0) {
			v->T().P() = sum[v]/((ScalarType)div[v]);
		}
	}
	
	/*div.Stop();
  sum.Stop();*/

}
// MIPSTexCoordFoldHealer
// ----------------------
// Uses MIPS optimization to attempt to remove folds.
// Acts only in proximity of foleded faces!
// Use "iterateUntilConvergence" to unfold faces (returns number of performed iterations, as usual)
// Use "maxStarSize" (direct access) to determine size of affected patch around folded face

// AUTO_SPEED metaparameter:
#define AUTO_SPEED 1
// if set to one, speed is reduced/increase automatically to avoid oscillating behaviour
// (consumes memory and CPU, but increase robustness with speed parameter and sometimes converge faster)

template<class MESH_TYPE> 
class MIPSTexCoordFoldHealer:public MIPSTexCoordOptimization<MESH_TYPE>{
public:
  
  int maxStarSize; // max star size that is affected around a folded face.. Defualt: 3
  
  typedef MESH_TYPE MeshType;
  typedef typename MESH_TYPE::VertexIterator VertexIterator;
  typedef typename MESH_TYPE::FaceIterator FaceIterator;
  typedef typename MESH_TYPE::VertexType VertexType;
  typedef typename MESH_TYPE::FaceType FaceType;
  typedef typename MESH_TYPE::ScalarType ScalarType;

  typedef MIPSTexCoordOptimization<MESH_TYPE> Super; // superclass (commodity)
protected:
  
  SimpleTempData<typename MESH_TYPE::FaceContainer, bool > foldf;
  SimpleTempData<typename MESH_TYPE::VertContainer, bool > foldv;
  
#if AUTO_SPEED
  SimpleTempData<typename MESH_TYPE::VertContainer, Point2<ScalarType> > lastDir;
  SimpleTempData<typename MESH_TYPE::VertContainer, ScalarType > lastSpeed;
#endif
  
  ScalarType sign;
  int nfolds;
  FaceType* aFoldedFace;

  void PropagateFoldF(){
    for (typename MeshType::FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++){
      if (foldv[ f->V(0)] || foldv[ f->V(1)] || foldv[ f->V(2) ] ) {
        foldf[f] = true;
      }
    }    
  }
    
  void PropagateFoldV(){
    for (typename MeshType::FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++) {
      if (foldf[f] ) {
        foldv[ f->V(0) ] = foldv[ f->V(1) ] = foldv[ f->V(2) ] = true;
      }
    }
  }
  
  bool FindFolds(){
  
    /*ScalarType lastsign=0;*/
    int npos=0, nneg=0;
    for (typename MESH_TYPE::FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++){
      ScalarType sign=((f->V(1)->T().P()-f->V(0)->T().P()) ^ (f->V(2)->T().P()-f->V(0)->T().P()));
      if (sign>0) { npos++; }
      if (sign<0) { nneg++; }
    }
    if (npos*nneg==0)     {sign=0; nfolds=0;} else
    if (npos>nneg) { sign=+1; nfolds=nneg; } else
    { sign=-1; nfolds=npos; };
    
    for (typename MeshType::FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++){
      ScalarType signf=((f->V(1)->T().P()-f->V(0)->T().P()) ^ (f->V(2)->T().P()-f->V(0)->T().P()));
	 /* if ((!Super::isFixed[f->V(0)])&&(!Super::isFixed[f->V(1)])&&(!Super::isFixed[f->V(2)]))*/
		foldf[f] = (signf*sign<0);
	 /* else
		foldf[f] =  false;*/
    }
    
    return true;
  }

public:
  
 int IterateUntilConvergence(ScalarType threshold=0.0001, int maxite=50){
	(void)threshold;
	  for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) foldv[v]=false;
    FindFolds();
    PropagateFoldV();
    PropagateFoldF();
   /* int i=0;*/
    int nite = 0, totIte=0, pass=0;
    while (Iterate()>0) { 
      totIte++; 
	  nite++;
      if (nite>=maxite) {
        PropagateFoldV();
        PropagateFoldF();
        nite=0;
        if (pass++>=maxStarSize) break; // number of passes
      } 
    }
    return totIte;
  }

   
  // constructor and destructor
  MIPSTexCoordFoldHealer(MeshType &_m):MIPSTexCoordOptimization<MeshType>(_m),foldf(_m.face),foldv(_m.vert)
#if AUTO_SPEED
                   ,lastDir(_m.vert),lastSpeed(_m.vert,1.0)
#endif
  {sign=0; nfolds=0;  maxStarSize=3; };
  
  ~MIPSTexCoordFoldHealer(){
   /* data.Stop();
    sum.Stop();
    Super::isFixed.Stop();*/
  }
  
  ScalarType Iterate(){
        
    #define v0 (f->V(0)->T().P())
    #define v1 (f->V(1)->T().P())
    #define v2 (f->V(2)->T().P())
    #define vi (f->V(i)->T().P())
    #define vj (f->V(j)->T().P())
    #define vk (f->V(k)->T().P())
	  for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) {
		  //sum[v].Zero();
		  Super::sum[v]=Point2<ScalarType>(0,0);
	  }

    ScalarType totProjArea=0;
    nfolds=0;
	  for (FaceIterator f=Super::m.face.begin(); f!=Super::m.face.end(); f++) {
      if (Super::isFixed[f->V(0)] && Super::isFixed[f->V(1)] && Super::isFixed[f->V(2)]) continue;
      if (!foldf[f]) continue;
      ScalarType area2 = ((v1-v0) ^ (v2-v0));
      if (area2*sign<0) nfolds++;
      totProjArea+=area2;
      ScalarType o[3] = { // (opposite edge)^2 
        (   v1-v2).SquaredNorm(),
        (v0   -v2).SquaredNorm(),
        (v0-v1   ).SquaredNorm(),
      };
      ScalarType e =( Super::data[f][0] * o[0] + 
                      Super::data[f][1] * o[1] + 
                      Super::data[f][2] * o[2] ) / (area2*area2);            

		  for (int i=0; i<3; i++){
        int j=(i+1)%3, k=(i+2)%3;                     
			  ScalarType p=(vj-vi)*(vk-vi);							  
				ScalarType gy= e*(o[k]-p) - 2*Super::data[f][j]; 
				ScalarType gx= e*(o[j]-p) - 2*Super::data[f][k];
				      
				// speed free mode: (a try!)
			  //sum[f->V(i)]+= ( (vj-vi) * gx + (vk-vi) * gy );// / area2; 
			  
			  // speed mode:
        Super::sum[f->V(i)]+= ( (vj-vi) * gx + (vk-vi) * gy ) / area2; 
		  }
	  }
 
    if (nfolds==0) return 0;
    
 	  for (VertexIterator v=Super::m.vert.begin(); v!=Super::m.vert.end(); v++) 
    if (  !Super::isFixed[v] && foldv[v] )
    {
      ScalarType n=Super::sum[v].Norm(); if ( n > 1 ) { Super::sum[v]/=n; n=1.0;}
#if AUTO_SPEED
      if (Super::sum[v]*lastDir[v] < 0.0) lastSpeed[v]*=(ScalarType)0.8; else lastSpeed[v]*=(ScalarType)1.1;
      lastDir[v]=Super::sum[v];
		  v->T().P()-=(Super::sum[v] ) * (Super::speed * lastSpeed[v] );
#else
		  v->T().P()-=(Super::sum[v] ) * Super::speed;
#endif
  	}
  	return (ScalarType)nfolds;
  	#undef v0
    #undef v1 
    #undef v2 
  	#undef vi
    #undef vj 
    #undef vk 
  	//printf("rejected %d\n",rejected);
  }
  
  
};


}	}	// End namespace vcg::tri

#endif //  __VCGLIB__TEXTCOOORD_OPTIMIZATION
