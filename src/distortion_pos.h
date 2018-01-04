#ifndef DISTORTION_POS_H
#define DISTORTION_POS_H

#include "mesh.h"
#include "vertex_position.h"
#include "math_utils.h"

enum ParameterizationGeometry { Model, Texture };

class DistortionMetric {

public:

    enum Type { Area, Angle };

    static void ComputeAreaScale(Mesh& m, double& areaScale, ParameterizationGeometry geometry)
    {
        double sumArea3D = 0;
        double sumAreaUV = 0;
        for (const auto& f : m.face) {
            sumArea3D = Area3D(m, f, geometry);
            sumAreaUV = AreaUV(f);
        }
        areaScale = sumArea3D / sumAreaUV;
    }

    static double AreaUV(const Mesh::FaceType& f)
    {
        Point2d u0 = f.cWT(0).P();
        Point2d u1 = f.cWT(1).P();
        Point2d u2 = f.cWT(2).P();
        return ((u1 - u0) ^ (u2 - u0)) / 2.0;
    }

    static double AngleUV(const Mesh::FaceType& f, int i)
    {
        Point2d u0 = f.cWT(i).P();
        Point2d u1 = f.cWT((i+1)%3).P();
        Point2d u2 = f.cWT((i+2)%3).P();
        return VecAngle(u1 - u0, u2 - u0);
    }

    static double Angle3D(Mesh& m, const Mesh::FaceType& f, int i, ParameterizationGeometry geometry)
    {
        int j = (i+1)%3;
        int k = (i+2)%3;
        if (geometry == Model) {
            return VecAngle(f.cP(j) - f.cP(i), f.cP(k) - f.cP(i));
        } else {
            WedgeTexCoordAttributePosition<Mesh> vpos{m, "WedgeTexCoordStorage"};
            return VecAngle<Point3d>(vpos(&f, j) - vpos(&f, i), vpos(&f, k) - vpos(&f, i));
        }
    }

    static double Area3D(Mesh& m, const Mesh::FaceType& f, ParameterizationGeometry geometry)
    {
        if (geometry == Model) {
            return ((f.cP(1) - f.cP(0)) ^ (f.cP(2) - f.cP(0))).Norm() / 2.0;
        } else {
            WedgeTexCoordAttributePosition<Mesh> vpos{m, "WedgeTexCoordStorage"};
            return ((vpos(&f, 1) - vpos(&f, 0)) ^ (vpos(&f, 2) - vpos(&f, 0))).Norm() / 2.0;
        }
    }

    static double AreaDistortion(Mesh& m, const Mesh::FaceType& f, double areaScale, ParameterizationGeometry geometry)
    {
        double parameterArea = AreaUV(f) * areaScale;
        double faceArea = Area3D(m, f, geometry);
        assert(std::isfinite(faceArea));
        return (parameterArea - faceArea) / faceArea;
    }

    static double AngleDistortion(Mesh& m, const Mesh::FaceType& f, ParameterizationGeometry geometry)
    {
        double d = 0;
        for (int i = 0; i < 3; ++i) {
            double parameterAngle = AngleUV(f, i);
            double faceAngle = Angle3D(m, f, i, geometry);
            d += std::abs(parameterAngle - faceAngle);
        }
        return d;
    }
};






#if 0

#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>
#include <vcg/complex/algorithms/parametrization/tangent_field_operators.h>


/*
 *  Energy types:
 *
 *    AreaDist : 0 for equiareal (equipotent) mappings
 *    EdgeDist (hack): 0 for isometric mappings (computed on edges only)
 *    AngleDist (hack): 0 for conformal mappings
 *    CrossDist : as above, but computed on tangent directions (not UVs)
 *    L2Stretch : 1 for isometric mappings (averaged case on the mesh),
 *                +inf on degenerate / folded cases
 *                Described in [1]
 *    LInfStretch : as above, but WORST case
 *                  (returns the worst stretch on any position and direction)
 *                  Described in [1]
 *    ARAPEnergy : 0 for isometric mappings
 *                 Described in [2]
 *
 * [1] Sander, P. V., Snyder, J., Gortler, S. J., & Hoppe, H.
 *     "Texture mapping progressive meshes."
 *      In Proc. ACM SIGGRAPH (pp. 409-416). 2001
 *
 * [2] Liu, L., Zhang, L., Xu, Y., Gotsman, C., & Gortler, S. J. (2008, July).
 *     A local/global approach to mesh parameterization.
 *     Computer Graphics Forum (Vol. 27, No. 5, pp. 1495-1504). Blackwell Publishing Ltd.
 */

template <class MeshType>
class DistortionPos
{
public:
    typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::VertexType VertexType;
    typedef typename MeshType::CoordType CoordType;
    typedef typename MeshType::ScalarType ScalarType;
    typedef typename MeshType::FaceType::CurVecType CurVecType;
    typedef UVHelper<MeshType, true>          UV;
    typedef typename UV::TexCoordType                 TexCoordType;
	typedef typename TexCoordType::ScalarType         TexScalarType;

	static TexCoordType UVCoord(const FaceType *f, int i)
	{
		return UV::Coord(f, i);
	}

    template <typename VP>
    static ScalarType Area3D(const FaceType *f, const VP& vp)
    {
        return 0.5f * ((vp(f, 1) - vp(f, 0)) ^ vp(f, 2) - vp(f, 0)).Norm();
    }

    static ScalarType AreaUV(const FaceType *f)
    {
		TexCoordType uv0 = UVCoord(f, 0);
		TexCoordType uv1 = UVCoord(f, 1);
		TexCoordType uv2 = UVCoord(f, 2);
        ScalarType AreaUV=((uv1-uv0)^(uv2-uv0))/2.0;
        return AreaUV;
    }

    template <typename VP>
    static ScalarType EdgeLenght3D(const FaceType *f,int e, const VP& vp)
    {
        assert((e>=0)&&(e<3));
        ScalarType length=(vp(f,(e+1)%3) - vp(f,e)).Norm();
        return (length);
    }

    static ScalarType EdgeLenghtUV(const FaceType *f,int e)
    {
        assert((e>=0)&&(e<3));
		TexCoordType uv0 = UVCoord(f, e+0);
		TexCoordType uv1 = UVCoord(f, (e+1)%3);
        ScalarType UVlength=Distance(uv0,uv1);
        return UVlength;
    }

    static ScalarType AngleCos3D(const FaceType *f,int e)
    {
        assert((e>=0)&&(e<3));
        CoordType p0=f->P((e+2)%3);
        CoordType p1=f->P(e);
        CoordType p2=f->P((e+1)%3);
        CoordType dir0=p2-p1;
        CoordType dir1=p0-p1;
        dir0.Normalize();
        dir1.Normalize();
        ScalarType angle=dir0*dir1;
        return angle;
    }

    static ScalarType AngleCosUV(const FaceType *f,int e)
    {
		TexCoordType uv0 = UVCoord(f, (e+2)%3);
		TexCoordType uv1 = UVCoord(f,  e);
		TexCoordType uv2 = UVCoord(f, (e+1)%3);
        vcg::Point2<ScalarType> dir0=uv2-uv1;
        vcg::Point2<ScalarType> dir1=uv0-uv1;
        dir0.Normalize();
        dir1.Normalize();
        ScalarType angle=dir0*dir1;
        return angle;
    }

    template <typename VP>
    static ScalarType AngleRad3D(const FaceType *f,int e, const VP& vp)
    {
        assert((e>=0)&&(e<3));
        CoordType p0=vp(f, (e+2)%3);
        CoordType p1=vp(f, e);
        CoordType p2=vp(f, (e+1)%3);
        CoordType dir0=p2-p1;
        CoordType dir1=p0-p1;
        return Angle(dir0,dir1);
    }

    static ScalarType AngleRadUV(const FaceType *f,int e)
    {
		TexCoordType uv0 = UVCoord(f, (e+2)%3);
		TexCoordType uv1 = UVCoord(f,  e);
		TexCoordType uv2 = UVCoord(f, (e+1)%3);
        vcg::Point2<TexScalarType> dir0=uv2-uv1;
        vcg::Point2<TexScalarType> dir1=uv0-uv1;
        dir0.Normalize();
        dir1.Normalize();
        ScalarType t=dir0*dir1;
        if(t>1) t = 1;
        else if(t<-1) t = -1;
        return acos(t);
    }


public:
    enum DistType{AreaDist,EdgeDist,AngleDist,CrossDist,L2Stretch,LInfStretch,ARAPDist};

    ///return the absolute difference between angle in 3D space and texture space
    ///Actually the difference in cos space
    static ScalarType AngleCosDistortion(const FaceType *f,int e)
    {
        ScalarType Angle_3D=AngleCos3D(f,e);
        ScalarType Angle_UV=AngleCosUV(f,e);
        ScalarType diff=fabs(Angle_3D-Angle_UV);///Angle_3D;
        return diff;
    }
    ///return the absolute difference between angle in 3D space and texture space
    ///Actually the difference in cos space
    template <typename VP>
    static ScalarType AngleRadDistortion(const FaceType *f,int e, const VP& vp)
    {
        ScalarType Angle_3D=AngleRad3D(f,e, vp);
        ScalarType Angle_UV=AngleRadUV(f,e);
        ScalarType diff=fabs(Angle_3D-Angle_UV)/Angle_3D;///Angle_3D;
        return diff;
    }

    ///return the variance of angle, normalized
    ///in absolute value
    template <typename VP>
    static ScalarType AngleDistortion(const FaceType *f, const VP& vp)
    {
        return  (AngleRadDistortion(f,0, vp) +
                AngleRadDistortion(f,1, vp) +
                AngleRadDistortion(f,2, vp))/3.0;
    }

    ///return the global scaling factors  from 3D to UV
    template <typename VP>
    static void MeshScalingFactor(const MeshType &m,
                                  ScalarType &AreaScale,
                                  ScalarType &EdgeScale, const VP& vp)
    {
        ScalarType SumArea3D=0;
        ScalarType SumArea2D=0;
        ScalarType SumEdge3D=0;
        ScalarType SumEdge2D=0;
        for (size_t i=0;i<m.face.size();i++)
        {
            SumArea3D+=Area3D(&m.face[i], vp);
            SumArea2D+=AreaUV(&m.face[i]);
            for (int j=0;j<3;j++)
            {
                SumEdge3D+=EdgeLenght3D(&m.face[i],j, vp);
                SumEdge2D+=EdgeLenghtUV(&m.face[i],j);
            }
        }
        AreaScale=SumArea3D/SumArea2D;
        EdgeScale=SumEdge3D/SumEdge2D;
    }

    ///return the variance of edge length, normalized in absolute value,
    ///the needed scaling factor EdgeScaleVal may be calculated
    ///by using the ScalingFactor function
    template <typename VP>
    static ScalarType EdgeDistortion(const FaceType *f,int e,
                                     ScalarType EdgeScaleVal, const VP& vp)
    {
        ScalarType edgeUV=EdgeLenghtUV(f,e)*EdgeScaleVal;
        ScalarType edge3D=EdgeLenght3D(f,e, vp);
        assert(edge3D > 0);
        ScalarType diff=fabs(edge3D-edgeUV)/edge3D;
        assert(!math::IsNAN(diff));
        return diff;
    }

    ///return the variance of area, normalized
    ///in absolute value, the scalar AreaScaleVal may be calculated
    ///by using the ScalingFactor function
    static ScalarType AreaDistortion(const FaceType *f,
                                     ScalarType AreaScaleVal)
    {
        ScalarType areaUV=AreaUV(f)*AreaScaleVal;
        ScalarType area3D=Area3D(f);
        assert(area3D > 0);
        ScalarType diff=fabs(areaUV-area3D)/area3D;
        assert(!math::IsNAN(diff));
        return diff;
    }

    static ScalarType L2StretchEnergySquared(const FaceType *f,
                                             ScalarType AreaScaleVal)
    {
		TexCoordType p0 = UVCoord(f, 0);
		TexCoordType p1 = UVCoord(f, 1);
		TexCoordType p2 = UVCoord(f, 2);

        CoordType q0 = f->cP(0);
        CoordType q1 = f->cP(1);
        CoordType q2 = f->cP(2);

        TexScalarType A2 = ((p1-p0)^(p2-p0));

        if (A2<0) A2 = 0; // will be NAN, +infinity

        CoordType Ss = ( q0 * ( p1[1]-p2[1] ) + q1 * (p2[1]-p0[1]) + q2 * (p0[1]-p1[1]) ) / A2;
        CoordType St = ( q0 * ( p2[0]-p1[0] ) + q1 * (p0[0]-p2[0]) + q2 * (p1[0]-p0[0]) ) / A2;

        ScalarType a = Ss.SquaredNorm() / AreaScaleVal;
        ScalarType c = St.SquaredNorm() / AreaScaleVal;

        return ((a+c)/2);
    }



    static ScalarType LInfStretchEnergy(const FaceType *f,  ScalarType AreaScaleVal)
    {
		TexCoordType p0 = UVCoord(f, 0);
		TexCoordType p1 = UVCoord(f, 1);
		TexCoordType p2 = UVCoord(f, 2);

        CoordType q0 = f->cP(0);
        CoordType q1 = f->cP(1);
        CoordType q2 = f->cP(2);

        TexScalarType A2 = ((p1-p0)^(p2-p0));

        if (A2<0) A2 = 0; // will be NAN, +infinity

        CoordType Ss = ( q0 * ( p1[1]-p2[1] ) + q1 * (p2[1]-p0[1]) + q2 * (p0[1]-p1[1]) ) / A2;
        CoordType St = ( q0 * ( p2[0]-p1[0] ) + q1 * (p0[0]-p2[0]) + q2 * (p1[0]-p0[0]) ) / A2;

        ScalarType a = Ss.SquaredNorm() / AreaScaleVal;
        ScalarType b = Ss*St / AreaScaleVal;
        ScalarType c = St.SquaredNorm() / AreaScaleVal;

        ScalarType delta = sqrt((a-c)*(a-c)+4*b*b);
        ScalarType G =  sqrt( (a+c+delta)/2 );
        //ScalarType g = sqrt( (a+c-delta)/2 ); // not needed
        return G;
    }

	static ScalarType ARAPEnergy(const FaceType *f)
	{
		if (f == NULL)
		{
			return std::numeric_limits<ScalarType>::infinity();
		}

		const Eigen::Matrix2d F = mappingTransform2D(*f);
		const Eigen::Vector2d singular = svd2x2(F);
		const double a = singular(0) - 1;
		const double b = singular(1) - 1;
		return ScalarType(0.5 * (a*a + b*b));
	}

	static Eigen::Matrix2d mappingTransform2D(const FaceType & triangle)
	{
		typedef  Eigen::Matrix<double, 3, 2> Matrix32;
		typedef  Eigen::Matrix2d             Matrix22;


		Matrix22 param3d, param2d;
		// 3D
		{
			Matrix32 edges3D, P3D;
			Eigen::Vector3d e0, e1;
			(triangle.cP(1) - triangle.cP(0)).ToEigenVector(e0); // 0->1
			(triangle.cP(2) - triangle.cP(0)).ToEigenVector(e1); // 0->2
			edges3D.col(0) = e0;
			edges3D.col(1) = e1;

			// Projection/frame change matrix
			P3D.col(0) = edges3D.col(0).normalized(); // 0->1 normalized                                 e0 basis
			P3D.col(1) = (edges3D.col(1) - edges3D.col(1).dot(P3D.col(0)) * P3D.col(0)).normalized(); // e1 basis orthogonal to e0

			param3d = (P3D.transpose() * edges3D);
		}

		// 2D
		{
			Matrix22 edges2D, P2D;
			TexCoordType uv0 = UVCoord(&triangle, 0);
			TexCoordType uv1 = UVCoord(&triangle, 1);
			TexCoordType uv2 = UVCoord(&triangle, 2);

			const TexCoordType e0 = (uv1 - uv0); // 0->1
			const TexCoordType e1 = (uv2 - uv0); // 0->2
			param2d << e0.X(), e1.X(),
			           e0.Y(), e1.Y();
		}

		return param2d * param3d.inverse(); // transf mapping
	}

	// svd 2x2 matrix (singular values only)
	static Eigen::Vector2d svd2x2(const Eigen::Matrix2d & M)
	{
		const double a=M(0,0), b=M(0,1), c=M(1,0), d=M(1,1);
		const double tmp1 = a*a + b*b;
		const double tmp2 = c*c + d*d;
		const double s1 = tmp1 + tmp2;
		const double s2 = std::sqrt(std::pow((tmp1 -tmp2), 2.0) + 4 * std::pow(a*c + b*d, 2.0));
		return Eigen::Vector2d(std::sqrt((s1+s2)/2.0), std::sqrt((s1-s2)/2.0));
	}


    ///return the number of folded faces
    static bool Folded(const FaceType *f)
    {
        ScalarType areaUV=AreaUV(f);
        /*if (areaUV<0)
                    printf("area %5.5f \n",areaUV);*/
        return (areaUV<0);
    }

    static int Folded(const MeshType &m)
    {
        int folded=0;
        for (size_t i=0;i<m.face.size();i++)
        {
            if (m.face[i].IsD())continue;
            if(Folded(&m.face[i]))folded++;
        }
        return folded;
    }

    static bool GloballyUnFolded(const MeshType &m)
    {
        int num=Folded(m);
        return (num>(m.fn)/2);
    }

    static ScalarType MeshAngleDistortion(const MeshType &m)
    {
        ScalarType UDdist=0;
        for (size_t i=0;i<m.face.size();i++)
        {
            if (m.face[i].IsD())continue;
            const FaceType *f=&(m.face[i]);
            UDdist+=AngleDistortion(f)*Area3D(f);
        }
        return UDdist;
    }

    static ScalarType SetFQAsCrossDirDistortion(MeshType &m)
    {
        //first save the old UV dir
        std::vector<CurVecType> Dir1,Dir2;
        for (size_t i=0;i<m.face.size();i++)
        {
            Dir1.push_back(m.face[i].PD1());
            Dir2.push_back(m.face[i].PD2());
        }
        vcg::tri::CrossField<MeshType>::InitDirFromWEdgeUV(m);

        ScalarType tot = 0, totA = 0;

        //then compute angle deficit
        for (size_t i=0;i<m.face.size();i++)
        {

            FaceType &f( m.face[i] );
            CoordType transfPD1=vcg::tri::CrossField<MeshType>::K_PI(CoordType::Construct( Dir1[i] ),
                                                                     CoordType::Construct( f.PD1() ),
                                                                     f.N());
            transfPD1.Normalize();
            ScalarType AngleDeficit=vcg::Angle(transfPD1,CoordType::Construct( f.PD1() ));
            AngleDeficit=math::ToDeg(AngleDeficit);
            if ((AngleDeficit>45)||(AngleDeficit<0))
            {
                std::cout<<"Warnign A Deficit "<<AngleDeficit<<std::endl;
            }
//            assert(AngleDeficit<45);
//            assert(AngleDeficit>=0);

            ScalarType doubleArea = vcg::DoubleArea( f );
            ScalarType distortion = (AngleDeficit)/ 45 ;

            m.face[i].Q()= distortion;
            tot += distortion * doubleArea;
            totA += doubleArea;
        }

        //finally restore the original directions
        for (size_t i=0;i<m.face.size();i++)
        {
            m.face[i].PD1()=Dir1[i];
            m.face[i].PD2()=Dir2[i];
        }

        return tot / totA;
    }

    static ScalarType SetQasDistorsion(MeshType &m, DistType DType=AreaDist)
    {
        if (DType==CrossDist)
        {
            ScalarType res = SetFQAsCrossDirDistortion(m);

            vcg::tri::UpdateQuality<MeshType>::VertexFromFace(m,true);
            return res;
        }

        ScalarType edge_scale,area_scale;
        MeshScalingFactor(m,area_scale,edge_scale);

        ScalarType tot = 0;
        ScalarType totA = 0;

        for (size_t i=0;i<m.face.size();i++)
        {
            if (m.face[i].IsD())continue;
            ScalarType q;
            switch (DType) {
            case CrossDist:
                // make compiler happy
                q = 0;
                break;
            case AreaDist:
                q = AreaDistortion(&m.face[i],area_scale);
                break;
            case AngleDist:
                q = AngleDistortion(&m.face[i]);
                break;
            case EdgeDist:
                q =( EdgeDistortion(&m.face[i],0,edge_scale)+
                     EdgeDistortion(&m.face[i],1,edge_scale)+
                     EdgeDistortion(&m.face[i],2,edge_scale) )/3;
                break;
            case L2Stretch:
                q = L2StretchEnergySquared( &m.face[i],area_scale );
                break;
            case LInfStretch:
                q = LInfStretchEnergy( &m.face[i],area_scale );
                break;
            case ARAPDist:
                q = ARAPEnergy(&m.face[i]);
                break;
            }

            m.face[i].Q() = q; // note: for L2Stretch, we are puttning E^2 on Q

            // aggregate:
            if (DType==LInfStretch) {
                tot = std::max( tot, q );
            } else {
                ScalarType a = Area3D(&m.face[i]);
                tot += q*a;
                totA += a;
            }

        }

        vcg::tri::UpdateQuality<MeshType>::VertexFromFace(m,true);

        switch (DType) {
            case L2Stretch: return sqrt(tot/totA);
            case LInfStretch: return tot;
            default:  return tot/totA;
        }
    }
};
#endif


#endif
