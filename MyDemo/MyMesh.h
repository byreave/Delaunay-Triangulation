#ifndef _MY_MESH_
#define _MY_MESH_


#include "Mesh\Vertex.h"
#include "Mesh\Edge.h"
#include "Mesh\Face.h"
#include "Mesh\HalfEdge.h"
#include "Mesh\BaseMesh.h"

#include "Mesh\boundary.h"
#include "Mesh\iterators.h"
#include "Parser\parser.h"

#include <math.h>
#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

namespace MeshLib
{
	class CMyVertex;
	class CMyEdge;
	class CMyFace;
	class CMyHalfEdge;

	
	class CMyVertex : public CVertex
	{
	public:
		CMyVertex() : m_rgb(1,1,1) {};
		~CMyVertex() {};

		void _from_string() ;
		void _to_string();

		CPoint & rgb() { return m_rgb; };
	protected:
		CPoint m_rgb;
	};

	inline void CMyVertex::_from_string()
	{
		CParser parser(m_string);
		for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			CToken * token = *iter;
			if (token->m_key == "uv") //CPoint2
			{
				token->m_value >> m_uv;
			}
			if (token->m_key == "rgb") // CPoint
			{
				token->m_value >> m_rgb;
			}
		}
	}

	inline void CMyVertex::_to_string()
	{
		CParser parser(m_string);
		parser._removeToken("uv");

		parser._toString(m_string);
		std::stringstream iss;

		iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";

		if (m_string.length() > 0)
		{
			m_string += " ";
		}
		m_string += iss.str();
	}
	
	class CMyEdge : public CEdge
	{
	public:
		CMyEdge() :m_sharp(false) {};
		~CMyEdge() {};

		void _from_string();

		bool & sharp() { return m_sharp; };
	protected:
		bool m_sharp;
	};

	inline void CMyEdge::_from_string()
	{
		CParser parser(m_string);
		for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
		{
			CToken * token = *iter;
			if (token->m_key == "sharp") // bool
			{
				m_sharp = true;
			}
		}
	}

	class CMyFace : public CFace
	{
	public:

		CPoint & normal() { return m_normal; };
	protected:
		CPoint m_normal;
	};

	class CMyHalfEdge : public CHalfEdge
	{
	public:
		double length();
	};

	template<typename V, typename E, typename F, typename H>
	class MyMesh : public CBaseMesh<V, E, F, H>
	{
	public:
		typedef CBoundary<V, E, F, H>					CBoundary;
		typedef CLoop<V, E, F, H>						CLoop;

		typedef MeshVertexIterator<V, E, F, H>			MeshVertexIterator;
		typedef MeshEdgeIterator<V, E, F, H>			MeshEdgeIterator;
		typedef MeshFaceIterator<V, E, F, H>			MeshFaceIterator;
		typedef MeshHalfEdgeIterator<V, E, F, H>		MeshHalfEdgeIterator;

		typedef VertexVertexIterator<V, E, F, H>		VertexVertexIterator;
		typedef VertexEdgeIterator<V, E, F, H>			VertexEdgeIterator;
		typedef VertexFaceIterator<V, E, F, H>			VertexFaceIterator;
		typedef VertexInHalfedgeIterator<V, E, F, H>	VertexInHalfedgeIterator;
		typedef VertexOutHalfedgeIterator<V, E, F, H>	VertexOutHalfedgeIterator;

		typedef FaceVertexIterator<V, E, F, H>			FaceVertexIterator;
		typedef FaceEdgeIterator<V, E, F, H>			FaceEdgeIterator;
		typedef FaceHalfedgeIterator<V, E, F, H>		FaceHalfedgeIterator;

		void output_mesh_info();
		void test_iterator();
		double cosine_law(double a, double b, double c);
		void verify_gauss_bonnet();
	};

	typedef MyMesh<CMyVertex, CMyEdge, CMyFace, CMyHalfEdge> CMyMesh;

	double MeshLib::CMyHalfEdge::length()
	{
		CVertex *v1 = source();
		CVertex *v2 = target();
		CPoint c = v1->point() - v2->point();
		return c.norm();
	}

	/*
	the law of cosines
	*/
	template<typename V, typename E, typename F, typename H>
	double MeshLib::MyMesh<V, E, F, H>::cosine_law(double a, double b, double c)
	{
		return acos((b*b + c*c - a*a) /(2 * b * c));
	}

	template<typename V, typename E, typename F, typename H>
	void MeshLib::MyMesh<V, E, F, H>::output_mesh_info()
	{
		int nv = this->numVertices();
		int ne = this->numEdges();
		int nf = this->numFaces();

		std::cout << "#V=" << nv << "  ";
		std::cout << "#E=" << ne << "  ";
		std::cout << "#F=" << nf << "  ";

		int euler_char= nv - ne + nf;
		std::cout << "Euler's characteristic=" << euler_char << "  ";

		CBoundary boundary(this);
		std::vector<CLoop*> & loops = boundary.loops();
		int nb = loops.size();

		int genus = (2 - (euler_char + nb)) / 2;
		std::cout << "genus=" << genus << std::endl;
	}

	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::test_iterator()
	{
		for (MeshVertexIterator viter(this); !viter.end(); ++viter)
		{
			double sum_alpha = 0;
			V * pV = *viter;
			// you can do something to the vertex here
			// ...
			H * pE1 = (CMyHalfEdge *)pV->halfedge(); // get first income halfedge
			for (int i = 0; i < 6; ++i)
			{
				H * pE2 = (CMyHalfEdge *)pE1->he_next();
				H * pE3 = (CMyHalfEdge *)pE2->he_next();
				double a = pE1->length();
				double b = pE2->length();
				double c = pE3->length();
				sum_alpha += cosine_law(a, b, c);
				std::cout << a << b << c << std::endl;
				//std::cout << "b*b-c*c" << (b*b + c*c - a*a) / (2 * b * c) << std::endl;
				std::cout << i << ":" << cosine_law(a, b, c) << std::endl;

				pE1 = (CMyHalfEdge *)pE2->he_sym(); //compute next triangle
			}
			std::cout << "sum :" << sum_alpha<< std::endl;
			break;

			for (VertexEdgeIterator veiter(pV); !veiter.end(); ++veiter)
			{
				E * pE = *veiter;
				// you can do something to the neighboring edges with CCW
				// ...
			}

			for (VertexFaceIterator vfiter(pV); !vfiter.end(); ++vfiter)
			{
				F * pF = *vfiter;
				// you can do something to the neighboring faces with CCW
				// ...
			}

			for (VertexInHalfedgeIterator vhiter(this, pV); !vhiter.end(); ++vhiter)
			{
				H * pH = *vhiter;
				// you can do something to the incoming halfedges with CCW
				// ...
			}
		}

		for (MeshEdgeIterator eiter(this); !eiter.end(); ++eiter)
		{
			E * pE = *eiter;
			// you can do something to the edge here
			// ...
		}

		for (MeshFaceIterator fiter(this); !fiter.end(); ++fiter)
		{
			F * pF = *fiter;
			// you can do something to the face here
			// ...
		}

		//there are some other iterators which you can find them in class MyMesh

		std::cout << "Iterators test OK.\n";
		//std::cout << i << std::endl;
	}
	/*
	Verify Gauss-Bonnet Theorem
	*/
	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::verify_gauss_bonnet()
	{
		int nv = this->numVertices();
		int ne = this->numEdges();
		int nf = this->numFaces();
		int euler_char = nv - ne + nf;
		double result = 0;
		double sum_alpha = 0;
		for (MeshVertexIterator viter(this); !viter.end(); ++viter)
		{
			sum_alpha = 0;
			V * pV = *viter;
			if (pV->boundary())
			{
				H * pE1 = (CMyHalfEdge *)pV->halfedge();// get first income halfedge
				while (!pE1->edge()->boundary()) //loop until the halfedge is on boundary for convenient purpose
				{
					pE1 = (CMyHalfEdge *)pE1->he_sym()->he_next()->he_next();
				}
				//H * tmp = pE1;
				for (int i = 0;; ++i)
				{
					H * pE2 = (CMyHalfEdge *)pE1->he_next();
					H * pE3 = (CMyHalfEdge *)pE2->he_next();
					double a = pE1->length();
					double b = pE2->length();
					double c = pE3->length();
					sum_alpha += cosine_law(a, b, c); //compute angles via tha law of cosines

					pE1 = (CMyHalfEdge *)pE2->he_sym(); //compute next triangle
					if (pE1 == NULL)//if the pointer goes to another boundary, end loop 
						break;
				}
				result += M_PI - sum_alpha;
			}
			else
			{
				H * pE1 = (CMyHalfEdge *)pV->halfedge(); // get first income halfedge
				H * tmp = pE1;

				for (int i = 0; ; ++i)
				{
					H * pE2 = (CMyHalfEdge *)pE1->he_next();
					H * pE3 = (CMyHalfEdge *)pE2->he_next();
					double a = pE1->length();
					double b = pE2->length();
					double c = pE3->length();
					sum_alpha += cosine_law(a, b, c);//compute angles via tha law of cosines

					pE1 = (CMyHalfEdge *)pE2->he_sym(); //compute next triangle
					if (pE1 == tmp)//if the pointer goes back to the start, end loop
						break;
				}
				result += 2 * M_PI - sum_alpha;
			}
		}
		std::cout <<"Equation left = "<< result << std::endl;
		std::cout << "2 PI multiply Euler Characteristic = " << 2 * M_PI * euler_char << std::endl;
	}
}



#endif // !_MY_MESH_
