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

#include "Eigen\Dense"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#ifndef M_PI
#define M_PI 3.141592653589793238
#endif

#ifndef RAN_PRECISION
#define RAN_PRECISION 9999
#endif

using Eigen::MatrixXd;
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
		void delete_face_2d(F * f);
		void edge_swap(E * e);
		void face_split(F * f, V * v);
		double triangle_space(V * v1, V * v2, V * v3);
		bool barycentric_coordinate(H * f, V * v);
		F * locate_face(V * v);
		bool outside_circle(V * a, V * b, V * c, V * v);
		bool legalize_edge(V * v, E * e);
		void insert_vertex(V * v);
		double get_random();
		void init_mesh();
		void generate_mesh(int num);
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

	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::delete_face_2d(F * f)
	{
		std::map<int, tFace>::iterator fiter = m_map_face.find(f->id());
		if (fiter != m_map_face.end())
		{
			m_map_face.erase(fiter);
		}
		m_faces.remove(f);

		H * pH = (H *)f->halfedge();
		E * pE = (E *)pH->edge();

		for (int i = 0; i < 3; ++i) //set edge's halfedge to null
		{
			if (pE->halfedge(0) == pH)
			{
				if (pE->halfedge(1) == NULL) //if both halfedges are set to null, in other words , the edge is on boundary, delete the edge
				{
					V * tmpV1 = (V *)pE->halfedge(0)->source();
					V * tmpV2 = (V *)pE->halfedge(0)->target();
					std::list<CEdge *>list = (tmpV1->id() > tmpV2->id()) ? tmpV1->edges():tmpV2->edges();
					std::cout << list.size();

					list.remove(pE);
					std::cout << list.size();
					m_edges.remove(pE);
					delete pE->halfedge(0);
					pE->halfedge(0) = NULL;
					delete pE;
					pE = NULL;
				}
				else
				{
					pE->halfedge(0) = pE->halfedge(1);
					pE->halfedge(1) = NULL;
				}

			}
			else
			{
				pE->halfedge(1) = NULL;
			}
			pH = (H *)pH->he_next();
			pE = (E *)pH->edge();
		}

		delete f;
	}
	/*
	swap edge cwly
	*/
	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::edge_swap(E * e)
	{
		//get two halfedges of the edge
		H * h1 = (H *)e->halfedge(0);
		H * h2 = (H *)e->halfedge(1);
		//get two vertexes of the edge
		V * v1 = (V *)h1->source();
		V * v2 = (V *)h1->target();
		V * v3 = (V *)h1->he_next()->target();
		V * v4 = (V *)h2->he_next()->target();
		//changes on faces
		F * f1 = (F *)h1->face();
		if (faceHalfedge(f1) == h1->he_prev()) //if halfedge on the face needs changing
		{
			f1->halfedge() = h1; //simply set to h1
		}
		F * f2 = (F *)h2->face();
		if (faceHalfedge(f2) == h1->he_prev()) //if halfedge on the face needs changing
		{
			f2->halfedge() = h2; //simply set to h2
		}
		//changes on halfedges
		//target vertex on halfedges
		h1->vertex() = h2->he_next()->vertex();
		h2->vertex() = h1->he_next()->vertex(); //note that next pointer has not been changed yet
		//next & prev
		H * h1_next = (H *)h1->he_next();
		H * h1_prev = (H *)h1->he_prev();
		H * h2_next = (H *)h2->he_next();
		H * h2_prev = (H *)h2->he_prev();

		h1->he_next() = h2_prev;
		h1->he_prev() = h1_next;
		h2->he_next() = h1_prev;
		h2->he_prev() = h2_next;

		h1_next->he_next() = h1;
		h1_next->he_prev() = h2_prev;
		h1_prev->he_next() = h2_next;
		h1_prev->he_prev() = h2;
		
		h2_next->he_next() = h2;
		h2_next->he_prev() = h1_prev;
		h2_prev->he_next() = h1_next;
		h2_prev->he_prev() = h1;
		//faces
		h1_prev->face() = h2->face();
		h2_prev->face() = h1->face();

		//changes on vertexes
		//m_edges
		std::list<CEdge *> * tmpList = (v1->id() < v2->id()) ? v2->edges() : v1->edges();
		tmpList->remove(h1->edge());

		tmpList = (v3->id() < v4->id()) ? v4->edges() : v3->edges();
		tmpList->push_back(h1->edge());

		//halfedge
		if ((H *)v2->halfedge() == h1) v2->halfedge() = h2_prev;
		if ((H *)v1->halfedge() == h2) v1->halfedge() = h1_prev;
	}

	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::face_split(F * f, V * v)
	{
		V * fvarr[3];//get the three vertexes on face f
		fvarr[0] = (V *)f->halfedge()->source();
		fvarr[1] = (V *)f->halfedge()->target();
		fvarr[2] = (V *)f->halfedge()->he_next()->target();

		delete_face_2d(f); //delete old face

		for (int i = 0; i < 3; ++i) //construct three new faces
		{
			V * varr[3] = { fvarr[i], fvarr[(i + 1) % 3], v };
			createFace(varr, numFaces() + 1);
		}
	}

	template<typename V, typename E, typename F, typename H>
	double MeshLib::MyMesh<V, E, F, H>::triangle_space(V * a, V * b, V * c)
	{
		MatrixXd mat(3, 3);
		mat << a->point()[0], a->point()[1], 1,
			b->point()[0], b->point()[1], 1,
			c->point()[0], c->point()[1], 1;
		return 0.5 * mat.determinant();
	}

	template<typename V, typename E, typename F, typename H>
	bool MyMesh<V, E, F, H>::barycentric_coordinate(H * h, V * v)
	{
		return ((triangle_space(v, (V *)h->source(), (V *)h->target()) / triangle_space((V *)h->source(), (V *)h->target(), (V *) h->he_next()->target())) > 0);
	}

	template<typename V, typename E, typename F, typename H>
	F * MyMesh<V, E, F, H>::locate_face(V * v)
	{
		//Arbitrarily choose one face
		F * f = (F *)*faces().begin();
		H * h = (H *)f->halfedge();
		for (;;)
		{
			int i;
			for (i = 0; i < 3; ++i) //traverse three halfedges
			{
				if (!barycentric_coordinate(h, v)) //if bartcentric coordinate is nagative
				{
					f = (F *)h->he_sym()->face(); //get the face adjacent to the current face sharing h
					if (NULL == f)
						return NULL;
					h = (H *)f->halfedge();
					break;
				}
				else
				{
					h = (H *)h->he_next();
				}
			}
			if (i == 3) //if all edges are positive
				return f;
		}
	}
	//if v is outside the circle through a, b, c
	template<typename V, typename E, typename F, typename H>
	bool MyMesh<V, E, F, H>::outside_circle(V * a, V * b, V * c, V * v)
	{
		MatrixXd mat(4, 4);
		mat << a->point()[0], a->point()[1], a->point()[0] * a->point()[0] + a->point()[1] * a->point()[1], 1,
			b->point()[0], b->point()[1], b->point()[0] * b->point()[0] + b->point()[1] * b->point()[1], 1,
			c->point()[0], c->point()[1], c->point()[0] * c->point()[0] + c->point()[1] * c->point()[1], 1,
			v->point()[0], v->point()[1], v->point()[0] * v->point()[0] + v->point()[1] * v->point()[1], 1;
		return (mat.determinant() < 0);
	}

	template<typename V, typename E, typename F, typename H>
	bool MyMesh<V, E, F, H>::legalize_edge(V * v, E * e)
	{
		if (e->boundary())
			return true;
		//get the two vertex of the edge
		V * v0 = (V *)e->halfedge(0)->source();
		V * v1 = (V *)e->halfedge(0)->target();
		assert(e->halfedge(1) != NULL);
		V * v2;
		//find v2 oppsite v
		if (e->halfedge(0)->he_next()->target() == v)
		{
			v2 = (V *)e->halfedge(1)->he_next()->target();
		}
		else
		{
			v2 = (V *)e->halfedge(0)->he_next()->target();
		}
		if (outside_circle(v0, v1, v, v2))
		{
			return false;
		}
		legalize_edge(v, vertexEdge(v0, v2));
		legalize_edge(v, vertexEdge(v1, v2));
		return true;
	}

	//insert a vertex
	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::insert_vertex(V * v)
	{
		F * Fp = locate_face(v);
		V * v1 = (V *)Fp->halfedge()->source();
		V * v2 = (V *)Fp->halfedge()->target();
		face_split(Fp, v);
		//because the three original halfedges have been deleted, we get them by vertex
		H * Hp = vertexHalfedge(v1, v2);
		legalize_edge(v, (E *)Hp->edge());
		Hp = (H *)Hp->he_next()->he_sym()->he_next(); //second halfedge of Fp
		legalize_edge(v, (E *)Hp->edge());
		Hp = (H *)Hp->he_next()->he_sym()->he_next(); //third halfedge of Fp
		legalize_edge(v, (E *)Hp->edge());
	}

	/*
	get a random number between -1 and 1
	*/
	template<typename V, typename E, typename F, typename H>
	inline double  MyMesh<V, E, F, H>::get_random()
	{
		int r = rand() % (RAN_PRECISION + 1);
		r -= ((RAN_PRECISION + 1) / 2);

		return r /(0.5 * (double)(RAN_PRECISION + 1));
	}

	/*
	init the first large face 
	*/
	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::init_mesh()
	{
		// first large face  value choose arbitrarily
		CPoint p1, p2, p3;
		p1[0] = 0.103, p1[1] = 2.970;
		p2[0] = -2.402, p2[1] = -1.051;
		p3[0] = 2.511, p3[1] = 1.032;

		V * varr[3];
		varr[0] = createVertex(1);
		varr[0]->point() = p1;
		varr[1] = createVertex(2);
		varr[1]->point() = p2;
		varr[2] = createVertex(3);
		varr[2]->point() = p3;

		createFace(varr, 1);
	}
	/*
	generate a mesh with num vertexes
	*/
	template<typename V, typename E, typename F, typename H>
	void MyMesh<V, E, F, H>::generate_mesh(int num)
	{
		int id = 4;
		for (int i = 0; i < num; ++i)
		{
			CPoint p;
			for (int j = 0; j < 2; ++j)
			{
				p[j] = get_random();
			}
			p[2] = 0;
			id = i + 4;
			V * pV = createVertex(id);
			pV->point() = p;

			insert_vertex(pV);
		}
	}

}



#endif // !_MY_MESH_
