/*
 *   Boolean_operations.h
 *
 *   Author Cyril Leconte and Giles Bathgate
 *
 *   Original work Copyright (C) 2010 Cyril Leconte
 *   Modified work Copyright (C) 2015 Giles Bathgate
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef BOOLEAN_OPERATIONS_H
#define BOOLEAN_OPERATIONS_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel Enriched_kernel;

typedef unsigned long Id;
typedef Id VertexId;
typedef Id HalfedgeId;
typedef Id FacetId;
typedef Id InterId;

template <class Refs, class T>
class Enriched_facet :
	public CGAL::HalfedgeDS_face_base<Refs, T>
{
public:
	bool IsExt;
	bool IsOK;
	FacetId Label;
};

template <class Refs, class Tprev, class Tvertex, class Tface>
class Enriched_halfedge :
	public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
public:
	HalfedgeId Label;
};

// a refined vertex with a label
template <class Refs, class T, class P>
class Enriched_vertex :
	public CGAL::HalfedgeDS_vertex_base<Refs, T, P>
{
public:
	VertexId Label;
	Enriched_vertex() {}
	Enriched_vertex(const P& pt) : CGAL::HalfedgeDS_vertex_base<Refs, T, P>(pt)
	{
		this->point() = pt;
	}

};

struct Enriched_items : public CGAL::Polyhedron_items_3 {
	// wrap vertex
	template <class Refs, class Traits>
	struct Vertex_wrapper {
		typedef typename Traits::Point_3  Point;
		typedef Enriched_vertex<Refs,
				CGAL::Tag_true,
				Point> Vertex;
	};

	// wrap face
	template <class Refs, class Traits>
	struct Face_wrapper {
		typedef Enriched_facet<Refs,
				CGAL::Tag_true> Face;
	};

	// wrap halfedge
	template <class Refs, class Traits>
	struct Halfedge_wrapper {
		typedef Enriched_halfedge<Refs,
				CGAL::Tag_true,
				CGAL::Tag_true,
				CGAL::Tag_true> Halfedge;
	};
};

template <class kernel, class items>
class Enriched_polyhedron :
	public CGAL::Polyhedron_3<kernel,items>
{
public:
	typedef typename Enriched_polyhedron<kernel,items>::Facet_iterator Facet_iterator;
	typedef typename Enriched_polyhedron<kernel,items>::Halfedge_handle Halfedge_handle;

	void triangulate()
	{
		Facet_iterator f = this->facets_begin();
		Facet_iterator f2 = this->facets_begin();
		do { //for (; f != this->facets_end(); f++)
			f = f2;
			if(f == this->facets_end()) {
				break;
			}
			f2++;

			if(!(f->is_triangle())) {
				int num = (int)(f->facet_degree() - 3);
				Halfedge_handle h = f->halfedge();

				h = this->make_hole(h);

				Halfedge_handle g = h->next();
				g = g->next();
				Halfedge_handle new_he = this->add_facet_to_border(h, g);
//				new_he->texture_coordinates(h->texture_coordinates(0),h->texture_coordinates(1));
//				new_he->opposite()->texture_coordinates(g->texture_coordinates(0),g->texture_coordinates(1));
				g=new_he;

				num--;
				while(num != 0) {
					g = g->opposite();
					g = g->next();
					Halfedge_handle new_he = this->add_facet_to_border(h, g);
//					new_he->texture_coordinates(h->texture_coordinates(0),h->texture_coordinates(1));
//					new_he->opposite()->texture_coordinates(g->texture_coordinates(0),g->texture_coordinates(1));
					g=new_he;

					num--;
				}

				this->fill_hole(h);
			}

		} while(true);

		//this->compute_normals();
		//this->compute_type();
	}
};

typedef Enriched_polyhedron<Enriched_kernel, Enriched_items>	MEPP_Polyhedron;

/*!
 * \def BOOLEAN_OPERATIONS_DEBUG
 * \brief Enables computation time measuring
 */

//#define BOOLEAN_OPERATIONS_DEBUG

/*!
 * \enum Bool_Op
 * \brief The three Boolean operations
 */
enum Bool_Op {UNION, INTER, MINUS};


#ifdef BOOLEAN_OPERATIONS_DEBUG

/**
 * \fn inline double tr(double &n)
 * \brief Truncate a number to 1/1000
 *        (only if BOOLEAN_OPERATIONS_DEBUG is enable)
 * \param n : The input number in double
 * \return The truncation in double.
 */
inline double tr(double& n)
{
	return floor(n*1000)/1000;
}
#endif // BOOLEAN_OPERATIONS_DEBUG


#include <CGAL/Polyhedron_incremental_builder_3.h>

/*!
 * \class CPolyhedron_from_polygon_builder_3
 * \brief A polyhedron incremental builder
 */
template <class HDS>
class CPolyhedron_from_polygon_builder_3 : public CGAL::Modifier_base<HDS>
{

public:
	/*!
	 * \typedef typename Point_3
	 * \brief 3d point of an halfedge data structure
	 */
	typedef typename HDS::Traits::Point_3									Point_3;

	/*!
	 * \typedef typename Indices
	 * \brief A list of indices (unsigned long) to describe a facet
	 */
	typedef typename std::vector<unsigned long>								Indices;

	/*!
	 * \typedef typename Builder
	 * \brief The polyhedron incremental builder
	 */
	typedef typename CGAL::Polyhedron_incremental_builder_3<HDS>			Builder;

	typedef typename Builder::Halfedge_handle Halfedge_handle;

	typedef typename Builder::Facet_handle Facet_handle;

private:
	// Member variables
	/*! \brief List of the vertices*/
	std::vector<Point_3>													m_Sorted_vertices;
	/*! \brief List of the facets*/
	std::vector<Indices>													m_Facets_indices;

public:
	// Constructors
	/*!
	 * \brief Constructor
	 */
	CPolyhedron_from_polygon_builder_3() {}

	/*!
	 * \brief Adds a triangular facet from a facet of a polyhedron
	 * \param f : The facet handle
	 * \param invert : must be true if the orientation of the facet must be inverted
	 */
	void add_triangle(Facet_handle& f, bool invert)
	{
		//initially, the label of the vertices is 0xFFFFFFFF. if a vertex is added to the result, the tag is set to
		//the number of vertices added.

		//creation of a list of indices
		Indices	vi;

		//adding the first vertex to the result and adding its label to the list of indices.
		//if the vertex is already added (its label is not 0xFFFFFFFF) we only need to add
		//this label to "vi" without adding the vertex.
		Halfedge_handle he = f->facet_begin();
		if(he->vertex()->Label == 0xFFFFFFFF) add_vertex(he->vertex()->point(), he->vertex()->Label);
		vi.push_back(he->vertex()->Label);

		//the order of the two other vertices depends on the orientation of the facet
		if(!invert) {
			if(he->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->vertex()->point(), he->next()->vertex()->Label);
			vi.push_back(he->next()->vertex()->Label);
			if(he->next()->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->next()->vertex()->point(), he->next()->next()->vertex()->Label);
			vi.push_back(he->next()->next()->vertex()->Label);
		} else {
			if(he->next()->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->next()->vertex()->point(), he->next()->next()->vertex()->Label);
			vi.push_back(he->next()->next()->vertex()->Label);
			if(he->next()->vertex()->Label == 0xFFFFFFFF) add_vertex(he->next()->vertex()->point(), he->next()->vertex()->Label);
			vi.push_back(he->next()->vertex()->Label);
		}

		//finally, "vi" is added to the list of the facets
		m_Facets_indices.push_back(vi);
	}

	/*!
	 * \brief Adds a list of triangular facet from an intersected facet of a polyhedron
	 * \param T : The list of triangle to add. Each triangle is described as a list of three indices
	 * \param he : First halfedge handle of the facet
	 */
	void add_triangle(std::vector<std::vector<unsigned long> >& T, Halfedge_handle& he)
	{
		//For each triangle of the vector T...
		for(unsigned int i = 0; i != T.size(); ++i) {
			//we verify that the indices are valid.
			//if one of the indices equals to 0xFFFFFFFF, 0xFFFFFFFE or 0XFFFFFFFD, it means that
			//the corresponding vertex is respectively the first, the second or the third vertex
			//of the facet.
			for(unsigned int j = 0 ; j != 3 ; ++j) {
				switch(T[i][j]) {
				case 0xFFFFFFFF:
					if(he->vertex()->Label != 0xFFFFFFFF) {
						T[i][j] = he->vertex()->Label;
					} else {
						T[i][j] = m_Sorted_vertices.size();
						he->vertex()->Label = T[i][j];
						m_Sorted_vertices.push_back(he->vertex()->point());
					}
					break;
				case 0xFFFFFFFE:
					if(he->next()->vertex()->Label != 0xFFFFFFFF) {
						T[i][j] = he->next()->vertex()->Label;
					} else {
						T[i][j] = m_Sorted_vertices.size();
						he->next()->vertex()->Label = T[i][j];
						m_Sorted_vertices.push_back(he->next()->vertex()->point());
					}
					break;
				case 0xFFFFFFFD:
					if(he->next()->next()->vertex()->Label != 0xFFFFFFFF) {
						T[i][j] = he->next()->next()->vertex()->Label;
					} else {
						T[i][j] = m_Sorted_vertices.size();
						he->next()->next()->vertex()->Label = T[i][j];
						m_Sorted_vertices.push_back(he->next()->next()->vertex()->point());
					}
					break;
				}
			}

			//finally, the facet is added to the list of the facets
			m_Facets_indices.push_back(T[i]);
		}
	}

	/*!
	 * \brief Adds a Vertex
	 * \param p : The point to add
	 * \param l : The corresponding label
	 */
	void add_vertex(Point_3 p, unsigned long& l)    // MT: suppression référence
	{
		//The value of the label is updated
		l = m_Sorted_vertices.size();
		//The vertex is added
		m_Sorted_vertices.push_back(p);
	}

	/*!
	 * \brief this method builds the polyhedron, using the vertices and the facets stored
	 * \param hds : The halfedge data structure
	 */
	void operator()(HDS& hds)
	{
		Builder B(hds, true);
		B.begin_surface(3,1);
		add_vertices(B);
		add_facets(B);
		B.end_surface();
	}

private:

	/*!
	 * \brief Used to build the vertices of the polyhedron
	 * \param B : The builder
	 */
	void add_vertices(Builder& B)
	{
		for(int i = 0; i != (int)this->m_Sorted_vertices.size(); i++) {
			B.add_vertex(this->m_Sorted_vertices[i]);
		}
	}

	/*!
	 * \brief Used to build the facets of the polyhedron
	 * \param B : The builder
	 */
	void add_facets(Builder& B)
	{
		for(int i = 0; i != (int)this->m_Facets_indices.size(); i++) {
			B.begin_facet();
			for(int j = 0; j != (int)this->m_Facets_indices[i].size(); j++) {
				B.add_vertex_to_facet(this->m_Facets_indices[i][j]);
			}
			B.end_facet();
		}
	}
};

#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_triangulation_2.h>

/*!
 * \class Enriched_vertex_base
 * \brief Enriches the vertices of a triangulation
 */
template <class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class Enriched_vertex_base : public Vb
{
public:
	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
		typedef Enriched_vertex_base<Gt,Vb2> Other;
	};

private:
	/*! \brief An Id for the vertex*/
	unsigned long m_Label;

public:
	/*! \brief Accessor
	 * \param Label : The value to assign*/
	void set_Label(unsigned long Label)
	{
		m_Label = Label;
	}
	/*! \brief Accessor
	 * \return The Label of the vertex*/
	unsigned long get_Label()
	{
		return m_Label;
	}
};

/*!
 * \class Enriched_face_base
 * \brief Enriches the faces of a triangulation
 */
template <class Gt, class Fb = CGAL::Constrained_triangulation_face_base_2<Gt> >
class Enriched_face_base : public Fb
{
public:
	/*!
	 * \typedef typename Vertex_handle
	 * \brief Handle for a vertex of a triangulation
	 */
	typedef typename Fb::Triangulation_data_structure::Vertex_handle          Vertex_handle;

	/*!
	 * \typedef typename Face_handle
	 * \brief Handle for a face of a triangulation
	 */
	typedef typename Fb::Triangulation_data_structure::Face_handle            Face_handle;

	template < typename TDS2 >
	struct Rebind_TDS {
		typedef typename Fb::template Rebind_TDS<TDS2>::Other	Fb2;
		typedef Enriched_face_base<Gt,Fb2>						Other;
	};

private:
	/*! \brief True if the vertex has been processed*/
	bool m_OK;
	/*! \brief True if the vertex belongs to the result*/
	bool m_Ext;

public:
	Enriched_face_base() : Fb() {}
	Enriched_face_base(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2) : Fb(v0,v1,v2) {}
	Enriched_face_base(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
					   Face_handle n0, Face_handle n1, Face_handle n2) : Fb(v0,v1,v2,n0,n1,n2) {}
	Enriched_face_base(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Face_handle n0, Face_handle n1, Face_handle n2,
					   bool c0, bool c1, bool c2) : Fb(v0,v1,v2,n0,n1,n2) {}
	/*! \brief Accessor
	 * \param Ext : The value to assign*/
	void set_Ext(bool Ext)
	{
		m_Ext = Ext;
	}
	/*! \brief Accessor
	 * \return true if the triangle of the triangulation belongs to the result*/
	bool get_Ext()
	{
		return m_Ext;
	}
	/*! \brief Accessor
	 * \param OK : The value to assign*/
	void set_OK(bool OK)
	{
		m_OK = OK;
	}
	/*! \brief Accessor
	 * \return true if the parameter m_Ext has been determined*/
	bool get_OK()
	{
		return m_OK;
	}
};

/*!
 * \class Triangulation
 * \brief To subdivide a facet. (the kernel K must be exact)
 */
template <class K>
class Triangulation
{
	typedef typename K::Vector_3	Vector_exact;
	/*!
	 * \typedef typename Point_3
	 * \brief 3d point using exact number type
	 */
	typedef typename K::Point_3														Point_3;

	/*!
	 * \typedef typename Tri_vb
	 * \brief Vertex base
	 */
	typedef /*typename*/ Enriched_vertex_base<K>										Tri_vb;

	/*!
	 * \typedef typename Tri_fb
	 * \brief Face base
	 */
	typedef /*typename*/ Enriched_face_base<K>											Tri_fb;

	/*!
	 * \typedef typename Tri_DS
	 * \brief Data structure of the triangulation
	 */
	typedef typename CGAL::Triangulation_data_structure_2<Tri_vb,Tri_fb>			Tri_DS;

	/*!
	 * \typedef typename Itag
	 * \brief No intersection tag
	 */
	typedef typename CGAL::No_intersection_tag										Itag;

	/*!
	 * \typedef typename Constrained_tri
	 * \brief 2d constrained triangulation
	 */
	typedef typename CGAL::Constrained_triangulation_2<K, Tri_DS, Itag>	Constrained_tri;

	/*!
	 * \typedef typename Vertex_handle_tri
	 * \brief Vertex handle for the triangulation
	 */
	typedef typename Constrained_tri::Vertex_handle						Vertex_handle_tri;

	/*!
	 * \typedef typename Face_handle_tri
	 * \brief Face handle for the triangulation
	 */
	typedef typename Constrained_tri::Face_handle							Face_handle_tri;

	/*!
	 * \typedef typename Point_tri
	 * \brief 2d point for the triangulation
	 */
	typedef typename Constrained_tri::Point								Point_tri;

	/*!
	 * \typedef typename Face_iterator_tri
	 * \brief Iterator for the faces of the triangulation
	 */
	typedef typename Constrained_tri::Face_iterator						Face_iterator_tri;


	typedef typename CGAL::Polyhedron_3<K,Enriched_items> M;
	typedef typename M::Halfedge_handle Halfedge_handle;
public:
	/*!
	 * \brief Constructor
	 * \param he : A halfedge incident to the facet
	 * \param norm_dir : The vector directing the normal of the facet
	 */
	Triangulation(Halfedge_handle& he, Vector_exact& norm_dir)
	{
		//find the longest coordinate of the normal vector, and its sign
		double x = to_double(norm_dir.x());
		double y = to_double(norm_dir.y());
		double z = to_double(norm_dir.z());
		double absx = std::abs(x);
		double absy = std::abs(y);
		double absz = std::abs(z);

		//this information is stored using a code :
		//0 : The coordinate X is the longest, and positive
		//1 : The coordinate Y is the longest, and positive
		//2 : The coordinate Z is the longest, and positive
		//3 : The coordinate X is the longest, and negative
		//4 : The coordinate Y is the longest, and negative
		//5 : The coordinate Z is the longest, and negative
		if(absx >= absy && absx >= absz) max_coordinate = (x>0)?0:3;
		else if(absy >= absx && absy >= absz) max_coordinate = (y>0)?1:4;
		else if(absz >= absx && absz >= absy) max_coordinate = (z>0)?2:5;

		//we add the three vertices of the facet to the triangulation
		//The Label of these vertices is set for the corresponding point in the triangulation
		v1 = add_new_pt(he->vertex()->point(), he->vertex()->Label);
		v2 = add_new_pt(he->next()->vertex()->point(), he->next()->vertex()->Label);
		v3 = add_new_pt(he->next()->next()->vertex()->point(), he->next()->next()->vertex()->Label);

		//if the vertices does not have an Id (Label = OxFFFFFFFF), the labels
		//of the points in the triangulation is set as follows :
		//0xFFFFFFFF for the first point
		//0xFFFFFFFE for the second point
		//0xFFFFFFFD for the third point
		if(v2->get_Label() == 0xFFFFFFFF) v2->set_Label(0xFFFFFFFE);
		if(v3->get_Label() == 0xFFFFFFFF) v3->set_Label(0xFFFFFFFD);
	}

	/*!
	 * \brief Compute the orthogonal projection of a point to the plane defined by the longest coordinate of the normal vector
	 * \param p : the point (in 3d)
	 * \return The projection as a 2d point
	 */
	Point_tri get_minvar_point_2(Point_3& p)
	{
		switch(max_coordinate) {
		case 0:
			return Point_tri(p.y(),p.z());
			break;
		case 1:
			return Point_tri(p.z(),p.x());
			break;
		case 2:
			return Point_tri(p.x(),p.y());
			break;
		case 3:
			return Point_tri(p.z(),p.y());
			break;
		case 4:
			return Point_tri(p.x(),p.z());
			break;
		case 5:
			return Point_tri(p.y(),p.x());
			break;
		default:
			return Point_tri(p.y(),p.z());
		}
	}

	/*!
	 * \brief Adds a point in the triangulation
	 * \param p : The point (in 3d : the projection in 2d is done automatically)
	 * \param Label : The label of the point
	 * \return The Vertex_handle of the point added
	 */
	Vertex_handle_tri add_new_pt(Point_3 p, unsigned long& Label)   // MT: suppression référence
	{
		//if the point is not a new one, we verify that the point has not already been added
		if(Label != 0xFFFFFFFF)
			for(unsigned int i = 0; i != pts_point.size(); ++i)
				if(Label == pts_point[i])
					//if the point is already in the triangulation, we return its handle
					return pts_vertex[i];
		Vertex_handle_tri v;
		v = ct.insert(get_minvar_point_2(p));
		v->set_Label(Label);
		pts_point.push_back(Label);
		pts_vertex.push_back(v);
		return v;
	}

	/*!
	 * \brief Adds a constrained segment in the triangulation
	 * \param p1 : The first point (in 3d)
	 * \param p2 : The second point (in 3d)
	 * \param Label1 : The label of the first point
	 * \param Label2 : The label of the second point
	 */
	void add_segment(Point_3& p1, Point_3& p2, unsigned long& Label1, unsigned long& Label2)
	{
		// we add the two points in the triangulation and store their handles in c1 and c2
		c1 = add_new_pt(p1, Label1);
		c2 = add_new_pt(p2, Label2);
		// we set a constrained segment between these two points
		ct.insert_constraint(c1, c2);
		// if an other segment is added, we will overwrite c1 and c2
		// what is important is to memorize the handles of the last segment added
	}

	/*!
	 * \brief Gets the triangles of the triangulation that belongs to the result
	 * and deduce for the three neighboring facets if they belong to the result or not
	 * \param inv_triangles : must be true if the orientation of the triangles must be inverted
	 * \param IsExt : Pointer on a three-case boolean table.
	 * \return The list of the triangles belonging to the result.
	 * each triangle is defined by a list of three labels
	 */
	std::vector<std::vector<unsigned long> > get_triangles(bool inv_triangles, bool* IsExt)
	{
		//init
		IsExt[0] = false;
		IsExt[1] = false;
		IsExt[2] = false;
		std::vector<std::vector<unsigned long> > tris;
		for(Face_iterator_tri fi = ct.faces_begin(); fi != ct.faces_end(); fi++)
			fi->set_OK(false);

		//the constrained segments are oriented, we search the triangle (c1, c2, X), (X, c1, c2) or (c2, X, c1)
		//where c1 and c2 are the two points related to the last constrained segment added, and X another point
		//this triangle belongs to the result (thanks to the orientation of the segments)
		Face_handle_tri f, f2 = c1->face();
		int i=0; // MT
		do {
			f = f2;
			f->has_vertex(c1,i);
			f2 = f->neighbor(f->ccw(i));
		} while(!(f->has_vertex(c2) && f2->has_vertex(c2)));

		//dans le cas particulier ou la frontiere se trouve exactement sur un bord de la triangulation,
		//et que ce triangle n'appartient pas a la triangulation, on démarrera avec l'autre triangle
		//incluant le segment c1, c2 (et donc, n'appartenant pas au résultat

		//if the segment is exactly on the border of the triangulation, the triangle could be outside the triangulation
		//in that case, we will search the other triangle including the points c1 and c2
		//this triangle does not belong to the result
		if(f->has_vertex(ct.infinite_vertex())) {
			f = f2;
			f->set_Ext(false);
		} else {
			f->set_Ext(true);
		}

		std::stack<Face_handle_tri> sfh;
		f->set_OK(true);
		sfh.push(f);

		//we decide for all the triangles, if they belongs to the result, starting from the first triangle f,
		//by moving on the triangulation using the connectivity between the triangles.
		//If a constrained segment is crossed, the value of the tag "isext" is inverted
		while(!sfh.empty()) {
			f = sfh.top();
			sfh.pop();

			if(f->get_Ext()) {
				std::vector<unsigned long> tri;
				int i;
				tri.push_back(f->vertex(0)->get_Label());

				//verify if the neighboring facets belongs to the result or not
				if(f->has_vertex(v1,i) && f->neighbor(f->ccw(i))->has_vertex(ct.infinite_vertex())) IsExt[0] = true;
				if(f->has_vertex(v2,i) && f->neighbor(f->ccw(i))->has_vertex(ct.infinite_vertex())) IsExt[1] = true;
				if(f->has_vertex(v3,i) && f->neighbor(f->ccw(i))->has_vertex(ct.infinite_vertex())) IsExt[2] = true;

				if(inv_triangles) {
					tri.push_back(f->vertex(2)->get_Label());
					tri.push_back(f->vertex(1)->get_Label());
				} else {
					tri.push_back(f->vertex(1)->get_Label());
					tri.push_back(f->vertex(2)->get_Label());
				}
				tris.push_back(tri);
			}
			for(i = 0; i!=3; i++) {
				if(!(f->neighbor(i)->get_OK() || f->neighbor(i)->has_vertex(ct.infinite_vertex()))) {
					f->neighbor(i)->set_OK(true);
					f->neighbor(i)->set_Ext((f->is_constrained(i))?!f->get_Ext():f->get_Ext());
					sfh.push(f->neighbor(i));
				}
			}
		}
		return tris;
	}

	/*!
	 * \brief Gets all the triangles of the triangulation
	 * \param inv_triangles : must be true if the orientation of the triangles must be inverted
	 * \return The list of the triangles belonging to the result.
	 * each triangle is defined by a list of three labels
	 */
	std::vector<std::vector<unsigned long> > get_all_triangles(bool inv_triangles)
	{
		std::vector<std::vector<unsigned long> > tris;
		for(Face_iterator_tri f = ct.faces_begin(); f != ct.faces_end(); f++) {
			std::vector<unsigned long> tri;
			tri.push_back(f->vertex(0)->get_Label());
			if(inv_triangles) {
				tri.push_back(f->vertex(2)->get_Label());
				tri.push_back(f->vertex(1)->get_Label());
			} else {
				tri.push_back(f->vertex(1)->get_Label());
				tri.push_back(f->vertex(2)->get_Label());
			}
			tris.push_back(tri);
		}
		return tris;
	}

private:
	/*! \brief The triangulation*/
	Constrained_tri ct;
	/*! \brief List of the id of the points added in the triangulation*/
	std::vector<InterId> pts_point;
	/*! \brief List of the handles of the points added in the triangulation*/
	std::vector<Vertex_handle_tri> pts_vertex;
	/*! \brief Handle of the point corresponding to the first vertex of the facet*/
	Vertex_handle_tri v1;
	/*! \brief Handle of the point corresponding to the second vertex of the facet*/
	Vertex_handle_tri v2;
	/*! \brief Handle of the point corresponding to the third vertex of the facet*/
	Vertex_handle_tri v3;
	/*! \brief Handle of the point corresponding to the first vertex of the last segment added*/
	Vertex_handle_tri c1;
	/*! \brief Handle of the point corresponding to the second vertex of the last segment added*/
	Vertex_handle_tri c2;
	/*! \brief Code identifying the plane where the triangulation is done \n
	 * 0 : Plane (y, z) \n
	 * 1 : Plane (z, x) \n
	 * 2 : Plane (x, y) \n
	 * 3 : Plane (z, y) \n
	 * 4 : Plane (x, z) \n
	 * 5 : Plane (y, x)
	 */
	int max_coordinate;
};

#ifdef BOOLEAN_OPERATIONS_DEBUG
#include "Time_measure.h"
#endif // BOOLEAN_OPERATIONS_DEBUG



/*! \class Enriched_Triangle
 * \brief An enriched triangle*/
template <typename Kernel,typename AABB_Kernel>
class Enriched_Triangle : public AABB_Kernel::Triangle_3
{
public:
	typedef typename AABB_Kernel::Point_3 Point_3;
	typedef typename CGAL::Polyhedron_3<Kernel,Enriched_items>::Facet_handle Facet_handle;

	Enriched_Triangle(Facet_handle& _f)
		: AABB_Kernel::Triangle_3(to_K(_f->facet_begin()->vertex()->point() + (_f->facet_begin()->vertex()->point() - _f->facet_begin()->next()->vertex()->point()) / 1000),
								  to_K(_f->facet_begin()->next()->vertex()->point() + (_f->facet_begin()->next()->vertex()->point() - _f->facet_begin()->next()->next()->vertex()->point()) / 1000),
								  to_K(_f->facet_begin()->next()->next()->vertex()->point() + (_f->facet_begin()->next()->next()->vertex()->point() - _f->facet_begin()->vertex()->point()) / 1000)
								 ), f(_f) {}


	/*! \brief Accessor
	 * \return The handle of the facet used to build the triangle
	 */
	Facet_handle facet()
	{
		return f;
	}

	/*! \brief convert any 3d point in the kernel used by the AABB-tree
	 * \param p : The point to convert
	 * \return The 3d point converted
	 */
	inline Point_3 to_K(typename Kernel::Point_3 p)
	{
		return Point_3(to_double(p.x()),to_double(p.y()),to_double(p.z()));
	}  // MT: suppression référence

private:
	/*! \brief The handle of the facet used to build the triangle*/
	Facet_handle f;
};

/*! \class BoolPolyhedra
 * \brief The class that compute a Boolean operation*/
template <typename Kernel, typename Items>
class BoolPolyhedra
{
	typedef CGAL::Simple_cartesian<double> AABB_Kernel;
	typedef Enriched_Triangle<Enriched_kernel,AABB_Kernel> Triangle;
	typedef CGAL::AABB_triangle_primitive<AABB_Kernel,std::list<Triangle>::iterator> AABB_Primitive;
	typedef CGAL::AABB_traits<AABB_Kernel, AABB_Primitive> AABB_Traits;
	typedef CGAL::AABB_tree<AABB_Traits> AABB_Tree;

	typedef typename Kernel::FT	FT;
	typedef typename Kernel::Point_3 Point_3;
	typedef typename Kernel::Vector_3 Vector_3;

	typedef typename CGAL::Polyhedron_3<Kernel,Items> Polyhedron_3;
	typedef typename Polyhedron_3::Vertex Vertex;
	typedef typename Polyhedron_3::Facet Facet;
	typedef typename Polyhedron_3::HalfedgeDS HDS;

	typedef typename Polyhedron_3::Vertex_iterator Vertex_iterator;
	typedef typename Polyhedron_3::Facet_iterator Facet_iterator;

	typedef typename Polyhedron_3::Halfedge_handle Halfedge_handle;
	typedef typename Polyhedron_3::Facet_handle Facet_handle;

	typedef typename Vertex::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
	typedef typename Facet::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
private:
	/*! \struct Triangle_Cut
	 * \brief A structure containing informations about an intersected facet*/
	struct Triangle_Cut {
		/*! \brief true if the facet belongs to the first polyhedron*/
		bool								Facet_from_A;
		/*! \brief An exact vector giving the direction of the normal*/
		Vector_3						norm_dir;
		/*! \brief A list of segments (the intersections with the facets of the other polyhedron)*/
		std::vector<std::vector<InterId> >	CutList;
		/*! \brief A list of points (when the intersection is a point)*/
		std::set<InterId>					PtList;
		/*! \brief The list of the intersections*/
		std::map<HalfedgeId, InterId>		RefInter;

		/*! \brief Default constructor*/
		Triangle_Cut() {}
		/*! \brief Constructor
		 \param V : The normal direction
		 \param ffA : Must be true if the facet belongs to the first polyhedron*/
		Triangle_Cut(Vector_3 V, bool ffA)
		{
			norm_dir=V;    // MT
			Facet_from_A=ffA;
		}
	};

	/*! \struct Info_Inter
	 * \brief Contains informations about an intersection between a facet and a halfedge*/
	struct Info_Inter {
		/*! \brief The facet*/
		Facet_handle		f;
		/*! \brief The halfedge*/
		Halfedge_handle		he;
		/*! \brief true if the intersection is exactly on the vertex pointed by he*/
		bool				IsOnVertex;
		/*! \brief A code for the location of the intersection :\n\n
		 * 0 : Intersection is strictly in the facet\n
		 * 1 : Intersection is on the first edge of the facet\n
		 * 2 : Intersection is on the second edge of the facet\n
		 * 3 : Intersection is exactly on the first vertex of the facet\n
		 * 4 : Intersection is on the third edge of the facet\n
		 * 5 : Intersection is exactly on the third vertex of the facet\n
		 * 6 : Intersection is exactly on the second vertex of the facet\n
		 * 7 : There is no intersection */
		unsigned short		res;
		/*! \brief The intersection point (exact)*/
		Point_3		pt;
		/*! \brief The Id of the intersection point*/
		InterId				Id;
	};

public:
	/*! \brief Constructor.
	 * \brief Computes a boolean operation
	 * \param pMA : The first polyhedron
	 * \param pMB : The second polyhedron
	 * \param pMout : The result polyhedron
	 * \param BOOP : The Boolean operator. Must be UNION, INTER or MINUS*/
	BoolPolyhedra(MEPP_Polyhedron*& pMA, MEPP_Polyhedron*& pMB, MEPP_Polyhedron*& pMout, Bool_Op BOOP) : m_BOOP(BOOP)
	{

#ifdef BOOLEAN_OPERATIONS_DEBUG
		std::ofstream ofstrMA("input_A_boolsum.off");
		ofstrMA << *pMA;
		std::ofstream ofstrMB("input_B_boolsum.off");
		ofstrMB << *pMB;
		Time_measure Timer_total, Timer;

		Timer_total.Start();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		Init(pMA, pMB);

#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_Init = Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		FindCouples();

#ifdef BOOLEAN_OPERATIONS_DEBUG
		duration_FindCouples = Timer.GetDiff();
		Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

		if(!m_Couples.empty()) {
			ComputeIntersections();

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_ComputeIntersections = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			CutIntersectedFacets();

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_CutIntersectedFacets = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			PropagateFacets();

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_PropagateFacets = Timer.GetDiff();
			Timer.Start();
#endif // BOOLEAN_OPERATIONS_DEBUG

			pMout->delegate(ppbuilder);

#ifdef BOOLEAN_OPERATIONS_DEBUG
			duration_delegate = Timer.GetDiff();
			duration_total = Timer_total.GetDiff();
			std::ofstream ofstrMS("output_boolsum.off");
			ofstrMS << *pMout;
			ofstrtime.open("output_time.txt");
			WriteData(pMout);
			ColorType();
#endif // BOOLEAN_OPERATIONS_DEBUG

		}
	}

	/*! \brief Destructor*/
	~BoolPolyhedra() {}

private:
	/**
	 * \fn inline Vector_exact Compute_Normal_direction(Halfedge_handle &he)
	 * \brief Compute a vector in the same direction as the normal vector
	 *
	 * \param he : A Halfedge incident to the facet
	 * \return The normal direction (exact).
	 */
	inline Vector_3 Compute_Normal_direction(Halfedge_handle he)   // MT: suppression référence
	{
		return CGAL::cross_product(he->next()->vertex()->point() - he->vertex()->point(),
								   he->next()->next()->vertex()->point() - he->vertex()->point());
	}
	/*! \brief Initialisation of the tags, and triangulation of the two input polyhedra
	 * \param pMA : The first polyhedron
	 * \param pMB : The second polyhedron*/
	void Init(MEPP_Polyhedron*& pMA, MEPP_Polyhedron*& pMB)
	{
		m_pA = pMA;
		m_pB = pMB;

		//triangulation of the two input polyhedra
		//this is necessary for the AABB-tree, and simplify the computation of the intersections
		if(!m_pA->is_pure_triangle()) m_pA->triangulate();
		if(!m_pB->is_pure_triangle()) m_pB->triangulate();

		//initialize the tags
		for(Vertex_iterator pVertex = m_pA->vertices_begin(); pVertex != m_pA->vertices_end(); ++pVertex) {
			pVertex->Label = 0xFFFFFFFF;
		}
		for(Vertex_iterator pVertex = m_pB->vertices_begin(); pVertex != m_pB->vertices_end(); ++pVertex) {
			pVertex->Label = 0xFFFFFFFF;
		}
		for(Facet_iterator pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); ++pFacet) {
			pFacet->Label = 0xFFFFFFFF;
			pFacet->IsExt = false;
			pFacet->IsOK = false;
		}
		for(Facet_iterator pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); ++pFacet) {
			pFacet->Label = 0xFFFFFFFF;
			pFacet->IsExt = false;
			pFacet->IsOK = false;
		}
	}

	/*! \brief Finds every couple of facets between the two input polyhedra that intersects
	 * \brief Each couple is stored in the member m_Couples*/
	void FindCouples()
	{
		//A AABB-tree is built on the facets of one of the polyhedra. A collision test is done with each facet of the other polyhedron.
		Facet_iterator pFacet =	NULL;
		std::list<AABB_Tree::Primitive_id> primitives;
		std::list<Triangle> triangles;

		HalfedgeId i = 0;
		FacetId j = 0;

		//The AABB-tree is built on the polyhedron with the less number of facets
		if(m_pA->size_of_facets() < m_pB->size_of_facets()) {
			//Building the AABB-tree on the first polyhedron
			for(pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++) triangles.push_back(Triangle(pFacet));
			tree.rebuild(triangles.begin(),triangles.end());

			//collision test with each facet of the second polyhedron
			for(pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++) {
				//"primitives" is the list of the triangles intersected (as a list of triangles)
				tree.all_intersected_primitives(Triangle(pFacet), std::back_inserter(primitives));
				if(primitives.size() !=0) {
					Facet_Handle.push_back(pFacet);
					//update of the tags (the facet and the three incidents halfedges
					pFacet->Label = j++;
					pFacet->facet_begin()->Label = i++;
					pFacet->facet_begin()->next()->Label = i++;
					pFacet->facet_begin()->next()->next()->Label = i++;
					//creation of a Triangle_Cut structure to store the informations about the intersections
					Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(pFacet->facet_begin()), false));
					do {
						//same operations for the intersected primitives (only one time)
						if(primitives.back()->facet()->Label == 0xFFFFFFFF) {
							Facet_Handle.push_back(primitives.back()->facet());
							primitives.back()->facet()->Label = j++;
							primitives.back()->facet()->facet_begin()->Label = i++;
							primitives.back()->facet()->facet_begin()->next()->Label = i++;
							primitives.back()->facet()->facet_begin()->next()->next()->Label = i++;
							Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(primitives.back()->facet()->facet_begin()), true));
						}
						//store every couple of intersected facet
						m_Couples[primitives.back()->facet()->Label].insert(pFacet->Label);
						primitives.pop_back();
					} while(primitives.size() != 0);
				}
			}
		} else {
			//Building the AABB-tree on the second polyhedron
			for(pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++) triangles.push_back(Triangle(pFacet));
			tree.rebuild(triangles.begin(),triangles.end());

			//collision test with each facet of the first polyhedron
			for(pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++) {
				//"primitives" is the list of the triangles intersected (as a list of triangles)
				tree.all_intersected_primitives(Triangle(pFacet), std::back_inserter(primitives));
				if(primitives.size() !=0) {
					Facet_Handle.push_back(pFacet);
					//update of the tags (the facet and the three incidents halfedges
					pFacet->Label = j++;
					pFacet->facet_begin()->Label = i++;
					pFacet->facet_begin()->next()->Label = i++;
					pFacet->facet_begin()->next()->next()->Label = i++;
					//creation of a Triangle_Cut structure to store the informations about the intersections
					Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(pFacet->facet_begin()), true));
					do {
						//same operations for the intersected primitives (only one time)
						if(primitives.back()->facet()->Label == 0xFFFFFFFF) {
							Facet_Handle.push_back(primitives.back()->facet());
							primitives.back()->facet()->Label = j++;
							primitives.back()->facet()->facet_begin()->Label = i++;
							primitives.back()->facet()->facet_begin()->next()->Label = i++;
							primitives.back()->facet()->facet_begin()->next()->next()->Label = i++;
							Inter_tri.push_back(Triangle_Cut(Compute_Normal_direction(primitives.back()->facet()->facet_begin()), false));
						}
						//store every couple of intersected facet
						m_Couples[pFacet->Label].insert(primitives.back()->facet()->Label);
						primitives.pop_back();
					} while(primitives.size() != 0);
				}
			}
		}
	}

	/*! \brief Compute the intersections*/
	void ComputeIntersections()
	{
		while(!m_Couples.empty()) {
			FacetId fA, fB;
			fA = m_Couples.begin()->first;
			fB = *m_Couples[fA].begin();
			InterTriangleTriangle(fA, fB);
			rmCouple(fA, fB);
		}
	}


	/*! \brief Cuts the intersected facets and starts to build the result*/
	void CutIntersectedFacets()
	{
		Triangle_Cut TriCut;
		Halfedge_handle he;

		//every intersected facet is triangulated if at least one of the intersections is a segment
		for(FacetId Facet = 0 ; Facet != Inter_tri.size() ; ++Facet) {
			if(!Inter_tri[Facet].CutList.empty()) {
				TriCut = Inter_tri[Facet];
				he = Facet_Handle[Facet]->facet_begin();
				bool IsExt[3];
				//creation of a triangulation
				Triangulation<Kernel> T(he, TriCut.norm_dir);
				//add the list of intersection points (only happens in case of intersection of two edges)
				for(std::set<InterId>::iterator i = TriCut.PtList.begin(); i != TriCut.PtList.end(); ++i) {
					T.add_new_pt(InterPts[*i], (unsigned long&)*i);     // MT: ajout cast
				}
				//add the intersection segments
				for(int i = 0; i!=(int)TriCut.CutList.size(); ++i) {
					T.add_segment(InterPts[TriCut.CutList[i][0]], InterPts[TriCut.CutList[i][1]], TriCut.CutList[i][0], TriCut.CutList[i][1]);
				}
				//get the triangles of the triangulation thay belong to the result
				//and determine if the three neighboring facets belongs to the result (using IsExt[3])
				std::vector<std::vector<unsigned long> > Tri_set = T.get_triangles((m_BOOP == MINUS && !TriCut.Facet_from_A)?true:false, IsExt);
				//add these triangles to the result
				ppbuilder.add_triangle(Tri_set, he);

				//update the tags
				Facet_Handle[Facet]->IsOK = true;
				if(IsExt[0]) he->opposite()->facet()->IsExt = true;
				if(IsExt[1]) he->next()->opposite()->facet()->IsExt = true;
				if(IsExt[2]) he->next()->next()->opposite()->facet()->IsExt = true;
			}
		}
	}

	/*! \brief Complete the building of the result*/
	void PropagateFacets()
	{
		Facet_handle pFacet = NULL, f = NULL, nf = NULL;
		std::stack<Facet_handle> tmpTriangles;

		//add to a stack the intersected facets that have been cut during CutIntersectedFacets
		for(pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++) {
			if(pFacet->IsOK) tmpTriangles.push(pFacet);
		}

		//while the stack is not empty, we look the three neighboring facets
		//if these facets has not been validated (IsOK == false), the facet is validated and added to the stack
		//if this facet is taged as a part of the result (IsExt == true), the facet is added to the result
		while(!tmpTriangles.empty()) {
			f = tmpTriangles.top();
			tmpTriangles.pop();
			nf = f->facet_begin()->opposite()->facet();
			if(!nf->IsOK) {
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt) add_facet_to_solution(nf, true);
			}
			nf = f->facet_begin()->next()->opposite()->facet();
			if(!nf->IsOK) {
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt) add_facet_to_solution(nf, true);
			}
			nf = f->facet_begin()->next()->next()->opposite()->facet();
			if(!nf->IsOK) {
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt) add_facet_to_solution(nf, true);
			}
		}

		//same process for the second polyhedron
		for(pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++) {
			if(pFacet->IsOK) tmpTriangles.push(pFacet);
		}

		while(!tmpTriangles.empty()) {
			f = tmpTriangles.top();
			tmpTriangles.pop();
			nf = f->facet_begin()->opposite()->facet();
			if(!nf->IsOK) {
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt) add_facet_to_solution(nf, false);
			}
			nf = f->facet_begin()->next()->opposite()->facet();
			if(!nf->IsOK) {
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt) add_facet_to_solution(nf, false);
			}
			nf = f->facet_begin()->next()->next()->opposite()->facet();
			if(!nf->IsOK) {
				nf->IsOK = true;
				tmpTriangles.push(nf);
				if(nf->IsExt) add_facet_to_solution(nf, false);
			}
		}
	}


	/*! \brief removes properly a couple from the list
	 * \param A : Id of the first facet
	 * \param B : Id of the second facet
	 */
	void rmCouple(FacetId& A, FacetId& B)
	{
		if(m_Couples[A].count(B) != 0) m_Couples[A].erase(B);
		if(m_Couples[A].empty()) m_Couples.erase(A);
	}

	/*! \brief Compute the intersection between two facets
	 * \param A : Facet Id of the first facet (from the first polyhedron)
	 * \param B : Facet Id of the second facet (from the second polyhedron)*/
	void InterTriangleTriangle(FacetId& A, FacetId& B)
	{
		Vector_3 nA, nB;
		nA = Inter_tri[A].norm_dir;
		nB = Inter_tri[B].norm_dir;

		Facet_handle fA, fB, fA2, fB2;
		fA = Facet_Handle[A];
		fB = Facet_Handle[B];

		Halfedge_handle heA[3], heB[3];
		heA[0] = fA->facet_begin();
		heA[1] = heA[0]->next();
		heA[2] = heA[1]->next();
		heB[0] = fB->facet_begin();
		heB[1] = heB[0]->next();
		heB[2] = heB[1]->next();

		Point_3 ptA[3], ptB[3];
		ptA[0] = heA[0]->vertex()->point();
		ptA[1] = heA[1]->vertex()->point();
		ptA[2] = heA[2]->vertex()->point();
		ptB[0] = heB[0]->vertex()->point();
		ptB[1] = heB[1]->vertex()->point();
		ptB[2] = heB[2]->vertex()->point();

		//compute the position of the three vertices of each triangle regarding the plane of the other
		//positive if the vertex is above
		//negative if the vertex is under
		//zero if the vertex is exactly on the triangle
		FT posA[3], posB[3];
		posA[0] = nB * (ptA[0] - ptB[0]);
		posA[1] = nB * (ptA[1] - ptB[0]);
		posA[2] = nB * (ptA[2] - ptB[0]);
		posB[0] = nA * (ptB[0] - ptA[0]);
		posB[1] = nA * (ptB[1] - ptA[0]);
		posB[2] = nA * (ptB[2] - ptA[0]);

		//a code is computed on 6 bits using these results (two bits for each point)
		//10 -> above ; 01 -> under ; 00 -> on the plane
		unsigned short posAbin, posBbin;
		posAbin = ((posA[0] > 0)? 32 : 0)
				  + ((posA[0] < 0)? 16 : 0)
				  + ((posA[1] > 0)? 8 : 0)
				  + ((posA[1] < 0)? 4 : 0)
				  + ((posA[2] > 0)? 2 : 0)
				  + ((posA[2] < 0)? 1 : 0);

		posBbin = ((posB[0] > 0)? 32 : 0)
				  + ((posB[0] < 0)? 16 : 0)
				  + ((posB[1] > 0)? 8 : 0)
				  + ((posB[1] < 0)? 4 : 0)
				  + ((posB[2] > 0)? 2 : 0)
				  + ((posB[2] < 0)? 1 : 0);

		//if the intersection is not a segment, the intersection is not computed
		//the triangles intersects on a point (one vertex on the plane and the two others under or above
		if(posAbin == 5 || posAbin == 10 || posAbin == 17 || posAbin == 34 || posAbin == 20 || posAbin == 40
		   || posBbin == 5 || posBbin == 10 || posBbin == 17 || posBbin == 34 || posBbin == 20 || posBbin == 40) return;
		//no possible intersection (one of the triangle is completely under or above the other
		if(posAbin == 42 || posAbin == 21
		   || posBbin == 42 || posBbin == 21) return;
		//the triangles are coplanar
		if(posAbin == 0) return;

		//if an edge of a triangle is on the plane of the other triangle, it is necessary to verify if the
		//two polyhedra are intersecting on these edges, or if it only is a contact to know if the intersection
		//between the triangles must be computed or not.
		//"edgeA" and "edgeB" are codes
		//0 : the first edge is on the plane
		//1 : the second edge is on the plane
		//2 : the third edge is on the plane
		//3 : there is no edge on the plane
		unsigned short edgeA = 3, edgeB = 3;
		if(posAbin == 1  || posAbin == 2) edgeA = 1;       //points 0 and 1 on the plane
		else if(posAbin == 16 || posAbin == 32) edgeA = 2; //points 1 and 2 on the plane
		else if(posAbin == 4  || posAbin == 8) edgeA = 0;  //points 2 and 0 on the plane
		if(posBbin == 1  || posBbin == 2) edgeB = 1;       //points 0 and 1 on the plane
		else if(posBbin == 16 || posBbin == 32) edgeB = 2; //points 1 and 2 on the plane
		else if(posBbin == 4  || posBbin == 8) edgeB = 0;  //points 2 and 0 on the plane

		Vector_3 nA2, nB2;
		FT p;
		bool invert_direction = false;
		bool stop = false;

		//if an edge of the first triangle is on the plane
		if(edgeA != 3 && edgeB == 3) {
			fA2 = heA[edgeA]->opposite()->facet();
			nA2 = Inter_tri[fA2->Label].norm_dir;
			p = CGAL::cross_product(nA, nB) * CGAL::cross_product(nA2, nB);
			//if p is negative, the two triangles of the first polyhedron (including edgeA) are on the same side
			//so there is no intersection
			if(p < 0) stop = true;
			//if p == 0, fA2 is coplanar with the plane of fB
			//in that case, it is necessary to consider the boolean
			//operator used to determine if there is a contact or not
			else if(p == 0) {
				switch(m_BOOP) {
				case UNION:
					if(posA[(edgeA+1)%3] * (nA2 * nB) > 0) stop = true;
					break;
				case INTER:
					if(posA[(edgeA+1)%3] > 0) stop = true;
					break;
				case MINUS:
					if(posA[(edgeA+1)%3] * (nA2 * nB) < 0) stop = true;
					break;
				}
			}
			//the intersection between fA2 and fB is the same so this couple is removed from the list
			rmCouple(fA2->Label, fB->Label);
		}
		//if an edge of the second triangle is on the plane
		else if(edgeA == 3 && edgeB != 3) {
			fB2 = heB[edgeB]->opposite()->facet();
			nB2 = Inter_tri[fB2->Label].norm_dir;
			p = CGAL::cross_product(nA, nB) * CGAL::cross_product(nA, nB2);
			//if p is negative, the two triangles of the second polyhedron (including edgeB) are on the same side
			//so there is no intersection
			if(p < 0) stop = true;
			//if p == 0, fB2 is coplanar with the plane of fA
			//in that case, it is necessary to consider the boolean
			//operator used to determine if there is a contact or not
			else if(p == 0) {
				switch(m_BOOP) {
				case UNION:
					if(posB[(edgeB+1)%3] < 0) stop = true;
					break;
				case INTER:
					if(posB[(edgeB+1)%3] * (nB2 * nA) < 0) stop = true;
					break;
				case MINUS:
					if(posB[(edgeB+1)%3] > 0) stop = true;
					break;
				}
			}
			//the intersection between fA and fB2 is the same so this couple is removed from the list
			rmCouple(fA->Label, fB2->Label);
		}
		//if an edge of each triangle is on the plane of the other
		else if(edgeA != 3 && edgeB != 3) {
			//in this case, four triangles are concerned by the intersection
			//fA2 and fB2 are the two other concerned facets
			//we try to determine if fA and fA2 are inside or outside the second polyhedron, using fB and fB2
			bool Intersection = false;
			Vector_3 nAcnB2, nA2cnB;
			FT nAnB2, nA2nB, nA2nB2;
			FT posA2_A, posB_A, posB2_A, posB_B2, posA_B, posB2_B, posB_A2, posB2_A2, posA2_B, posA2_B2;
			Point_3 ptA2, ptB2;

			fA2 = heA[edgeA]->opposite()->facet();
			fB2 = heB[edgeB]->opposite()->facet();
			nA2 = Inter_tri[fA2->Label].norm_dir;
			nB2 = Inter_tri[fB2->Label].norm_dir;

			nAcnB2 = CGAL::cross_product(nA, nB2);
			nA2cnB = CGAL::cross_product(nA2, nB);

			nAnB2 = nA * nB2;
			nA2nB = nA2 * nB;
			nA2nB2 = nA2 * nB2;

			ptA2 = heA[edgeA]->opposite()->next()->vertex()->point();
			ptB2 = heB[edgeB]->opposite()->next()->vertex()->point();

			posA_B = posA[(edgeA+1)%3];
			posB_A = posB[(edgeB+1)%3];
			posB_A2 = nA2 * (ptB[(edgeB+1)%3] - ptA[edgeA]);
			posB_B2 = nB2 * (ptB[(edgeB+1)%3] - ptA[edgeA]);
			posA2_A = nA * (ptA2 - ptA[edgeA]);
			posA2_B = nB * (ptA2 - ptA[edgeA]);
			posA2_B2 = nB2 * (ptA2 - ptA[edgeA]);
			posB2_A = nA * (ptB2 - ptA[edgeA]);
			posB2_A2 = nA2 * (ptB2 - ptA[edgeA]);
			posB2_B = nB * (ptB2 - ptA[edgeA]);

			if(nAcnB2 == CGAL::NULL_VECTOR && nA2cnB == CGAL::NULL_VECTOR
			   && nAnB2 * nA2nB > 0) stop = true;

			//firstly, we search the position of fA
			//if fA is inside the poyhedron, Intersection = true
			if(posB_A * posB2_A > 0) { //fB and fB2 on the same side
				if(posB_B2 > 0) Intersection = true;
			} else if(posB_A * posB2_A < 0) { //fB and fB2 on opposite side
				if(posA_B < 0) Intersection = true;
			} else { //fA and fB2 coplanar
				if(posA_B * posB2_B < 0) {
					if(posB_B2 > 0) Intersection = true;
				} else {
					if(nAnB2 < 0) {
						if(m_BOOP == UNION) Intersection = true;
					} else {
						if(m_BOOP == MINUS) Intersection = true;
					}
				}
			}

			//secondly, we search the position of fA2
			//if fA2 is inside the poyhedron, "Intersection" is inverted
			if(posB_A2 * posB2_A2 > 0) { //fB and fB2 on the same side
				if(posB_B2 > 0) Intersection = !Intersection;
			} else if(posB_A2 * posB2_A2 < 0) { //fB and fB2 on opposite side
				if(posA2_B < 0) Intersection = !Intersection;
			} else if(posB2_A2 == 0) { //fA2 and fB2 coplanar
				if(posA2_B * posB2_B < 0) {
					if(posB_B2 > 0) Intersection = !Intersection;
				} else {
					if(nA2nB2 < 0) {
						if(m_BOOP == UNION) Intersection = !Intersection;
					} else {
						if(m_BOOP == MINUS) Intersection = !Intersection;
					}
				}
			} else { //fA2 and fB coplanar
				if(posA2_B2 * posB_B2 < 0) {
					if(posB_B2 > 0) Intersection = !Intersection;
				} else {
					if(nA2nB < 0) {
						if(m_BOOP == UNION) Intersection = !Intersection;
					} else {
						if(m_BOOP == MINUS) Intersection = !Intersection;
					}
				}
			}

			//if Intersection == false, fA and fA2 are both inside or outside the second polyhedron.
			if(!Intersection) stop = true;

			//the intersection between (fA, fB2), (fA2, fB) and (fA2, fB2) are the same so these couples are removed from the list
			rmCouple(fA->Label, fB2->Label);
			rmCouple(fA2->Label, fB->Label);
			rmCouple(fA2->Label, fB2->Label);

			//it is possible that the direction of the intersection have to be inverted
			if(posB_A * posA2_A > 0 && posB_A * posB2_A >= 0 && posB2_B * posA_B > 0) invert_direction = true;
		}

		//if the intersection must not be compute
		if(stop) return;

		Info_Inter inter[4];
		inter[0].f = fA;
		inter[1].f = fA;
		inter[2].f = fB;
		inter[3].f = fB;

		//the two intersection points between the edges of a triangle and the
		//other triangle are computed for the two triangles
		switch(posBbin) {
		//common intersections : one point one one side of the plane and the two other points on the other side
		case 26:
		case 37:
			inter[0].he = heB[0];
			InterTriangleSegment(&inter[0]);
			inter[1].he = heB[1];
			InterTriangleSegment(&inter[1]);
			break;
		case 25:
		case 38:
			inter[0].he = heB[1];
			InterTriangleSegment(&inter[0]);
			inter[1].he = heB[2];
			InterTriangleSegment(&inter[1]);
			break;
		case 22:
		case 41:
			inter[0].he = heB[2];
			InterTriangleSegment(&inter[0]);
			inter[1].he = heB[0];
			InterTriangleSegment(&inter[1]);
			break;
		//particular cases : one point on the plane, one point one one side and one point on the other side
		case 6:
		case 9:
			inter[0].he = heB[2];
			InterTriangleSegment(&inter[0]);
			inter[1].he = heB[0];
			IsInTriangle(&inter[1]);
			break;
		case 18:
		case 33:
			inter[0].he = heB[0];
			InterTriangleSegment(&inter[0]);
			inter[1].he = heB[1];
			IsInTriangle(&inter[1]);
			break;
		case 24:
		case 36:
			inter[0].he = heB[1];
			InterTriangleSegment(&inter[0]);
			inter[1].he = heB[2];
			IsInTriangle(&inter[1]);
			break;
		//particular case : two points on the plane
		case 1:
		case 2:
			inter[0].he = heB[0];
			IsInTriangle(&inter[0]);
			inter[1].he = heB[2]->opposite();
			IsInTriangle(&inter[1]);
			break;
		case 16:
		case 32:
			inter[0].he = heB[1];
			IsInTriangle(&inter[0]);
			inter[1].he = heB[0]->opposite();
			IsInTriangle(&inter[1]);
			break;
		case 4:
		case 8:
			inter[0].he = heB[2];
			IsInTriangle(&inter[0]);
			inter[1].he = heB[1]->opposite();
			IsInTriangle(&inter[1]);
			break;
		default:
			return;
		}

		switch(posAbin) {
		//common intersections : one point one one side of the plane and the two other points on the other side
		case 26:
		case 37:
			inter[2].he = heA[0];
			InterTriangleSegment(&inter[2]);
			inter[3].he = heA[1];
			InterTriangleSegment(&inter[3]);
			break;
		case 25:
		case 38:
			inter[2].he = heA[1];
			InterTriangleSegment(&inter[2]);
			inter[3].he = heA[2];
			InterTriangleSegment(&inter[3]);
			break;
		case 22:
		case 41:
			inter[2].he = heA[2];
			InterTriangleSegment(&inter[2]);
			inter[3].he = heA[0];
			InterTriangleSegment(&inter[3]);
			break;
		//particular cases : one point on the plane, one point one one side and one point on the other side
		case 6:
		case 9:
			inter[2].he = heA[2];
			InterTriangleSegment(&inter[2]);
			inter[3].he = heA[0];
			IsInTriangle(&inter[3]);
			break;
		case 18:
		case 33:
			inter[2].he = heA[0];
			InterTriangleSegment(&inter[2]);
			inter[3].he = heA[1];
			IsInTriangle(&inter[3]);
			break;
		case 24:
		case 36:
			inter[2].he = heA[1];
			InterTriangleSegment(&inter[2]);
			inter[3].he = heA[2];
			IsInTriangle(&inter[3]);
			break;
		//particular case : two points on the plane
		case 1:
		case 2:
			inter[2].he = heA[0];
			IsInTriangle(&inter[2]);
			inter[3].he = heA[2]->opposite();
			IsInTriangle(&inter[3]);
			break;
		case 16:
		case 32:
			inter[2].he = heA[1];
			IsInTriangle(&inter[2]);
			inter[3].he = heA[0]->opposite();
			IsInTriangle(&inter[3]);
			break;
		case 4:
		case 8:
			inter[2].he = heA[2];
			IsInTriangle(&inter[2]);
			inter[3].he = heA[1]->opposite();
			IsInTriangle(&inter[3]);
			break;
		default:
			return;
		}

		//if two distincts points belongs to the two triangles
		if(IsSegment(inter)) {
			//we get this segment in ptInter
			std::vector<InterId> ptInter;
			Get_Segment(inter, ptInter);
			//and we build the opposite segment in ptInterInv
			std::vector<InterId> ptInterInv;
			ptInterInv.push_back(ptInter[1]);
			ptInterInv.push_back(ptInter[0]);

			//the segments are stored in the concerned triangles, and oriented
			if(CGAL::cross_product(nA, nB) * (InterPts[ptInter[1]] - InterPts[ptInter[0]]) * ((invert_direction == true)?-1:1) > 0) {
				switch(m_BOOP) {
				case UNION:
					Inter_tri[fA->Label].CutList.push_back(ptInter);
					if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInter);
					Inter_tri[fB->Label].CutList.push_back(ptInterInv);
					if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInterInv);
					break;
				case INTER:
					Inter_tri[fA->Label].CutList.push_back(ptInterInv);
					if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInterInv);
					Inter_tri[fB->Label].CutList.push_back(ptInter);
					if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInter);
					break;
				case MINUS:
					Inter_tri[fA->Label].CutList.push_back(ptInter);
					if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInter);
					Inter_tri[fB->Label].CutList.push_back(ptInter);
					if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInter);
					break;
				}
			} else {
				switch(m_BOOP) {
				case UNION:
					Inter_tri[fA->Label].CutList.push_back(ptInterInv);
					if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInterInv);
					Inter_tri[fB->Label].CutList.push_back(ptInter);
					if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInter);
					break;
				case INTER:
					Inter_tri[fA->Label].CutList.push_back(ptInter);
					if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInter);
					Inter_tri[fB->Label].CutList.push_back(ptInterInv);
					if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInterInv);
					break;
				case MINUS:
					Inter_tri[fA->Label].CutList.push_back(ptInterInv);
					if(edgeA != 3) Inter_tri[fA2->Label].CutList.push_back(ptInterInv);
					Inter_tri[fB->Label].CutList.push_back(ptInterInv);
					if(edgeB != 3) Inter_tri[fB2->Label].CutList.push_back(ptInterInv);
					break;
				}
			}
		}
	}

	/*! \brief Compute the intersection between a facet and a halfedge
	 * \param inter : A pointer to an Info_Inter structure.*/
	void InterTriangleSegment(Info_Inter* inter)
	{
		Facet_handle f = inter->f;
		Halfedge_handle he = inter->he;
		//if the intersection has been computed, the function returns directly the Id of the intersection
		if(Inter_tri[f->Label].RefInter.count(he->Label) != 0) {
			inter->Id = Inter_tri[f->Label].RefInter[he->Label];
			return;
		}
		//else, the calculation is done

		//this method is called when the intersection is not on the vertex pointed by the halfedge
		inter->IsOnVertex = false;
		//the intersection does not have an Id. 0xFFFFFFFF is set (this value means "no Id")
		inter->Id = 0xFFFFFFFF;

		Vector_3 e1, e2, dir, p, s, q;
		FT u, v, tmp;

		Point_3 s1 = he->opposite()->vertex()->point();
		Point_3 s2 = he->vertex()->point();
		Point_3 v0 = f->facet_begin()->vertex()->point();
		Point_3 v1 = f->facet_begin()->next()->vertex()->point();
		Point_3 v2 = f->facet_begin()->next()->next()->vertex()->point();

		//computation of the intersection (exact numbers)
		e1 = v1 - v0;
		e2 = v2 - v0;
		dir = s2 - s1;
		p = CGAL::cross_product(dir, e2);
		tmp = (FT)1/(p*e1);
		s = s1 - v0;
		u = tmp * s * p;
		if(u < 0 || u > 1) {
			//the intersection is not in the triangle
			inter->res = 7;
			return;
		}
		q = CGAL::cross_product(s, e1);
		v = tmp * dir * q;
		if(v < 0 || v > 1) {
			//the intersection is not in the triangle
			inter->res = 7;
			return;
		}
		if(u + v > 1) {
			//the intersection is not in the triangle
			inter->res = 7;
			return;
		}

		//the result is stored in inter->pt
		inter->pt = s1+(tmp*e2*q)*dir;

		//creation of the code for the location of the intersection
		inter->res = 0;
		if(u == 0) inter->res += 1;	//intersection on he(0)
		if(v == 0) inter->res += 2;	//intersection on he(1)
		if(u+v == 1) inter->res += 4;	//intersection on he(2)
	}


	/*! \brief Finds the position of a point in a 3d triangle
	 * \param inter : A pointer to an Info_Inter structure*/
	void IsInTriangle(Info_Inter* inter)
	{
		Facet_handle f = inter->f;
		Halfedge_handle he = inter->he;
		//if the intersection has been computed, the function returns directly the Id of the intersection
		if(Inter_tri[f->Label].RefInter.count(he->Label) != 0) {
			inter->Id = Inter_tri[f->Label].RefInter[he->Label];
			return;
		}
		//else, the calculation is done

		//this method is called when the intersection is exactly on the vertex pointed by the halfedge
		inter->IsOnVertex = true;
		//the intersection does not have an Id. 0xFFFFFFFF is set (this value means "no Id")
		inter->Id = 0xFFFFFFFF;

		Point_3 p = he->vertex()->point();
		Point_3 v0 = f->facet_begin()->vertex()->point();
		Point_3 v1 = f->facet_begin()->next()->vertex()->point();
		Point_3 v2 = f->facet_begin()->next()->next()->vertex()->point();

		Vector_3 N = Inter_tri[f->Label].norm_dir;
		FT u, v, w;

		u = N * CGAL::cross_product(v0 - v2, p - v2);
		if(u < 0) {
			//the intersection is not in the triangle
			inter->res = 7;
			return;
		}
		v = N * CGAL::cross_product(v1 - v0, p - v0);
		if(v < 0) {
			//the intersection is not in the triangle
			inter->res = 7;
			return;
		}
		w = N * CGAL::cross_product(v2 - v1, p - v1);
		if(w < 0) {
			//the intersection is not in the triangle
			inter->res = 7;
			return;
		}

		//the point is in the triangle
		inter->pt = p;

		//creation of the code for the location of the intersection
		inter->res = 0;
		if(u == 0) inter->res += 1;	//intersection on he(0)
		if(v == 0) inter->res += 2;	//intersection on he(1)
		if(w == 0) inter->res += 4;	//intersection on he(2)
	}

	/*! \brief Verify that the intersection is a segment
	 * \param inter : A pointer to four Info_Inter structures
	 * \return true if two distinct points are found in the four intersections computed*/
	bool IsSegment(Info_Inter* inter)
	{
		bool point = false; //true if a point is founded
		Point_3 pt; //the point founded
		bool id = false; //true if an Id is founded
		unsigned long Id = 0; //the Id founded // MT

		//each intersection is checked separately.
		//first intersection
		if(inter[0].Id != 0xFFFFFFFF) {
			//an Id different than 0xFFFFFFFF is founded
			//this intersection has already been computed and is valid
			id = true;
			Id = inter[0].Id;
		} else if(inter[0].res != 7) {
			//the intersection have no Id (0xFFFFFFFF)
			//but the intersection is in the triangle
			point = true;
			pt = inter[0].pt;
		}
		//second intersection
		if(inter[1].Id != 0xFFFFFFFF) {
			//an Id different than 0xFFFFFFFF is founded
			//this intersection has already been computed and is valid

			//if a point or an Id has already be founded, we founded the two distinct valid points (the intersection is a segment)
			//(it is not possible that the two first points are the same)
			if(point || id) return true;
			id = true;
			Id = inter[1].Id;
		} else if(inter[1].res != 7) {
			//the intersection have no Id (0xFFFFFFFF)
			//but the intersection is in the triangle

			//if a point or an Id has already be founded, we founded the two distinct valid points (the intersection is a segment)
			//(it is not possible that the two first points are the same)
			if(point || id) return true;
			point = true;
			pt = inter[1].pt;
		}
		//third intersection
		if(inter[2].Id != 0xFFFFFFFF) {
			//an Id different than 0xFFFFFFFF is founded
			//this intersection has already been computed and is valid

			//if a point or a different Id has already be founded, we founded the two distinct valid points (the intersection is a segment)
			//(it is not possible that the two first points are the same)
			if(point || (id && Id != inter[2].Id)) return true;
			id = true;
			Id = inter[2].Id;
		} else if(inter[2].res != 7) {
			//the intersection have no Id (0xFFFFFFFF)
			//but the intersection is in the triangle

			//if an Id or a different point has already be founded, we founded the two distinct valid points (the intersection is a segment)
			//(it is not possible that the two first points are the same)
			if((point && pt != inter[2].pt) || id) return true;
			point = true;
			pt = inter[2].pt;
		}
		//fourth intersection
		if(inter[3].Id != 0xFFFFFFFF) {
			//an Id different than 0xFFFFFFFF is founded
			//this intersection has already been computed and is valid

			//if a point or a different Id has already be founded, we founded the two distinct valid points (the intersection is a segment)
			//(it is not possible that the two first points are the same)
			if(point || (id && Id != inter[3].Id)) return true;
		} else if(inter[3].res != 7) {
			//the intersection have no Id (0xFFFFFFFF)
			//but the intersection is in the triangle

			//if an Id or a different point has already be founded, we founded the two distinct valid points (the intersection is a segment)
			//(it is not possible that the two first points are the same)
			if((point && pt != inter[3].pt) || id) return true;
		}
		return false;
	}

	/*! \brief Extracts the segment from a set of four intersection points and store these points in the list of intersecion points
	 * \n There must be two valid and distinct points in the set
	 * \param inter : A pointer to four Info_Inter structure
	 * \param I : A vector to store the Id of the two intersection points (the output segment)*/
	void Get_Segment(Info_Inter* inter, std::vector<InterId>& I)
	{
		for(unsigned int i = 0; i != 4; ++i) {
			//if the point have an Id
			if(inter[i].Id != 0xFFFFFFFF) {
				//the Id is stored if it is not already done
				if(I.size() == 0 || I[0] != inter[i].Id) I.push_back(inter[i].Id);
			}
			//else if the point is valid
			else if(inter[i].res != 7) {
				//the intersection point is stored in the list of the intersection points
				//and its new Id is stored in the output segment
				if(I.size() == 0 || InterPts[I[0]] != inter[i].pt) {
					Store_Intersection(&inter[i]);
					I.push_back(inter[i].Id);
				}
			}
			//return if the two points are founded
			if(I.size() == 2) return;
		}
	}

	/*! \brief Store the intersection and memorize it for every couples of facet-halfedge
	 * \param inter : A pointer to an Info_Inter structure*/
	void Store_Intersection(Info_Inter* inter)
	{
		Facet_handle f;
		Halfedge_handle he;
		f = inter->f;
		he = inter->he;
		InterId I;

		//store the point to the list of the intersections and store its new Id
		inter->Id = InterPts.size();
		I = inter->Id;
		InterPts.push_back(inter->pt);

		//add this point as a vertex of the result
		ppbuilder.add_vertex(inter->pt, inter->Id);

		//if the intersection is on the vertex pointed by the halfedge (he), we update the Id (Label) of this vertex
		if(inter->IsOnVertex) he->vertex()->Label = I;

		//the intersection is memorized for each possible couple of (facet, halfedge) concerned by the intersection
		//if the intersection is exactly on the vertex pointed by the halfedge (he), it is necessary to take account
		//of every halfedge pointing to this vertex
		switch(inter->res) {
		case 0: { //intersection on the facet
			if(!inter->IsOnVertex) {
				Inter_tri[f->Label].RefInter[he->Label] = I;
				Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
			} else {
				Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
				do {
					Inter_tri[f->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			}
		}
		break;
		case 1: { //Intersection on the first halfedge of the facet
			Inter_tri[f->Label].PtList.insert(I);
			Inter_tri[f->facet_begin()->opposite()->facet()->Label].PtList.insert(I);
			if(!inter->IsOnVertex) {
				Inter_tri[he->facet()->Label].PtList.insert(I);
				Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);
				Inter_tri[f->Label].RefInter[he->Label] = I;
				Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
				Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[he->Label] = I;
				Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[he->opposite()->Label] = I;
				Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->Label] = I;
				Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->opposite()->Label] = I;
				Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->Label] = I;
				Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->opposite()->Label] = I;
			} else {
				Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
				do {
					Inter_tri[f->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->facet_begin()->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			}
		}
		break;
		case 2: { //Intersection on the second halfedge of the facet
			Inter_tri[f->Label].PtList.insert(I);
			Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].PtList.insert(I);
			if(!inter->IsOnVertex) {
				Inter_tri[he->facet()->Label].PtList.insert(I);
				Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);
				Inter_tri[f->Label].RefInter[he->Label] = I;
				Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
				Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[he->Label] = I;
				Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[he->opposite()->Label] = I;
				Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->Label] = I;
				Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->opposite()->Label] = I;
				Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->Label] = I;
				Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->opposite()->Label] = I;
			} else {
				Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
				do {
					Inter_tri[f->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->facet_begin()->next()->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			}
		}
		break;
		case 3: { //Intersection on the first and second halfedge of the facet (vertex pointed by the first halfedge)
			//update the Id (Label) of the first vertex of the facet
			f->facet_begin()->vertex()->Label = I;
			if(!inter->IsOnVertex) {
				Inter_tri[he->facet()->Label].PtList.insert(I);
				Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);

				Halfedge_around_vertex_circulator	H_circ = f->facet_begin()->vertex_begin(),
													H_end = f->facet_begin()->vertex_begin();
				do {
					Inter_tri[H_circ->facet()->Label].RefInter[he->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[he->opposite()->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			} else {
				Halfedge_around_vertex_circulator	H_circ = he->vertex_begin(),
													H_end = he->vertex_begin();
				do {
					Halfedge_around_vertex_circulator	F_circ = f->facet_begin()->vertex_begin(),
														F_end = f->facet_begin()->vertex_begin();
					do {
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->opposite()->Label] = I;
						F_circ++;
					} while(F_circ != F_end);
					H_circ++;
				} while(H_circ != H_end);
			}
		}
		break;
		case 4: { //Intersection on the third halfedge of the facet
			Inter_tri[f->Label].PtList.insert(I);
			Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].PtList.insert(I);
			if(!inter->IsOnVertex) {
				Inter_tri[he->facet()->Label].PtList.insert(I);
				Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);
				Inter_tri[f->Label].RefInter[he->Label] = I;
				Inter_tri[f->Label].RefInter[he->opposite()->Label] = I;
				Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[he->Label] = I;
				Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[he->opposite()->Label] = I;
				Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->next()->Label] = I;
				Inter_tri[he->facet()->Label].RefInter[f->facet_begin()->next()->next()->opposite()->Label] = I;
				Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->next()->Label] = I;
				Inter_tri[he->opposite()->facet()->Label].RefInter[f->facet_begin()->next()->next()->opposite()->Label] = I;
			} else {
				Halfedge_around_vertex_circulator H_circ = he->vertex_begin(), H_end = he->vertex_begin();
				do {
					Inter_tri[f->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[f->facet_begin()->next()->next()->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->next()->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[f->facet_begin()->next()->next()->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			}
		}
		break;
		case 5: { //Intersection on the first and third halfedge of the facet (vertex pointed by the third halfedge)
			//update the Id (Label) of the third vertex of the facet
			f->facet_begin()->next()->next()->vertex()->Label = I;
			if(!inter->IsOnVertex) {
				Inter_tri[he->facet()->Label].PtList.insert(I);
				Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);

				Halfedge_around_vertex_circulator	H_circ = f->facet_begin()->next()->next()->vertex_begin(),
													H_end = f->facet_begin()->next()->next()->vertex_begin();
				do {
					Inter_tri[H_circ->facet()->Label].RefInter[he->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[he->opposite()->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			} else {
				Halfedge_around_vertex_circulator 	H_circ = he->vertex_begin(),
													H_end = he->vertex_begin();
				do {
					Halfedge_around_vertex_circulator 	F_circ = f->facet_begin()->next()->next()->vertex_begin(),
														F_end = f->facet_begin()->next()->next()->vertex_begin();
					do {
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->opposite()->Label] = I;
						F_circ++;
					} while(F_circ != F_end);
					H_circ++;
				} while(H_circ != H_end);
			}
		}
		break;
		case 6: { //Intersection on the second and third halfedge of the facet (vertex pointed by the second halfedge)
			//update the Id (Label) of the second vertex of the facet
			f->facet_begin()->next()->vertex()->Label = I;
			if(!inter->IsOnVertex) {
				Inter_tri[he->facet()->Label].PtList.insert(I);
				Inter_tri[he->opposite()->facet()->Label].PtList.insert(I);

				Halfedge_around_vertex_circulator	H_circ = f->facet_begin()->next()->vertex_begin(),
													H_end = f->facet_begin()->next()->vertex_begin();
				do {
					Inter_tri[H_circ->facet()->Label].RefInter[he->Label] = I;
					Inter_tri[H_circ->facet()->Label].RefInter[he->opposite()->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[he->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->Label] = I;
					Inter_tri[he->opposite()->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
					H_circ++;
				} while(H_circ != H_end);
			} else {
				Halfedge_around_vertex_circulator	H_circ = he->vertex_begin(),
													H_end = he->vertex_begin();
				do {
					Halfedge_around_vertex_circulator	F_circ = f->facet_begin()->next()->vertex_begin(),
														F_end = f->facet_begin()->next()->vertex_begin();
					do {
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->Label] = I;
						Inter_tri[F_circ->facet()->Label].RefInter[H_circ->opposite()->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->Label] = I;
						Inter_tri[H_circ->facet()->Label].RefInter[F_circ->opposite()->Label] = I;
						F_circ++;
					} while(F_circ != F_end);
					H_circ++;
				} while(H_circ != H_end);
			}
		}
		break;
		}

	}

	/*! \brief Add a facet to the result
	 * \param pFacet : A handle to the facet to add
	 * \param facet_from_A : must be true if the facet belongs to the first polyhedron*/
	void add_facet_to_solution(Facet_handle& pFacet, bool facet_from_A)
	{
		//if the facet contains an intersection point but no intersection segment, the facet must be triangulate before
		if(pFacet->Label < Inter_tri.size()) {
			Triangle_Cut TriCut = Inter_tri[pFacet->Label];
			Halfedge_handle he = pFacet->facet_begin();
			//creation of the triangulation
			Triangulation<Kernel> T(he, TriCut.norm_dir);
			//add the intersection points to the triangulation
			for(std::set<InterId>::iterator i = TriCut.PtList.begin(); i != TriCut.PtList.end(); ++i) {
				T.add_new_pt(InterPts[*i], (unsigned long&)*i);     // MT: ajout cast
			}
			//get all the triangles of the triangulation
			std::vector<std::vector<unsigned long> > Tri_set = T.get_all_triangles((m_BOOP == MINUS && !TriCut.Facet_from_A)?true:false);
			//add these triangles to the result
			ppbuilder.add_triangle(Tri_set, he);
		} else {
			//the facet is added to the result. If the facet belongs to the second polyhedron, and if the
			//Boolean operation is a Subtraction, it is necessary to invert the orientation of the facet.
			if(m_BOOP == MINUS && !facet_from_A) ppbuilder.add_triangle(pFacet, true);
			else ppbuilder.add_triangle(pFacet, false);
		}
		//the tag of the three neighboring facets is updated
		pFacet->facet_begin()->opposite()->facet()->IsExt = true;
		pFacet->facet_begin()->next()->opposite()->facet()->IsExt = true;
		pFacet->facet_begin()->next()->next()->opposite()->facet()->IsExt = true;
	}

#ifdef BOOLEAN_OPERATIONS_DEBUG

	/*! \brief Colors the facets of the input polyhedra
	 * \n The intersected facets are red
	 * \n The facets that belong to the result are in green
	 * \n The facets that does not belong to the result are in blue*/
	void ColorType()
	{
		Facet_iterator pFacet =	NULL;
		for(pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); pFacet++) {
			if(pFacet->Label < Inter_tri.size()) pFacet->color(1.0, 0.0 ,0.0);
			else if(pFacet->IsExt) pFacet->color(0.0, 1.0 ,0.0);
			else pFacet->color(0.0, 0.0 ,1.0);
		}
		for(pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); pFacet++) {
			if(pFacet->Label < Inter_tri.size()) pFacet->color(1.0, 0.0 ,0.0);
			else if(pFacet->IsExt) pFacet->color(0.0, 1.0 ,0.0);
			else pFacet->color(0.0, 0.0 ,1.0);
		}
	}


	/*! \brief Writes a report containing the computation time of the diffrent parts of the algorithm
	 * \param pMout : The result polyhedron*/
	void WriteData(PolyhedronPtr& pMout)
	{
		unsigned int N_IFA = 0;
		unsigned int N_FFA = 0;
		unsigned int N_LFFA = 0;
		unsigned int N_IFB = 0;
		unsigned int N_FFB = 0;
		unsigned int N_LFFB = 0;
		unsigned int N_CF;

		for(Facet_iterator pFacet = m_pA->facets_begin(); pFacet != m_pA->facets_end(); ++pFacet) {
			if(pFacet->Label < Inter_tri.size()) N_IFA++;
			else if(pFacet->IsExt) N_FFA++;
			else N_LFFA++;
		}
		for(Facet_iterator pFacet = m_pB->facets_begin(); pFacet != m_pB->facets_end(); ++pFacet) {
			if(pFacet->Label < Inter_tri.size()) N_IFB++;
			else if(pFacet->IsExt) N_FFB++;
			else N_LFFB++;
		}
		N_CF = pMout->size_of_facets() - N_FFA - N_FFB;

		ofstrtime << "Computation time :"																<< std::endl;
		ofstrtime																						<< std::endl;
		ofstrtime << "Algorithm"																		<< std::endl;
		ofstrtime << "   Initialization :             "			<< tr(duration_Init)					<< std::endl;
		ofstrtime << " + Finding the Intersections :  "			<< tr(duration_FindCouples)				<< std::endl;
		ofstrtime << " + Compute the Intersections :  "			<< tr(duration_ComputeIntersections)	<< std::endl;
		ofstrtime << " + Cut the Intersected Facets : "			<< tr(duration_CutIntersectedFacets)	<< std::endl;
		ofstrtime << " + Complete the result :        "			<< tr(duration_PropagateFacets)			<< std::endl;
		ofstrtime << " + Create the polyhedron :      "			<< tr(duration_delegate)				<< std::endl;
		ofstrtime << "---------------------------------------"											<< std::endl;
		ofstrtime << " Total :                        "			<< tr(duration_total)					<< std::endl;
		ofstrtime																						<< std::endl;
		ofstrtime																						<< std::endl;
		ofstrtime << "Details :"																		<< std::endl;
		ofstrtime																						<< std::endl;
		ofstrtime << "Polyedron A :"																	<< std::endl;
		ofstrtime << "Number of Facets :                   "	<< m_pA->size_of_facets()				<< std::endl;
		ofstrtime << "Number of Intersected Facets :       "	<< N_IFA								<< std::endl;
		ofstrtime << "Number of Facets not in the result : "	<< N_LFFA								<< std::endl;
		ofstrtime																						<< std::endl;
		ofstrtime << "Polyedron B :"																	<< std::endl;
		ofstrtime << "Number of Facets :                   "	<< m_pB->size_of_facets()				<< std::endl;
		ofstrtime << "Number of Intersected Facets :       "	<< N_IFB								<< std::endl;
		ofstrtime << "Number of Facets not in the result : "	<< N_LFFB								<< std::endl;
		ofstrtime																						<< std::endl;
		ofstrtime << "Result :"																			<< std::endl;
		ofstrtime << "Number of Facets :                   "	<< pMout->size_of_facets()				<< std::endl;
		ofstrtime << "Number of Facets from A :            "	<< N_FFA								<< std::endl;
		ofstrtime << "Number of Facets from B :            "	<< N_FFB								<< std::endl;
		ofstrtime << "Number of Created Facets :           "	<< N_CF									<< std::endl;
	}

#endif // BOOLEAN_OPERATIONS_DEBUG

	/*! \brief Boolean operation computed*/
	Bool_Op m_BOOP;
	/*! \brief The first input polyhedron*/
	MEPP_Polyhedron* m_pA;
	/*! \brief The second input polyhedron*/
	MEPP_Polyhedron* m_pB;
	/*! \brief The polyhedron builder*/
	CPolyhedron_from_polygon_builder_3<HDS> ppbuilder;

	/*! \brief Lists the couples of facets that intersect*/
	std::map<FacetId, std::set<FacetId> > m_Couples;
	/*! \brief Lists the exact intersection points computed*/
	std::vector<Point_3> InterPts;
	/*! \brief Informations about the intersected facets*/
	std::vector<Triangle_Cut> Inter_tri;
	/*! \brief Index to obtain the handle of a facet with its Id*/
	std::vector<Facet_handle> Facet_Handle;

	/*! \brief the AABB-tree*/
	AABB_Tree tree;

#ifdef BOOLEAN_OPERATIONS_DEBUG

	/*! \brief The output report file*/
	std::ofstream ofstrtime;
	/*! \brief Initialisation time*/
	double duration_Init;
	/*! \brief Time to find all the couples of facet that intersect*/
	double duration_FindCouples;
	/*! \brief Time to Compute all the intersections*/
	double duration_ComputeIntersections;
	/*! \brief Time to Cut the facets*/
	double duration_CutIntersectedFacets;
	/*! \brief Time to complete the result with the facets from the input polyhedra*/
	double duration_PropagateFacets;
	/*! \brief Time to create the result polyhedron*/
	double duration_delegate;
	/*! \brief Time to Compute a Boolean operation*/
	double duration_total;

#endif // BOOLEAN_OPERATIONS_DEBUG
};

#endif // BOOLEAN_OPERATIONS_H
