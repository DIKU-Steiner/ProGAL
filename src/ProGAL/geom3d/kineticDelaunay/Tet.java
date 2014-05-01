package ProGAL.geom3d.kineticDelaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import ProGAL.dataStructures.Queue;
import ProGAL.geom2d.TriangulationVertex;
import ProGAL.geom3d.Plane;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.kineticDelaunay.Hole.Face;
import ProGAL.geom3d.kineticDelaunay.KineticDelaunayTessellation.ProblemInstanceType;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;

public class Tet {
	private Vertex[] corners = new Vertex[4];
	private Set<Edge> edges = new LinkedHashSet<Edge>(); // Daisy
	private Tri[] tris = new Tri[4]; // Daisy
	Tet[] neighbors = new Tet[4];
	Sphere circumSphere = null;
	Integer count = null;
	boolean dAlive = true;
	boolean cAlive = true;
	Shape[] LSSs = new Shape[6];
	Shape[] faces = new Shape[4];
	Shape faceShape = null;
	boolean onStack = false;   // used when creating holes
	Face selectedFace = null;  // used when creating holes
	boolean[] oppositeInside = new boolean[4];
	boolean centerInside;
	boolean flag = false;
	Integer depth = null;
	int alph = 0;
	
	public Tet(Vertex[] corners){
		this.corners = corners;
		this.sortCorners();
		normalizePredicates();
	}
	
	public Tet(Vertex v0, Vertex v1, Vertex v2, Vertex v3) {
		corners[0] = v0;
		corners[1] = v1;
		corners[2] = v2;
		corners[3] = v3;
		this.sortCorners();
		normalizePredicates();
	}
	 
	public Tet(Tetrahedron tetra) {
		for (int i = 0; i < 4; i++) corners[i] = new Vertex(tetra.getCorner(i));
		this.sortCorners();		
		normalizePredicates();
	}

	public Tet clone() {
		Tet nt = new Tet(corners);
		for (int i = 0; i < 4; i++) nt.neighbors[i] = neighbors[i];
		return nt;
	}
	//Daisy
	public void setEdge(Edge e) {
		edges.add(e); 
	}
	public void setEdges(Edge[] es) {
		for (Edge e : es) {
			edges.add(e);
		}
	}
	public void removeEdge(Edge e) {
		edges.remove(e);
	}
	public void setTri(Tri t, int i) {
		tris[i] = t; 
	}
	public Tri getTri(int i) {
		return tris[i];
	}
	public Tri getTri(Vertex v) {
		for (int i = 0 ; i<4 ; i++) {
			if (corners[i].equals(v)) {
				return tris[i];
			}
		}
		return null;
	}
	public Tri[] getTris() {
		return tris;
	}
	
	public Tetrahedron makeTetrahedron() {
		return new Tetrahedron(corners[0], corners[1], corners[2], corners[3]);
	}
	
	public Vertex getCorner(int i) { return corners[i]; }
	
	
	public boolean getFlag() { return flag; }
	
	public void setFlag(boolean flag) { this.flag = flag; }
	
	/** Returns the pair of vertices of this tetrahedron other than a and b*/
	public Vertex[] getCornerPair(Vertex a, Vertex b) {
		Vertex[] pair = new Vertex[2];
		int k;
		for (k = 0; k < 4; k++) {	
			pair[0] = corners[k];
			if ((pair[0] != a) && (pair[0] != b)) break;
		}
		for (int l = k+1; l < 4; l++) {
			pair[1] = corners[l];
			if ((pair[1] != a) && (pair[1] != b)) break;
		}
		return pair;
	}
	

	/** Returns the sphere circumscribing this tetrahedron */
	public Sphere getCircumSphere() { 
//		if (circumSphere == null) {
			setCircumSphere();
			return circumSphere;
//		}
//		return circumSphere;
	}
	
	/** Returns the radius of the sphere circumscribing this tetrahedron */
	public double getCircumSphereRadius() { return getCircumSphere().getRadius(); }
	
	public void setCircumSphere() {
		Tetrahedron tetra = makeTetrahedron();
		circumSphere = tetra.circumSphere();
	}
	
	/** Calculate the radius of the circumsphere. */
	public double circumRadius(){
		Vector a = corners[3].vectorTo(corners[0]);
		Vector b = corners[3].vectorTo(corners[1]);
		Vector c = corners[3].vectorTo(corners[2]);
		Vector O = b.cross(c).multiplyThis(a.dot(a));
		O.addThis(c.cross(a).multiplyThis(b.dot(b)));
		O.addThis(a.cross(b).multiplyThis(c.dot(c)));
		O.multiplyThis(1.0/(2*a.dot(b.crossThis(c))));
		return O.length();
	}

	// Daisy
	public void setAlph(int i) {
		alph = i;
	}
	public int getAlph() {
		return alph;
	}
	
	public boolean hasVertex(Vertex v) {
		for (int i = 0; i < 4; i++) 
			if (corners[i] == v) return true;
		return false;
	}
	//Daisy :
	public boolean hasEdge(Edge e) {
		if (edges.contains(e)) return true;
		return false;
	}
	public Edge[] getCommonEdges(Tet t) {
		Edge[] es = new Edge[4];
		int num = 0;
		for (Edge e : edges) {
			if (t.hasEdge(e)) {
				es[num] = e;
				num++;
			}
		}
		return es;
	}

	public boolean hasNeighbor(Tet t) {
		for (int i = 0; i < 4; i++) 
			if (neighbors[i] == t) return true;
		return false;
	}

	public boolean isAlive() { return dAlive; }
	public void setAlive(boolean dAlive) {this.dAlive = dAlive; }
	
	/** Returns TRUE if one of the corners is a big point - there are 4 big points with id's 0, 1, 2 , 3*/
	public boolean isBig() {
		return ((corners[0].getId() < 4) || (corners[1].getId() < 4) || (corners[2].getId() < 4) || (corners[3].getId() < 4));
	}
	
	public boolean isFlat() { return Point.coplanar(corners[0], corners[1], corners[2], corners[3]); }
	
	/** returns true if this tetrahedron and its neighbor tet form a convex set. */
	public boolean isConvex(Tet tet) {
		Triangle tr = commonFace(tet);
		Point p = corners[apex(tet)];
		Point q = tet.corners[tet.apex(this)];
		return tr.getIntersection(p, q) != null;
	}
	
	/** returns TRUE if this terahedron is in the alpha complex */
	public boolean isAlpha(double alpha) { return (getCircumSphereRadius() < (alpha+Constants.EPSILON)); }
	
	private void normalizePredicates(){
		for(int i=0;i<4;i++){
			oppositeInside[i] = orient(i, corners[i]);
		}
		centerInside = inSphere(Point.getMidpoint(corners[0], corners[1]));
	}
	
	
	public int getCount() { 
		if (count != null) return count; 
		else {
			count = 0;
			for (int i = 0; i < 4; i++) { if (corners[i].getType() == Vertex.VertexType.R) count = count + (int)Math.pow(2,i);  }
			return count;
		}
	}
	
	//Daisy:
	public int  nrRotating() {
		int ret = 0;
		for (Vertex v : corners){
			if (v.getType()==Vertex.VertexType.R) ret += 1;
		}
		return ret;
	}
	
	/** Returns the vertex not in the facet-sharing tetrahedron tet */
	public Vertex getOppVertex(Tet tet) {
		for (int k = 0; k < 4; k++) {
			if (!hasVertex(tet.corners[k])) return tet.corners[k];
		}
		return null;
	}
	
	/** Returns face sharing tetrahedra */
	public List<Tet> getFaceSharingTetrahedra() {
		List<Tet> nList = new ArrayList<Tet>();
		for (int k = 0; k < 4; k++) {
			if (!neighbors[k].isBig())
			nList.add(neighbors[k]);
		}
		return nList;
	}
	
	//Daisy
	public Tet getNeighbour(int i) {
		return neighbors[i];
	}
	
	/** Returns all tetrahedra sharing the edge of this tetrahedron between vertices a and b */
	public void getEdgeSharingTetrahedra(Vertex a, Vertex b, HashSet<Tet> processedTets, HashSet<Vertex> processedVers) {
		Vertex[] pair = getCornerPair(a, b);
        getEdgeSharingTetrahedra(a, b, pair, processedTets, processedVers, this, this);
	}
	
	private void getEdgeSharingTetrahedra(Vertex a, Vertex b, Vertex[] pair, HashSet<Tet> processedTets, HashSet<Vertex> processedVers, Tet prev, Tet first) {
		if (!processedTets.contains(this) && !isBig()) processedTets.add(this);
		if (!processedVers.contains(pair[0])) processedVers.add(pair[0]);
		if (!processedVers.contains(pair[1])) processedVers.add(pair[1]);
		Tet tet = neighbors[indexOf(pair[0])];
		if ((tet != prev) && (tet != first)) tet.getEdgeSharingTetrahedra(a, b, tet.getCornerPair(a, b), processedTets, processedVers, this, first);
		else {
			tet = neighbors[indexOf(pair[1])];
			if ((tet != prev) && (tet != first)) tet.getEdgeSharingTetrahedra(a, b, tet.getCornerPair(a, b), processedTets, processedVers, this, first);
		}
	}
	
	
	/** Returns all tetrahedra sharing a vertex of this tetrahedron */
/*	public void getVertexSharingTetrahedra(double alpha, HashSet<Tet> processedTets, HashSet<Vertex> processedVers) {
		if (!processedTets.contains(this) && (getCircumSphereRadius() < alpha)) processedTets.add(this);
		for (int k = 0; k < 4; k++) getVertexSharingTetrahedra(corners[k], alpha, processedTets, processedVers);
	}
*/
	/** Returns all tetrahedra sharing a given vertex a of this tetrahedron */
	public void getVertexSharingTetrahedra(Vertex a, double alpha, HashSet<Tet> alphaTets, HashSet<Tet> processedTets, HashSet<Vertex>processedVers) {
		if (!processedTets.contains(this)) {
			processedTets.add(this);
			if (!isBig() && (getCircumSphereRadius() < alpha)) alphaTets.add(this);		
			Vertex b;
			int indxA = indexOf(a);
			for (int k = 0; k < 4; k++) {
				if (k != indxA) {
					b = corners[k];
				    if (!processedVers.contains(b) && !b.isBig() && (getCircumSphereRadius() < alpha) && (indxA < b.getId())) processedVers.add(b);
				    if (neighbors[k] != null)
				    	neighbors[k].getVertexSharingTetrahedra(a, alpha, alphaTets, processedTets, processedVers);
				}
			}
		}
	}
	
	
	private boolean orient(int face, Point p){
		Matrix m = new Matrix(4,4);
		for(int r=0;r<3;r++){
			for(int c=0;c<3;c++){
				m.set(r, c, corners[(r+(r>=face?1:0))%4].getCoord(c));
			}
			m.set(r, 3, 1);
		}
		for(int c=0;c<3;c++)
			m.set(3, c, p.getCoord(c));
		m.set(3,3,1);

		double det = m.determinant();
		if (Math.abs(det)<Constants.EPSILON) return true;
		return det<0;
	}
	
	// Daisy
	public Set<Edge> getEdges() { return edges; }
	public void killEdges() {
		for (Edge e : edges) {
			e.setAlive(false);
		}
	}
	public void reviveEdges() {
		for (Edge e : edges) {
			e.setAlive(true);
		}
	}
	
	private boolean inSphere(Point p){
		Matrix m = new Matrix(5,5);
		for(int r=0;r<4;r++){
			for(int c=0;c<3;c++){
				m.set(r, c, corners[r].getCoord(c));
			}
			m.set(r, 3, corners[r].dot(corners[r]));
			m.set(r, 4, 1);
		}
		for(int c=0;c<3;c++)
			m.set(4, c, p.getCoord(c));
		m.set(4, 3, p.dot(p));
		m.set(4, 4, 1);

		double det = m.determinant();
		if (Math.abs(det)<Constants.EPSILON) return true;
		return det<0;
	}
	
	public boolean insideFace(int face, Point p){
		return orient(face, p)==oppositeInside[face];
	}
	

	public boolean insideCircumsphere(Point p){
		return inSphere(p)==centerInside;
	}
	

	public int apex(Tet c){
		for(int i=0;i<4;i++) if(neighbors[i].equals(c)) return i;
		return -1;
	}
	
	public Triangle commonFace(Tet t) {
		Vertex[] common = new Vertex[3];
		int k = 0;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				if (corners[i] == t.corners[j]) {
					common[k++] = corners[i];
					continue;
				}
				if (k == 3) break;
			}
			if (k == 3) break;
		}
		return new Triangle(common);
	}
	
	public int indexOf(Vertex v) {
		for(int i = 0; i < 3; i++) if(corners[i].equals(v)) return i;
		return 3;
	}

	/* Returns index of vertex v */
	public int indexOf_slow(Vertex v){
		return Arrays.binarySearch(corners, v);
	}
	
	/* Returns index of fourth vertex */
	public int indexOf(Vertex a, Vertex b, Vertex c) {
		for (int i = 0; i < 4; i++)
			if ((!corners[i].equals(a)) && (!corners[i].equals(b)) && (!corners[i].equals(c))) return i;
		return -1;
	}
	

	private final void swap(int i, int j){
		Vertex tmpV = corners[i];
		Tet tmpT = neighbors[i];
		corners[i] = corners[j];
		neighbors[i] = neighbors[j];
		corners[j] = tmpV;
		neighbors[j] = tmpT;
	}
	
	
	/** 
	 * Sorts the corners while maintaining that corners[i] is not on the shared 
	 * face between <code>this</code> and neighbors[i]. Uses 5 comparisons which is optimal. 
	 */
	protected final void sortCorners(){
		if(corners[0].compareTo(corners[1])>0) swap(0,1);//First two elements are sorted
		if(corners[2].compareTo(corners[3])>0) swap(2,3);//Last two elements are sorted
		if(corners[0].compareTo(corners[2])>0) swap(0,2);//Smallest element is at 0
		if(corners[1].compareTo(corners[3])>0) swap(1,3);//Largest element is at 3
		if(corners[1].compareTo(corners[2])>0) swap(1,2);//Check middle elements
	}
	
	
	public String toString(){
		return Arrays.toString(corners);
	}
	
	public void toConsole() {
		System.out.println(toString());
	}
	
	public void toConsoleNeighbors() {
		System.out.print(toString() + " has neighbours: ");
		for (int i = 0; i < 4; i++ ) 
			if (neighbors[i] != null) {
				System.out.print(neighbors[i].toString()); 
				if (!neighbors[i].hasNeighbor(this)) System.out.print("?");
			}
			else System.out.print("[null]");
		System.out.println();
	}
	
	public void fromSceneEdges(J3DScene scene) {
		for (int i = 0; i < 6; i++) {
			scene.removeShape(LSSs[i]);
		}
		scene.repaint();
	}
	
	/* draws edges of the tetrahedron.*/
	public void toSceneEdges(J3DScene scene, Color clr, double width) {
		int k = 0;
		for (int i = 0; i < 3; i++)
			for (int j = i+1; j < 4; j++) {
				if (LSSs[k] == null) LSSs[k] = new LSS(corners[i], corners[j], width);
				scene.addShape(LSSs[k++], clr, 3);
			}
		scene.repaint();
	}

	/* draws edges of the tetrahedron. Special width is used if some corners are big points.*/
	public void toSceneEdges(J3DScene scene, Color clr, double width, double bigWidth) {
		double edgeWidth = width;
		for (int i = 0; i < 4; i++)  if (corners[i].getId() < 4) edgeWidth = bigWidth;
		int k = 0;
		for (int i = 0; i < 3; i++)
			for (int j = i+1; j < 4; j++) {
				if (LSSs[k] == null) LSSs[k] = new LSS(corners[i], corners[j], edgeWidth);
				scene.addShape(LSSs[k++], clr, 3);
			}
		scene.repaint();
	}

	
	public void toSceneFaces(J3DScene scene, Color clr) {
		if (faces[0] == null) faces[0] = new Triangle(corners[1], corners[2], corners[3]);
		if (faces[1] == null) faces[1] = new Triangle(corners[0], corners[2], corners[3]);
		if (faces[2] == null) faces[2] = new Triangle(corners[0], corners[1], corners[3]);
		if (faces[3] == null) faces[3] = new Triangle(corners[0], corners[1], corners[2]);
		for (int i = 0; i < 4; i++) scene.addShape(faces[i], clr, 1);
	}
	
	public Shape toSceneFace(J3DScene scene, int i, Color clr) {
		if (faces[i] == null) faces[i] = new Triangle(corners[(i+1)%4], corners[(i+2)%4], corners[(i+3)%4]);
		scene.addShape(faces[i], clr, 1);
		return faces[i];
	}
	
	public void fromSceneFaces(J3DScene scene) {
		for (int i = 0; i < 4; i++) {
			scene.removeShape(faces[i]);
		}
		scene.repaint();
	}
	
	public void fromSceneFace(J3DScene scene, int i) {
		scene.removeShape(faces[i]);
		scene.repaint();
	}

	public ArrayList<Tet> breadthFirstTetrahedra(int maxDepth) {
		Queue queue = new Queue();
		ArrayList<Tet> tets = new ArrayList<Tet>();
		depth = 0;
		queue.push(this);
		
		while (!queue.isEmpty()) {
			Tet tet = (Tet)queue.pop();
			tets.add(tet);
			for (int k = 0; k < 4; k++) {
				Tet nTet = tet.neighbors[k];
				if (nTet.depth == null) {
					nTet.depth = tet.depth + 1;
					if (nTet.depth < maxDepth) queue.push(nTet);
				}
			}
		}
		for (Tet tet : tets) tet.depth = null;
		return tets;
	}
	
	public List<Tri> getListOfTris() {
		List<Tri> ts = new ArrayList<Tri>(4);
		for (Tri t : tris) {
			ts.add(t);
		}
		return ts;
	}
	
	public ArrayList<Vertex> breadthFirstVertices(int maxDepth) {
		Queue queue = new Queue();
		Queue queue2 = new Queue();
		ArrayList<Tet> tets = new ArrayList<Tet>();
		ArrayList<Vertex> vertices = new ArrayList<Vertex>();
		Tet tet, nTet;
		Vertex v;
		depth = 0;
		queue.push(this);
		
		while (!queue.isEmpty()) {
			tet = (Tet)queue.pop();
			tets.add(tet);
			for (int k = 0; k < 4; k++) {
				v = tet.corners[k];
				if (!v.flag) {
					v.flag = true;
					vertices.add(v);
				}
				nTet = tet.neighbors[k];
				if ((!nTet.isBig()) && (nTet.depth == null)) {
					nTet.depth = tet.depth + 1;
					if (nTet.depth < maxDepth) queue.push(nTet);
				}
			}
		}
		for (Tet t : tets) {
			for (int k = 0; k < 4; k++) {
				v = t.corners[k];
			}
		}
		for (Tet t : tets) t.depth = null;
		return vertices;
	}

	public int hashCode() {
		return corners[0].getId()+corners[1].getId()*2+corners[2].getId()*4+corners[3].getId()*8;
	}
}
