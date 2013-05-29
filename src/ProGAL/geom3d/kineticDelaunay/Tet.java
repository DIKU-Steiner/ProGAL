package ProGAL.geom3d.kineticDelaunay;

import java.awt.Color;
import java.util.Arrays;

import ProGAL.geom2d.TriangulationVertex;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;

public class Tet {
	Vertex[] corners = new Vertex[4];
	Tet[] neighbors = new Tet[4];
	Sphere circumSphere = null;
	Integer count = null;
	boolean dAlive = true;
	boolean cAlive = true;
	Shape[] LSSs = new Shape[6];
	Shape[] faces = new Shape[4];
	
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
		for (int i = 0; i < 4; i++) {
			Vertex p = new Vertex(tetra.getCorner(i));
			corners[i] = p;
		}
		this.sortCorners();		
		normalizePredicates();
	}

	public Tet clone() {
		Tet nt = new Tet(corners);
		for (int i = 0; i < 4; i++) nt.neighbors[i] = neighbors[i];
		return nt;
	}
	
	public Tetrahedron makeTetrahedron() {
		return new Tetrahedron(corners[0], corners[1], corners[2], corners[3]);
	}
	
	boolean[] oppositeInside = new boolean[4];
	boolean centerInside;

	public Sphere getCircumSphere() { return circumSphere; }
	
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

	
	public boolean hasVertex(Vertex v) {
		for (int i = 0; i < 4; i++) 
			if (corners[i] == v) return true;
		return false;
	}

	public boolean hasNeighbor(Tet t) {
		for (int i = 0; i < 4; i++) 
			if (neighbors[i] == t) return true;
		return false;
	}

	public boolean isAlive() { return dAlive; }
	public void setAlive(boolean dAlive) {this.dAlive = dAlive; }
	
	public boolean isFlat() { return Point.coplanar(corners[0], corners[1], corners[2], corners[3]); }
	
	/* returns true if this tetrahedron and its neighbor tet form a convex set. */
	public boolean isConvex(Tet tet) {
		Triangle tr = commonFace(tet);
		Point p = corners[apex(tet)];
		Point q = tet.corners[tet.apex(this)];
		return tr.getIntersection(p, new Vector(p, q)) != null;
	}
	
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
		for(int i=0;i<4;i++) if(neighbors[i]==c) return i;
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
	
	public int indexOf_slow(Vertex v){
		for(int i=0;i<4;i++) if(corners[i]==v) return i;
		return -1;
	}

	public int indexOf(Vertex v){
		return Arrays.binarySearch(corners, v);
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
//		StringBuilder sb = new StringBuilder();
//		sb.append("[");
//		for(Vertex v: corners){
//			sb.append(v.index);
//			sb.append(',');
//		}
//		sb.append("]");
//		return sb.toString();
		return Arrays.toString(corners);
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
		for (int i = 0; i < 6; i++) scene.removeShape(LSSs[i]);
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
		for (int i = 0; i < 4; i++) scene.addShape(faces[i], clr);
		scene.repaint();
	}
	
	public void toSceneFace(J3DScene scene, int i, Color clr) {
		if (faces[i] == null) faces[i] = new Triangle(corners[(i+1)%4], corners[(i+2)%4], corners[(i+3)%4]);
		scene.addShape(faces[i], clr);
	}
	
	public void fromSceneFaces(J3DScene scene) {
		for (int i = 0; i < 4; i++) scene.removeShape(faces[i]);
		scene.repaint();
	}
	
	public void fromSceneFaces(J3DScene scene, int i) {
		scene.removeShape(faces[i]);
		scene.repaint();
	}

}
