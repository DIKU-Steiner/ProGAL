package ProGAL.geom3d.kineticDelaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ProGAL.geom3d.Circle;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;

/**
 * 
 * @author Daisy
 */
public class Tri {
	private Vertex[] corners = new Vertex[3];
//	private List<Edge> edges = new ArrayList<Edge>(3); // Daisy
	Integer count = null;
	boolean dAlive = true;
	boolean cAlive = true;
//	private List<Tet> tets = new ArrayList<Tet>(2); // Daisy
	private Shape[] LSSs = new Shape[3];
	private Triangle triangle;
	int alph = 0; // Daisy: 0 and 2 = not in alpha complex
	
	public Tri(Vertex[] corners){
		this.corners = corners;
		this.sortCorners();
	}
	
	public Tri(Vertex v0, Vertex v1, Vertex v2) {
		corners[0] = v0;
		corners[1] = v1;
		corners[2] = v2;
		this.sortCorners();
	}
	
	private void sortCorners() {
		int id0 = corners[0].getId();
		int id1 = corners[1].getId();
		int id2 = corners[2].getId();
		Vertex tmp;
		if (id0>id1) {
			if (id2>id1) {
				tmp = corners[0];
				corners[0] = corners[1];
				if (id0>id2) {
					corners[1] = corners[2];
					corners[2] = tmp;
				} else corners[1] = tmp;
			} else {
				tmp = corners[0];
				corners[0] = corners[2];
				corners[2] = tmp;
			}
		}
		if (id1>id2) {
			if (id0>id2){
				tmp = corners[0];
				corners[0] = corners[2];
				corners[2] = corners[1];
				corners[1] = tmp;
			} else {
				tmp = corners[1];
				corners[1] = corners[2];
				corners[2] = tmp;
			}
		}
	}
	
	public Vertex getCorner(int i) { 
		return corners[i];
	}
	
//	public void setTet(Tet t) {
//		tets.add(t);
//	}
	
//	public List<Tet> getTets() { return tets; }
	
//	public boolean hasBigNeighbour() {
//		return (tets.get(0).isBig() || tets.get(1).isBig());
//	}
	
	// Daisy
		public void setAlph(int i) {
			alph = i;
		}
		public int getAlph() {
			return alph;
		}
	
	/** Returns the radius  of the smallest sphere circumscribing this triangle */
	public double getCircumRadius() {
		return (new Circle(corners[0], corners[1], corners[2])).getRadius();
	}
	
	public boolean hasVertex(Vertex v) {
		for (int i = 0; i < 3; i++) 
			if (corners[i] == v) return true;
		return false;
	}
	
	public boolean isAlive() { return dAlive; }
	public void setAlive(boolean dAlive) {this.dAlive = dAlive; }
	
	/** Returns TRUE if one of the corners is a big point - there are 4 big points with id's 0, 1, 2 , 3*/
	public boolean isBig() {
		return ((corners[0].getId() < 4) || (corners[1].getId() < 4));
	}
	
	/** returns TRUE if this triangle is in the alpha complex for the given alpha*/
	public boolean isAlpha(double alpha) { return (getCircumRadius() < alpha); }
	
	public int getCount() { 
		if (count != null) return count; 
		else {
			count = 0;
			for (int i = 0; i < 3; i++) { if (corners[i].getType() == Vertex.VertexType.R) count = count + (int)Math.pow(2,i);  }
			return count;
		}
	}
	
	/* Returns index of vertex v */
	public int indexOf_slow(Vertex v){
		return Arrays.binarySearch(corners, v);
	}
	
	public String toString(){
		return Arrays.toString(corners);
	}
	
	//Daisy
	public boolean equals(Object o){
		if(!(o instanceof Tri)) return false;
		return ((Tri)o).getCorner(0).equals(corners[0])&&((Tri)o).getCorner(1).equals(corners[1])&&((Tri)o).getCorner(2).equals(corners[2]);
	}
	public int hashCode() {
		return corners[0].getId()+corners[1].getId()*2+corners[2].getId()*4;
	}
	
	public void toConsole() {
		System.out.println(toString());
	}
	
	public void fromSceneEdges(J3DScene scene) {
		for (int i = 0; i < 3; i++) {
			scene.removeShape(LSSs[i]);
		}
		scene.repaint();
	}
	
	public void fromSceneEdge(J3DScene scene) {
		for (int i = 0; i < 3; i++) {
			scene.removeShape(LSSs[i]);
		}
		scene.repaint();
	}
	
	public void toSceneEdges(J3DScene scene, Color clr, double width) {
		int k = 0;
		for (int i = 0; i < 2; i++)
			for (int j = i+1; j < 3; j++) {
				if (LSSs[k] == null) LSSs[k] = new LSS(corners[i], corners[j], width);
				scene.addShape(LSSs[k++], clr, 3);
			}
		scene.repaint();
	}
	
	public void toScene(J3DScene scene, Color clr) {
		triangle = new Triangle(corners[0], corners[1], corners[2]);
		scene.addShape(triangle, clr);
		scene.repaint();
	}
	
	public void fromScene(J3DScene scene) {
		scene.removeShape(triangle);
		scene.repaint();
	}
	
	/* draws the edges */
	public void toSceneEdges(J3DScene scene, Color clr, double width, double bigWidth) {
		double edgeWidth = width;
		for (int i = 0; i < 3; i++)  if (corners[i].getId() < 4) edgeWidth = bigWidth;
		int k = 0;
		for (int i = 0; i < 2; i++)
			for (int j = i+1; j < 3; j++) {
				if (LSSs[k] == null) LSSs[k] = new LSS(corners[i], corners[j], edgeWidth);
				scene.addShape(LSSs[k++], clr, 3);
			}
		scene.repaint();
	}
}
