package ProGAL.geom3d.kineticDelaunay;

import java.awt.Color;
import java.util.Arrays;

import ProGAL.geom3d.Shape;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;

public class Edge {
	private Vertex[] corners = new Vertex[2];
	Integer count = null;
	boolean dAlive = true;
	boolean cAlive = true;
//	Tet tet = null;
	private Shape LSSs;
	boolean alph = false;
	
	public Edge(Vertex[] corners){
		this.corners = corners;
		this.sortCorners();
	}
	
	public Edge(Vertex v0, Vertex v1) {
		corners[0] = v0;
		corners[1] = v1;
		this.sortCorners();
	}
	
	private void sortCorners() {
		if (corners[0].getId()>corners[1].getId()) {
			Vertex tmp = corners[0];
			corners[0] = corners[1];
			corners[1] = tmp;
		}
	}
	
	public Vertex getCorner(int i) { 
		return corners[i];
	}
	
//	public void setTet(Tet t) {
//		tet = t;
//	}
	
//	public Tet getTet() { return tet; }
	
	public void setAlph(boolean b) {
		alph = b;
	}
	public boolean getAlph() {
		return alph;
	}
	
	public int hashCode() {
		return corners[0].getId()+corners[1].getId()*2;
	}
	
	/** Returns the radius  of the smallest sphere circumscribing this edge */
	public double getCircumRadius() { 
		return corners[0].distance(corners[1])/2;
	}
	
	public boolean hasVertex(Vertex v) {
		for (int i = 0; i < 2; i++) 
			if (corners[i] == v) return true;
		return false;
	}
	
	public boolean isAlive() { return dAlive; }
	public void setAlive(boolean dAlive) {this.dAlive = dAlive; }
	
	/** Returns TRUE if one of the corners is a big point - there are 4 big points with id's 0, 1, 2 , 3*/
	public boolean isBig() {
		return ((corners[0].getId() < 4) || (corners[1].getId() < 4));
	}
	
	/** returns TRUE if this tetrahedron is in the alpha complex for the given alpha*/
	public boolean isAlpha(double alpha) { return (getCircumRadius() < (alpha+ProGAL.math.Constants.EPSILON)); }
	
	public int getCount() { 
		if (count != null) return count; 
		else {
			count = 0;
			for (int i = 0; i < 2; i++) { if (corners[i].getType() == Vertex.VertexType.R) count = count + (int)Math.pow(2,i);  }
			return count;
		}
	}
	
	public double getLength() {
		return corners[0].distance(corners[1]);
	}
	
	public double getLengthSquared() {
		return corners[0].distanceSquared(corners[1]);
	}
	
	/* Returns index of vertex v */
	public int indexOf_slow(Vertex v){
		return Arrays.binarySearch(corners, v);
	}
	
	public String toString(){
		return Arrays.toString(corners);
	}
	
	public void toConsole() {
		System.out.println(toString());
	}
	
	public void fromSceneEdge(J3DScene scene) {
		scene.removeShape(LSSs);
		scene.repaint();
	}
	
	public boolean equals(Object o){
		if(!(o instanceof Edge)) return false;
		return ((Edge)o).getCorner(0).equals(corners[0])&&((Edge)o).getCorner(1).equals(corners[1]);
	}
	
	/* draws the edge */
	public void toSceneEdge(J3DScene scene, Color clr, double width) {
		LSSs = new LSS(corners[0], corners[1], width);
		scene.addShape(LSSs, clr, 3);
		scene.repaint();
	}

}
