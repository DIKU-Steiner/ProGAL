package ProGAL.geom3d.kineticDelaunay;

import java.awt.Color;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Stack;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;

public class Vertex extends Point implements Comparable<Vertex>{
	private static final long serialVersionUID = 1L;
	public static enum VertexType { S, R };
	private VertexType type;
	private double polarAngle;
	private double initAngle = -1.0;
	private double cosAngle;
	private double sinAngle;
	private double polarRadius;
	private double squaredPolarRadius;
	private Tet tet;
	private Sphere sphere;
	private Integer depth = null;
	public boolean flag = false;
	
	public HashSet<Vertex> adjacentVertices = null; 
	private ArrayList<Vertex> adjacentVerticesFast = new ArrayList<Vertex>();
	private ArrayList<Tet> processedTetsFast = new ArrayList<Tet>(); 
	
	private int index;
	private static int indexCounter = 0;
	
	public Vertex(Point p) {
		super(p);
		this.index = indexCounter++;
	}

    public int getId()                               { return index; }
	public void setId(int index)                     { this.index = index; }

	public VertexType getType() { return type; }
	public void setType(VertexType type) { this.type = type; }
	public double getPolarRadius()                   { return polarRadius; }
	public void   setPolarRadius(double polarRadius) { this.polarRadius = polarRadius; }
	public double getSquaredPolarRadius()            { return squaredPolarRadius; }
	public void   setSquaredPolarRadius(double sqPR) { this.squaredPolarRadius = sqPR; }
	public double getPolarAngle()                    { return polarAngle; }
    public void   setPolarAngle(double polarAngle)   { this.polarAngle = polarAngle; }
    public double getInitAngle()                    { return initAngle; }
    public void   setInitAngle(double initAngle)   {
    	if (this.initAngle==-1.0) {
    		this.initAngle = initAngle;
    	}
    }
    public Point returnAsPoint() {
    	return new Point(this.toVector());
    }
    public double getCosAngle()                      { return cosAngle; }
    public void setCosAngle(double cosAngle)         { this.cosAngle = cosAngle; }
    public double getSinAngle()                      { return sinAngle; }
    public void setSinAngle(double sinAngle)         { this.sinAngle = sinAngle; }
    public Integer getDepth() { return depth; } 
    public void setDepth(Integer depth) { this.depth = depth; }
  /*  public HashSet<Vertex> getAdjacentVertices(double alpha2) {
    	if (adjacentVertices == null) computeAdjacentVertices(alpha2); 
    	return adjacentVertices;
    }
  */  public ArrayList<Vertex> getAdjacentVerticesFast() { return adjacentVerticesFast; }
    
    /** Returns a tetrahedron incident with this vertex */
    public Tet getTet() { return tet; }
    
    public void setTet(Tet tet) { this.tet = tet; }
    
    /** Returns a tetrahedron incident with this vertex */
    public Tet getTetrahedron(KineticAlphaComplex kDT) {
		for (Tet tet : kDT.getTetrahedra()) if (tet.hasVertex(this)) return tet;
		return null;
	}
    
 
    public void computeAdjacentVerticesFast(double alpha2) {
    	processedTetsFast.clear();
    	computeAdjacentVerticesFast(tet, alpha2);
    	for (Tet tet : processedTetsFast) tet.setFlag(false);
    	for (Vertex v : adjacentVerticesFast) v.setDepth(null);
    }
    
    public void computeAdjacentVerticesFast(Tet tet, double alpha2) {
       	Vertex b;
    	Tet nTet, nnTet;
    	Stack<Tet> stack = new Stack<Tet>();
    	stack.push(tet); tet.setFlag(true);
    	while (!stack.isEmpty()) {
    		nTet = stack.pop();
    		processedTetsFast.add(nTet);
    		for (int k = 0; k < 4; k++) {
    			b = nTet.getCorner(k);
    			if (b != this) {
    				if (b.getDepth() == null) {
    					if (distanceSquared(b) < alpha2) {
    						adjacentVerticesFast.add(b);
    						b.setDepth(1);
    					}
    				}
    				nnTet = nTet.neighbors[k];
    				if (!nnTet.getFlag()) { stack.push(nnTet); nnTet.setFlag(true); }
    			}
    		}
    	}
    }

      
	public boolean isBig() { return index < 4; }
    
	public int compareTo(Vertex arg0) {
		return index-arg0.index;
	}
	
	public void toScene(J3DScene scene, Color clr, double size) {
		scene.removeShape(sphere);
		sphere = new Sphere(this, size);
		scene.addShape(sphere, clr);
		scene.repaint();
	}
	
	public String toString(){
		return Integer.toString(index);
	}
}
