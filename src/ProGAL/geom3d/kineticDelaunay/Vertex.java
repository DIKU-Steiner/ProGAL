package ProGAL.geom3d.kineticDelaunay;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import ProGAL.geom2d.TriangulationFace;
import ProGAL.geom3d.Point;

public class Vertex extends Point implements Comparable<Vertex>{
	private static final long serialVersionUID = 1L;
	public static enum VertexType { S, R };
	private VertexType type;
	private double polarAngle;
	private double cosAngle;
	private double sinAngle;
	private double polarRadius;
	private double squaredPolarRadius;


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
    public double getCosAngle()                      { return cosAngle; }
    public void setCosAngle(double cosAngle)         { this.cosAngle = cosAngle; }
    public double getSinAngle()                      { return sinAngle; }
    public void setSinAngle(double sinAngle)         { this.sinAngle = sinAngle; }

    
    
    public List<Tet> getTets() {
    	Tet tet = getFirstTetrahedron();
    	return getTets(tet);
    }
    
    /* returns list of tetrahedra having this vertex as a corner given one tetrahedron with this vertex */
    public List<Tet> getTets(Tet tet) {
    	List<Tet> incidentTetrahedra = new ArrayList<Tet>();
    	Stack<Tet> stack = new Stack<Tet>();
    	incidentTetrahedra.add(tet);
    	stack.push(tet);
    	Tet cTet, nTet;
    	while (!stack.isEmpty()) {
    		cTet = stack.pop();
    		for (int i = 0; i < 4; i++) {
    			if (cTet.corners[i] != this) {
    				nTet = cTet.neighbors[i];
    				if (!incidentTetrahedra.contains(nTet)) {
    					incidentTetrahedra.add(nTet);
    					stack.push(nTet);
    				}
    			}
    		}
    	}
    	return incidentTetrahedra;
    }
    
	public int compareTo(Vertex arg0) {
		return index-arg0.index;
	}
	
	public String toString(){
		return Integer.toString(index);
	}
}
