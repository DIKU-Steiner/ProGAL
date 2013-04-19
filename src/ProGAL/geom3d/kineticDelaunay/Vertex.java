package ProGAL.geom3d.kineticDelaunay;

import java.util.ArrayList;
import java.util.List;

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
		List<Tet> tets = new ArrayList<Tet>();
		TriangulationFace firstTet = getTet();
		tets.add(firstTet);
		Tet face = getNextTet(firstTet);
    }
    
	public int compareTo(Vertex arg0) {
		return index-arg0.index;
	}
	
	public String toString(){
		return Integer.toString(index);
	}
}
