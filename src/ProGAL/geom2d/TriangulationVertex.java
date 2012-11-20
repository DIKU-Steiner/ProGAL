package ProGAL.geom2d;

import java.util.ArrayList;
import java.util.List;

public class TriangulationVertex extends Point {
	public static enum VertexType { S, R };

	
	protected int id;
	protected TriangulationFace face;
	protected Circle orbit = null;
	private double polarAngle;
	private double cosAngle;
	private double sinAngle;
	private double polarRadius;
	private double squaredPolarRadius;
	private VertexType type;
	
	public TriangulationVertex(double x, double y) {
		super(x, y);
	}
	
	public TriangulationFace getFace() { return face; }
	public void setFace(TriangulationFace face) { this.face = face; }

	public TriangulationFace getNextFace(TriangulationFace f) {
		int indx = f.getIndex(this);
		return f.getNeighbor((indx+1)%3);
	}
	
	public TriangulationFace getLastFace() {
		TriangulationFace currFace = face;
		TriangulationFace nextFace = getNextFace(currFace);
		while ((nextFace != null) && (nextFace != face)) { 
			currFace = nextFace;
			nextFace = getNextFace(currFace);
		}
		if (nextFace == null) return currFace; else return face;
	}
	
	public List<TriangulationFace> getFaces() {
		List<TriangulationFace> faces = new ArrayList<TriangulationFace>();
		TriangulationFace firstFace = getFace();
		faces.add(firstFace);
		TriangulationFace face = getNextFace(firstFace);
		while ((face != null) && (face != firstFace)) {
			faces.add(face);
			face = getNextFace(face);
		}
		return faces;
	}
	
	public Circle getOrbit(Point p) { 
		if (orbit == null) orbit = new Circle(p, distance(p)); 
		return orbit;
	}

    public int getId()                               { return id; }
	public void setId(int id)                        { this.id = id; }
	public VertexType getType()                      { return type; }
	public void       setType(VertexType type)       { this.type = type; }
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

    
}
