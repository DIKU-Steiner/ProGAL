package ProGAL.geom2d;

import java.util.ArrayList;
import java.util.List;

import ProGAL.dataStructures.DLCyclicList;
import ProGAL.geom3d.volumes.Sphere;

public class TriangulationVertex extends Point {
	private static final long serialVersionUID = 1L;
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
	protected boolean burried = false; // used in connection with deletion
	protected ProGAL.geom3d.Point liftedPoint;
	protected ProGAL.geom3d.Point groundPoint;
	protected Sphere liftedSphere;
	protected Sphere groundSphere;
	protected boolean bigPoint;

	public TriangulationVertex(double x, double y) {
		super(x, y);
		bigPoint = false;
	}
	
	public TriangulationVertex(Point p) {
		super(p.x(), p.y());
		bigPoint = false;
	}
	
	public TriangulationFace getFace() { return face; }
	
	public void setFace(TriangulationFace face) { this.face = face; }

	public boolean isBigPoint() { return bigPoint; }
	public void setBigPoint(boolean bigPoint) { this.bigPoint = bigPoint; }
	
	public TriangulationFace getNextFace(TriangulationFace f) {
		int indx = f.getIndex(this);
		return f.getNeighbor((indx+1)%3);
	}
	
	public TriangulationFace getPrevFace(TriangulationFace f) {
		int indx = f.getIndex(this);
		return f.getNeighbor((indx+2)%3);
		
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
	
	/* returns in counterclockwise order the faces incident with this vertex */
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
	
	/* returns counterclockwise list of vertices adjacent to this vertex */
	public DLCyclicList<TriangulationVertex> getNeighboringVertices() {
		DLCyclicList<TriangulationVertex> vertices = new DLCyclicList<TriangulationVertex>();
		TriangulationFace currentFace = getFace();
		TriangulationVertex u = currentFace.getCorner((face.getIndex(this)+1)%3);
		TriangulationVertex v = currentFace.getCorner((face.getIndex(this)+2)%3);
		vertices.pushBefore(v);
		while (v != u) {
			currentFace = getNextFace(currentFace);
			v = currentFace.getThirdVertex(this, v);
			vertices.pushBefore(v);
		}
		return vertices;
	}
	
	public Circle getOrbit() { 
		if (orbit == null) orbit = new Circle(Point.origo, distance()); 
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
