package ProGAL.geom3d.complex;

import java.util.ArrayList;
import java.util.List;

import ProGAL.geom3d.LineSegment;
import ProGAL.geom3d.Point;

public class CEdge extends LineSegment{

	private final List<CTriangle> adjacentTriangles = new ArrayList<CTriangle>();

	public CEdge(CVertex p0, CVertex p1){
		super(p0,p1);
	}

	public void addTriangle(CTriangle tri){
		adjacentTriangles.add(tri);
	}
	
	public List<CTriangle> getAdjacentTriangles() {
		return adjacentTriangles;
	}
	
	public CVertex getPoint(int i){
		return (CVertex)super.getPoint(i);
	}

	public boolean containsPoint(Point p){
		return a==p||b==p;
	}
	
	//Equals and hashCode used in DelaunayComplex to uniquely specify an edge
	public boolean equals(Object o){
		if(!(o instanceof CEdge)) return false;
		return ( ((CEdge)o).a==a&&((CEdge)o).b==b )||( ((CEdge)o).a==b&&((CEdge)o).b==a );
	}
	public int hashCode(){
		return a.hashCode()^b.hashCode();
	}

	public CVertex opposite(CVertex v) {
		if(v==a) return (CVertex)b;
		if(v==b) return (CVertex)a;
		throw new Error("Vertex is not an end-point of this edge");
	}
}

