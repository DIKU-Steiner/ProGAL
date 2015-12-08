package ProGAL.geom2d.delaunay;

import ProGAL.geom2d.Point;

public class Vertex extends ProGAL.geom2d.Point{
	private static final long serialVersionUID = 1L;
	int id;
	static int idCounter = 0;
	
	Triangle first;
	Triangle last;
	
	Vertex(Point p){
		super(p.x(), p.y());
		id = idCounter++;
	}
	
	public boolean onBoundary(){
		return first!=last || first.neighbors[ (first.indexOf(this)+2)%3 ] == null;
	}
	
	public String toString(){
		return String.format("V[%d]",id);
	}
}