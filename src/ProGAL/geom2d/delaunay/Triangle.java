package ProGAL.geom2d.delaunay;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.Vector;

/** 
 * All triangles should be counter-clock-wise (TODO: add an assert)
 * 
 * @author ras
 */
public class Triangle extends ProGAL.geom2d.Triangle{
	/** neighbors[i] is always opposite to points[i] */
	protected Triangle[] neighbors;
	Circle circumsphere;
	private Vector ab,bc,ca;
	
	private final int id;
	private static int idCounter = 0;
	
	Triangle(Point p1, Point p2, Point p3){
		super(p1,p2,p3);
		neighbors = new Triangle[3];
		this.circumsphere = new Circle(p1,p2,p3);
		id = idCounter++;
	}

	/** [Ericcson05 p. 205] */
	public boolean contains(Point p){
		if(ab==null){
			ab = points[0].vectorTo(points[1]);
			bc = points[1].vectorTo(points[2]);
			ca = points[2].vectorTo(points[0]);
		}
		Vector ap = points[0].vectorTo(p);
		if(cross2D(ap,ab)>0) return false;
		
		Vector bp = points[1].vectorTo(p);
		if(cross2D(bp,bc)>0) return false;
		
		Vector cp = points[2].vectorTo(p);
		if(cross2D(cp,ca)>0) return false;
		
		return true;
	}
	private static double cross2D(Vector u, Vector v){
		return u.x()*v.y() - u.y()*v.x();
	}
	
	public Triangle getNeighbor(int n){
		return neighbors[n];
	}
	public Point getOpposite(int n){
		int m = 0;
		while(neighbors[n].neighbors[m]!=this) m++;
		return neighbors[n].getCorner(m);
	}
	
	void setNeighbors(Triangle n0, Triangle n1, Triangle n2){
		neighbors[0] = n0;
		neighbors[1] = n1;
		neighbors[2] = n2;
		for(int n=0;n<3;n++){
			if(neighbors[n]!=null){
				//Find the corner of neighbors[n] opposite any two points in this triangle
				int oppositeIdx = -1;
				for(int i=0;i<3;i++) {
					Point opposite= neighbors[n].points[i];
					if(opposite!=points[0] && opposite!=points[1] && opposite!=points[2]) {
						oppositeIdx = i;
//						opposites[i] = opposite;
						break;
					}
				}
				//Neighbor-indices are opposite to corner-indices
				neighbors[n].neighbors[oppositeIdx] = this;
			}
		}
	}
	
	public String toString(){ return toString(2); }
	public String toString(int dec){
		return String.format("delaunay.Triangle[id=%d, %s,%s,%s, neighbors=[%d,%d,%d]]",
				id, 
				points[0].toString(dec),points[1].toString(dec),points[2].toString(dec), 
				neighbors[0]==null?-1:neighbors[0].id, 
				neighbors[1]==null?-1:neighbors[1].id, 
				neighbors[2]==null?-1:neighbors[2].id 
		);
	}
}
