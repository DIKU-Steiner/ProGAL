package ProGAL.geom2d.delaunay;

import ProGAL.geom2d.Point;

/** 
 * All triangles should be counter-clock-wise (TODO: add an assert)
 * 
 * @author ras
 */
public class Triangle extends ProGAL.geom2d.Triangle{
	
	Triangle[] neighbors = new Triangle[3];
	Vertex[] corners = new Vertex[3];
	
	public Triangle(Vertex p0, Vertex p1, Vertex p2){
		super(p0,p1,p2);
		corners[0] = p0;
		if(Point.rightTurn(p0, p1, p2)){
			corners[1] = p2;
			corners[2] = p1;
			super.setCorner(p2, 1);
			super.setCorner(p1, 2);
		}else{
			corners[1] = p1;
			corners[2] = p2;
		}
	}
	
	int indexOf(Vertex v){
		for(int i=0;i<3;i++) if(corners[i]==v) return i;
		return -1;
	}
	
	public String toString(){
		return String.format("Triangle[%d,%d,%d]",corners[0].id, corners[1].id, corners[2].id);
	}
	
//	/** neighbors[i] is always opposite to points[i] */
//	protected Triangle[] neighbors;
//	Circle circumsphere;
//	private Vector ab,bc,ca;
//	
//	private final int id;
//	private static int idCounter = 0;
//	
//	Triangle(Point p1, Point p2, Point p3){
//		super(p1,p2,p3);
//		neighbors = new Triangle[3];
//		this.circumsphere = new Circle(p1,p2,p3);
//		id = idCounter++;
//	}
//
//	/** [Ericcson05 p. 205] */
//	public boolean contains(Point p){
//		if(ab==null){
//			ab = points[0].vectorTo(points[1]);
//			bc = points[1].vectorTo(points[2]);
//			ca = points[2].vectorTo(points[0]);
//		}
//		Vector ap = points[0].vectorTo(p);
//		if(cross2D(ap,ab)>0) return false;
//		
//		Vector bp = points[1].vectorTo(p);
//		if(cross2D(bp,bc)>0) return false;
//		
//		Vector cp = points[2].vectorTo(p);
//		if(cross2D(cp,ca)>0) return false;
//		
//		return true;
//	}
//	private static double cross2D(Vector u, Vector v){
//		return u.x()*v.y() - u.y()*v.x();
//	}
//	
//	public int indexOf(Point p){
//		if(super.points[0]==p) return 0;
//		if(super.points[1]==p) return 1;
//		if(super.points[2]==p) return 2;
//		return -1;
//	}
//	
//	public Triangle getNeighbor(int n){
//		return neighbors[n];
//	}
//	public Point getOpposite(int n){
//		int m = 0;
//		while(neighbors[n].neighbors[m]!=this) m++;
//		return neighbors[n].getCorner(m);
//	}
//	
//	void setNeighbors(Triangle n0, Triangle n1, Triangle n2){
//		neighbors[0] = n0;
//		neighbors[1] = n1;
//		neighbors[2] = n2;
//		for(int n=0;n<3;n++){
//			if(neighbors[n]!=null){
//				int apexIdx = (neighbors[n].indexOf(points[(n+1)%1])+1)%3;
//				neighbors[n].neighbors[apexIdx] = this;
//			}
//		}
//	}
//	
//	public String toString(){
//		return String.format("delaunay.Triangle[id=%d, %s,%s,%s, neighbors=[%d,%d,%d]]",
//				id, 
//				points[0].toString(),points[1].toString(),points[2].toString(), 
//				neighbors[0]==null?-1:neighbors[0].id, 
//				neighbors[1]==null?-1:neighbors[1].id, 
//				neighbors[2]==null?-1:neighbors[2].id 
//		);
////		return String.format("delaunay.Triangle[id=%d, %s,%s,%s, neighbors=[%d,%d,%d]]",
////				id, 
////				points[0].toString(),points[1].toString(),points[2].toString(), 
////				neighbors[0]==null?-1:neighbors[0].id, 
////				neighbors[1]==null?-1:neighbors[1].id, 
////				neighbors[2]==null?-1:neighbors[2].id 
////		);
//	}
}
