package ProGAL.geom3d.kineticDelaunay;

import java.util.Arrays;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.volumes.Tetrahedron;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;

public class Tet {
	Vertex[] corners;// = new Vertex[4];
	Tet[] neighbors = new Tet[4];
	
	public Tet(Vertex[] corners){
//		super(corners[0], corners[1],corners[2],corners[3]);
		this.corners = corners;
		this.sortCorners();
		normalizePredicates();
	}
	
	
	boolean[] oppositeInside = new boolean[4];
	boolean centerInside;

	private void normalizePredicates(){
		for(int i=0;i<4;i++){
			oppositeInside[i] = orient(i, corners[i]);
		}
		centerInside = inSphere(Point.getMidpoint(corners[0], corners[1]));
	}
	
	private boolean orient(int face, Point p){
		Matrix m = new Matrix(4,4);
		for(int r=0;r<3;r++){
			for(int c=0;c<3;c++){
				m.set(r, c, corners[(r+(r>=face?1:0))%4].getCoord(c));
			}
			m.set(r, 3, 1);
		}
		for(int c=0;c<3;c++)
			m.set(3, c, p.getCoord(c));
		m.set(3,3,1);

		double det = m.determinant();
		if (Math.abs(det)<Constants.EPSILON) return true;
		return det<0;
	}
	
	private boolean inSphere(Point p){
		Matrix m = new Matrix(5,5);
		for(int r=0;r<4;r++){
			for(int c=0;c<3;c++){
				m.set(r, c, corners[r].getCoord(c));
			}
			m.set(r, 3, corners[r].dot(corners[r]));
			m.set(r, 4, 1);
		}
		for(int c=0;c<3;c++)
			m.set(4, c, p.getCoord(c));
		m.set(4, 3, p.dot(p));
		m.set(4, 4, 1);

		double det = m.determinant();
		if (Math.abs(det)<Constants.EPSILON) return true;
		return det<0;
	}
	
	public boolean insideFace(int face, Point p){
		return orient(face, p)==oppositeInside[face];
	}
	

	public boolean insideCircumsphere(Point p){
		return inSphere(p)==centerInside;
	}
	

	public int apex(Tet c){
		for(int i=0;i<4;i++) if(neighbors[i]==c) return i;
		return -1;
	}
	
	public int indexOf_slow(Vertex v){
		for(int i=0;i<4;i++) if(corners[i]==v) return i;
		return -1;
	}

	public int indexOf(Vertex v){
		return Arrays.binarySearch(corners, v);
	}
	

	private final void swap(int i, int j){
		Vertex tmpV = corners[i];
		Tet tmpT = neighbors[i];
		corners[i] = corners[j];
		neighbors[i] = neighbors[j];
		corners[j] = tmpV;
		neighbors[j] = tmpT;
		
	}
	/** 
	 * Sorts the corners while maintaining that corners[i] is not on the shared 
	 * face between <code>this</code> and neighbors[i]. Uses 5 comparisons which is optimal. 
	 */
	protected final void sortCorners(){
		if(corners[0].compareTo(corners[1])>0) swap(0,1);//First two elements are sorted
		if(corners[2].compareTo(corners[3])>0) swap(2,3);//Last two elements are sorted
		if(corners[0].compareTo(corners[2])>0) swap(0,2);//Smallest element is at 0
		if(corners[1].compareTo(corners[3])>0) swap(1,3);//Largest element is at 3
		if(corners[1].compareTo(corners[2])>0) swap(1,2);//Check middle elements
	}
	
	
	public String toString(){
//		StringBuilder sb = new StringBuilder();
//		sb.append("[");
//		for(Vertex v: corners){
//			sb.append(v.index);
//			sb.append(',');
//		}
//		sb.append("]");
//		return sb.toString();
		return Arrays.toString(corners);
	}
}
