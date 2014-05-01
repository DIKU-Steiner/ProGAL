package ProGAL.geom3d.complex;

import java.awt.Color;

import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Tetrahedron;

public class CTriangle extends Triangle {
	private CTetrahedron[] adjacentTetrahedra = new CTetrahedron[2];
	private CEdge[] edges = new CEdge[3];


	public CTriangle(CVertex p0, CVertex p1, CVertex p2, CTetrahedron t1, CTetrahedron t2) {
		super(p0,p1,p2);
		orderPoints(p0,p1,p2);
		
		this.adjacentTetrahedra[0] = t1;
		this.adjacentTetrahedra[1] = t2;
	}

	public void setNeighbour(int index, CTetrahedron t){ 
		this.adjacentTetrahedra[index]=t;
	}

	public CTetrahedron getAdjacentTetrahedron(int index){ 
		return this.adjacentTetrahedra[index];
	}
	
	/*public CTriangle[] getNeighborTriangles() {
		CTriangle[] triangles = new CTriangle[3];
		return triangles;
	}*/

	public CEdge getEdge(int i){
		return edges[i];
	}

	public void setEdge(int i, CEdge e){
		edges[i] = e;
	}
	
	/** TODO: Copy to Triangle */
	public CVertex oppositeVertex(CEdge e){
		if(!e.containsPoint(p1)) 	return (CVertex)p1;
		if(!e.containsPoint(p2)) 	return (CVertex)p2;
		else						return (CVertex)p3;
	}


	private void orderPoints(CVertex a, CVertex b, CVertex c){
		CVertex p[] = new CVertex[3];
		if(a.dominates(b) && a.dominates(c)){
			p[0]=a;
			if(b.dominates(c)) { p[1]=b; p[2]=c; }
			else        { p[1]=c; p[2]=b; }
		}
		else if(b.dominates(a) && b.dominates(c)){
			p[0]=b;
			if(a.dominates(c)) { p[1]=a; p[2]=c; }
			else        { p[1]=c; p[2]=a; }
		}
		else{
			p[0]=c;
			if(a.dominates(b)) { p[1]=a; p[2]=b; }
			else        { p[1]=b; p[2]=a; }
		}
		super.p1 = p[0];
		super.p2 = p[1];
		super.p3 = p[2];
	}

	
	public boolean equals(Object o){
		if(!(o instanceof CTriangle)) return false;
		return ((CTriangle)o).p1==p1&&((CTriangle)o).p2==p2&&((CTriangle)o).p3==p3;
	}
	public int hashCode(){
		return p1.hashCode()^p2.hashCode()^p3.hashCode();
	}

	public boolean containsPoint(CVertex point) {
		return (p1==point || p2==point || p3==point);
	}
	public boolean containsBigPoint(){
		return ((CVertex)p1).isBigpoint() || ((CVertex)p2).isBigpoint() || ((CVertex)p3).isBigpoint(); 
	}
	
	public void toScene(J3DScene scene, double rad, Color clr) { }

	public Tetrahedron toScene(J3DScene scene, Color clr) {
		Tetrahedron ret = new Tetrahedron(getCorner(0), getCorner(1), getCorner(2), getCenter());
		scene.addShape(ret, clr);
		return ret;
	}
}
