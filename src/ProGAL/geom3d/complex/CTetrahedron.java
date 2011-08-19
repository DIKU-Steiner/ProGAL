package ProGAL.geom3d.complex;

import ProGAL.geom3d.volumes.Tetrahedron;


public class CTetrahedron extends Tetrahedron{
	private CTetrahedron[] neighbours = new CTetrahedron[4];
	private CTriangle[] triangles = new CTriangle[4];
	private boolean modified = false;
	private boolean flat = false;


	public CTetrahedron(CVertex p0, CVertex p1, CVertex p2, CVertex p3) {
		super(p0,p1,p2,p3);
	}
	protected CTetrahedron(){
		this(null,null,null,null);
	}

	
	
	public void setFlat(boolean flat) {				this.flat = flat;			}
	public void setModified(boolean modified) {		this.modified = modified;	}

	public void setPoint(CVertex p, int i){					super.corners[i] = p;	}
	public void setNeighbour(int index, CTetrahedron t){	neighbours[index] = t;	}
	public void setTriangle(int index, CTriangle t){		triangles[index] = t;	}
	
	public CVertex getPoint(int i){						return (CVertex)corners[i];	}
	public CTetrahedron getNeighbour(int index) {		return neighbours[index];	}
	public CTriangle getTriangle(int index){			return triangles[index];	}

	public boolean isModified() {			return modified;			}
	public boolean isFlat() {				return flat;				}
	

	public boolean containsBigPoint() {
		if(getPoint(0).isBigpoint() || getPoint(1).isBigpoint() || getPoint(2).isBigpoint() || getPoint(3).isBigpoint()) return true;
		return false;
	}

	
	
	public void updateNeighbour(CTetrahedron lookfor, CTetrahedron replacement){
		for(int i=0; i<4;i++){
			if(neighbours[i]==lookfor){
				neighbours[i]=replacement;
				break;
			}
		}
	}
	
	//find id of point
	public int findpoint(CVertex p){
		for(int i = 0; i<4; i++){
			if(getPoint(i)==p) {
				return i;
			}
		}
		System.out.println("Problemer med findpoint\n");
		//never happends:
		return -1;
	}

	public boolean containsTriangle(CTriangle t){
		for(int tp=0;tp<3;tp++){
			boolean found = false;
			for(int p=0;p<4;p++) if(this.getPoint(p)==t.getPoint(tp)) { found=true; break; }
			if(!found) return false;
		}
		return true;
	}
	
	/** TODO: Copy to Tetrahedron */
	public CVertex oppositeVertex(CTriangle base){
		for(int p=0;p<4;p++){
			if(!base.containsPoint(getPoint(p))) return getPoint(p); 
		}
		throw new RuntimeException("The triangle is not part of this tetrahedron");
	}

	public CTriangle oppositeTriangle(CVertex v) {
		for(CTriangle t: triangles){
			if(t!=null && !t.containsPoint(v)) return t;
		}
		throw new RuntimeException("The vertex is not part of this tetrahedron");
	}
	
	//given a point index this method finds the index of the apex - meaning the opposite point id that is in the tetrahedron opposite the given point id 
	//input: point index 
	//output: point index of the point opposite
	public int apexid(int index){
		//Point ap0,ap1,ap2,ap3;
		CTetrahedron apex_tet= getNeighbour(index);
		if(apex_tet!= null){
			for(int i=0;i<4;i++){
				if(apex_tet.getNeighbour(i)== this){
					return i;
				}
			}
		}
		//never happens:
		return -1;

	}
	

}
