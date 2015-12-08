package ProGAL.geom3d.complex.delaunayComplex;

import ProGAL.geom3d.complex.CTetrahedron;

class Flip {
	private int id;
	private CTetrahedron t;
	
	public int getId() {
		return id;
	}
	
	public void setId(int id) {
		this.id = id;
	}
	
	public CTetrahedron getT() {
		return t;
	}
	
	public void setT(CTetrahedron t) {
		this.t = t;
	}
	
	public Flip(int id, CTetrahedron t) {
		super();
		this.id = id;
		this.t = t;
	}
	
	
	
}
