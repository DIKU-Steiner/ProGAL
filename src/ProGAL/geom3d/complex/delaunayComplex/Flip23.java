package ProGAL.geom3d.complex.delaunayComplex;

import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CVertex.DegenerateCase;

class Flip23 {
	private int a1, b1, cid;
	private CTetrahedron t3=null;
	
	private Flips flips;
	

	public Flip23(Flips flips) {
		super();		
		this.flips=flips;		
		t3=null;
		
	}

	public CTetrahedron getT3() {
		return t3;
	}


	public void setT3(CTetrahedron t3) {
		this.t3 = t3;
	}


	public int getCid() {
		return cid;
	}


	public void setCid(int cid) {
		this.cid = cid;
	}


	public int getA1() {
		return a1;
	}


	public void setA1(int a1) {
		this.a1 = a1;
	}


	public int getB1() {
		return b1;
	}


	public void setB1(int b1) {
		this.b1 = b1;
	}


	public CTetrahedron flip23(CTetrahedron t, int pid, int did){
		CTetrahedron next=null;

		CVertex p = t.getPoint(pid);
		CTetrahedron t2 = t.getNeighbour(pid);
		CVertex d = t2.getPoint(did);

		int aid = (pid+1)%4;
		int bid = (pid+2)%4;
		int cid = (pid+3)%4;

		int a2id = findpoint(t2,t.getPoint(aid));
		int b2id = findpoint(t2,t.getPoint(bid));
		int c2id = findpoint(t2,t.getPoint(cid));
		
		CVertex a = t.getPoint(aid);
		CVertex b = t.getPoint(bid);
		CVertex c = t.getPoint(cid);
		
		CTetrahedron pabd = new CTetrahedron(p,a,b,d);
		CTetrahedron pacd = new CTetrahedron(p,a,c,d);
		CTetrahedron pbcd = new CTetrahedron(p,b,c,d);
		
		if(p.isDegenerate()){
			if(p.getDegCase()==DegenerateCase.ONEDGE){
				//TODO: There is some code missing from the C++ here
			}			
		}
		
		pabd.setNeighbour(0, t2.getNeighbour(c2id));
		pabd.setNeighbour(1, pbcd);
		pabd.setNeighbour(2, pacd);
		pabd.setNeighbour(3, t.getNeighbour(cid));

		pacd.setNeighbour(0, t2.getNeighbour(b2id));
		pacd.setNeighbour(1, pbcd);
		pacd.setNeighbour(2, pabd);
		pacd.setNeighbour(3, t.getNeighbour(bid));
		
		pbcd.setNeighbour(0, t2.getNeighbour(a2id));
		pbcd.setNeighbour(1, pacd);
		pbcd.setNeighbour(2, pabd);
		pbcd.setNeighbour(3, t.getNeighbour(aid));
		
		if(t.getNeighbour(aid)!=null) 	(t.getNeighbour(aid)).updateNeighbour(t, pbcd);
		if(t2.getNeighbour(a2id)!=null) (t2.getNeighbour(a2id)).updateNeighbour(t2, pbcd);
		if(t.getNeighbour(bid)!=null) 	(t.getNeighbour(bid)).updateNeighbour(t, pacd);
		if(t2.getNeighbour(b2id)!=null) (t2.getNeighbour(b2id)).updateNeighbour(t2, pacd);
		if(t.getNeighbour(cid)!=null) 	(t.getNeighbour(cid)).updateNeighbour(t, pabd);
		if(t2.getNeighbour(c2id)!=null) (t2.getNeighbour(c2id)).updateNeighbour(t2, pabd);
		
		flips.addFlip(0,pabd);
		flips.addFlip(0, pacd);
		flips.addFlip(0, pbcd);
		
		flips.addTetrahedron(pabd);
		flips.addTetrahedron(pacd);
		flips.addTetrahedron(pbcd);
		
		t.setModified(true);
		t2.setModified(true);
		
		if(next!=null) return next;
		
		return pabd;
	}


	public int findpoint(CTetrahedron t, CVertex p){
		for(int i = 0; i<4; i++){
			if(t.getPoint(i)==p) {
				return i;
			}
		}
		System.out.println("Problemer med findpoint\n");
		//never happends:
		return -1;
	}

}