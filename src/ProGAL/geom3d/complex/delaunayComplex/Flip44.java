package ProGAL.geom3d.complex.delaunayComplex;

import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.CTetrahedron;

class Flip44 {
	
	private Flips flips;
	private Flip23 f23;

	public Flip44(Flip23 f23, Flips flips) {
		super();
		this.flips = flips;
		this.f23 = f23;
	}

	public void flip44(CTetrahedron t1, CTetrahedron t2, CTetrahedron t3, CTetrahedron t4, CVertex p, CVertex d){
	
		int a1 = f23.getA1();
		int b1 = f23.getB1();
		int cid = f23.getCid();
		
		CVertex a = t1.getPoint(a1);
		CVertex b = t1.getPoint(b1);
		CVertex c = t1.getPoint(cid);
		
		int eid = findlastid(d,a,b,t3);
		CVertex e= t3.getPoint(eid);
		
		CTetrahedron dacp = new CTetrahedron(d, a, c, p);
		CTetrahedron dbcp = new CTetrahedron(d, b, c, p);
		CTetrahedron dbep = new CTetrahedron(d, b, e, p);
		CTetrahedron daep = new CTetrahedron(d, a, e, p);
		
		dacp.setNeighbour(0, t1.getNeighbour(b1));
		dacp.setNeighbour(1, dbcp);
		dacp.setNeighbour(2, daep);
		dacp.setNeighbour(3, t2.getNeighbour(t2.findpoint(b)));
		
		dbcp.setNeighbour(0, t1.getNeighbour(a1));
		dbcp.setNeighbour(1, dacp);
		dbcp.setNeighbour(2, dbep);
		dbcp.setNeighbour(3, t2.getNeighbour(t2.findpoint(a)));
		
		dbep.setNeighbour(0, t4.getNeighbour(t4.findpoint(a)));
		dbep.setNeighbour(1, daep);
		dbep.setNeighbour(2, dbcp);
		dbep.setNeighbour(3, t3.getNeighbour(t3.findpoint(a)));
		
		daep.setNeighbour(0, t4.getNeighbour(t4.findpoint(b)));
		daep.setNeighbour(1, dbep);
		daep.setNeighbour(2, dacp);
		daep.setNeighbour(3, t3.getNeighbour(t3.findpoint(b)));
		
		if(t1.getNeighbour(a1)!=null) 				t1.getNeighbour(a1).updateNeighbour(t1, dbcp);
		if(t1.getNeighbour(b1)!=null) 				t1.getNeighbour(b1).updateNeighbour(t1, dacp);
		if(t2.getNeighbour(t2.findpoint(a))!=null) 	t2.getNeighbour(t2.findpoint(a)).updateNeighbour(t2, dbcp);
		if(t2.getNeighbour(t2.findpoint(b))!=null) 	t2.getNeighbour(t2.findpoint(b)).updateNeighbour(t2, dacp);
		if(t3.getNeighbour(t3.findpoint(a))!=null) 	t3.getNeighbour(t3.findpoint(a)).updateNeighbour(t3, dbep);
		if(t3.getNeighbour(t3.findpoint(b))!=null) 	t3.getNeighbour(t3.findpoint(b)).updateNeighbour(t3, daep);
		if(t4.getNeighbour(t4.findpoint(a))!=null) 	t4.getNeighbour(t4.findpoint(a)).updateNeighbour(t4, dbep);
		if(t4.getNeighbour(t4.findpoint(b))!=null) 	t4.getNeighbour(t4.findpoint(b)).updateNeighbour(t4, daep);

		t1.setModified(true);
		t2.setModified(true);
		t3.setModified(true);
		t4.setModified(true);
		
		flips.addFlip(3, dacp);
		flips.addFlip(3, dbcp);
		flips.addFlip(3, dbep);
		flips.addFlip(3, daep);
		
		flips.addTetrahedron(dacp);
		flips.addTetrahedron(dbcp);
		flips.addTetrahedron(dbep);
		flips.addTetrahedron(daep);
		
	}

	private int findlastid(CVertex q1, CVertex q2, CVertex q3, CTetrahedron t) {
		CVertex tmp;
		for(int i=0;i<4;i++){
			tmp = t.getPoint(i);
			if(tmp!=q1 && tmp != q2 && tmp !=q3){
				return i;
			}
		}
		System.out.println("Failure is not an option");
		return -1;
	}
}
