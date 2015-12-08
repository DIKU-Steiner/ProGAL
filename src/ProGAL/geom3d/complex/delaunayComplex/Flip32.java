package ProGAL.geom3d.complex.delaunayComplex;

import ProGAL.geom3d.complex.CTetrahedron;

class Flip32 {
	
	private Flips flips;
	private Flip23 f23;
	
	public Flip32(Flip23 f23, Flips flips) {
		this.f23=f23;
		this.flips=flips;
		
	}
	
	public CTetrahedron flip32(CTetrahedron t1, CTetrahedron t2, CTetrahedron t3, int pid, int did){
		int i;
		
		int c1=-1;
		
		int a1 = f23.getA1();
		int b1 = f23.getB1();		
		
		for(i=0;i<4;i++){
			if(i!=a1 && i!=pid && i!=b1){
				c1=i;
				break;
			}				
		}
		
		if(c1==-1){
			System.out.println("Problemer med c1");
		}
		
		CTetrahedron pacd = new CTetrahedron(t1.getPoint(pid), t1.getPoint(a1), t1.getPoint(c1), t2.getPoint(did));
		CTetrahedron pbcd = new CTetrahedron(t1.getPoint(pid), t1.getPoint(b1), t1.getPoint(c1), t2.getPoint(did));
		
		int a2= t2.findpoint(t1.getPoint(a1));
		int a3= t3.findpoint(t1.getPoint(a1));
		
		int b2 = t2.findpoint(t1.getPoint(b1));
		int b3 = t3.findpoint(t1.getPoint(b1));
		
		pacd.setNeighbour(0, t2.getNeighbour(b2));
		pacd.setNeighbour(1, pbcd);
		pacd.setNeighbour(2, t3.getNeighbour(b3));
		pacd.setNeighbour(3, t1.getNeighbour(b1));
		
		pbcd.setNeighbour(0, t2.getNeighbour(a2));
		pbcd.setNeighbour(1, pacd);
		pbcd.setNeighbour(2, t3.getNeighbour(a3));
		pbcd.setNeighbour(3, t1.getNeighbour(a1));
		
		if(t1.getNeighbour(a1)!=null) t1.getNeighbour(a1).updateNeighbour(t1, pbcd);
		if(t1.getNeighbour(b1)!=null) t1.getNeighbour(b1).updateNeighbour(t1, pacd);
		if(t2.getNeighbour(a2)!=null) t2.getNeighbour(a2).updateNeighbour(t2, pbcd);
		if(t2.getNeighbour(b2)!=null) t2.getNeighbour(b2).updateNeighbour(t2, pacd);
		if(t3.getNeighbour(a3)!=null) t3.getNeighbour(a3).updateNeighbour(t3, pbcd);
		if(t3.getNeighbour(b3)!=null) t3.getNeighbour(b3).updateNeighbour(t3, pacd);
		
		flips.addFlip(0, pacd);
		flips.addFlip(0, pbcd);
		
		flips.addTetrahedron(pacd);
		flips.addTetrahedron(pbcd);
		
		t1.setModified(true);
		t2.setModified(true);
		t3.setModified(true);
		
		return pacd;		
		
	}
	


}
