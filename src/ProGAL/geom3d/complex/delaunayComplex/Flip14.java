package ProGAL.geom3d.complex.delaunayComplex;

import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.CTetrahedron;


class Flip14 {

	private Flips flips;
	
	
	
	public Flip14(Flips flips) {
		super();
		this.flips=flips;
	}



	public CTetrahedron flip14(CTetrahedron t, CVertex p){
		CVertex p0, p1, p2, p3;
		
		p0=t.getPoint(0);
		p1=t.getPoint(1);
		p2=t.getPoint(2);
		p3=t.getPoint(3);
		
		CTetrahedron p0p1p2p = new CTetrahedron(p0,p1,p2,p);
		CTetrahedron p0p1p3p = new CTetrahedron(p0,p1,p3,p);
		CTetrahedron p0p2p3p = new CTetrahedron(p0,p2,p3,p);
		CTetrahedron p1p2p3p = new CTetrahedron(p1,p2,p3,p);
		
		p0p1p2p.setNeighbour(0, p1p2p3p);
		p0p1p2p.setNeighbour(1, p0p2p3p);
		p0p1p2p.setNeighbour(2, p0p1p3p);
		p0p1p2p.setNeighbour(3, t.getNeighbour(3));
		
		p0p1p3p.setNeighbour(0, p1p2p3p);
		p0p1p3p.setNeighbour(1, p0p2p3p);
		p0p1p3p.setNeighbour(2, p0p1p2p);
		p0p1p3p.setNeighbour(3, t.getNeighbour(2));
				
		p0p2p3p.setNeighbour(0,p1p2p3p);
		p0p2p3p.setNeighbour(1,p0p1p3p);
		p0p2p3p.setNeighbour(2,p0p1p2p);
		p0p2p3p.setNeighbour(3,t.getNeighbour(1));
		
		p1p2p3p.setNeighbour(0, p0p2p3p);
		p1p2p3p.setNeighbour(1, p0p1p3p);
		p1p2p3p.setNeighbour(2, p0p1p2p);
		p1p2p3p.setNeighbour(3, t.getNeighbour(0));
		
		if(t.getNeighbour(0)!=null)	t.getNeighbour(0).updateNeighbour(t, p1p2p3p);
		if(t.getNeighbour(1)!=null) t.getNeighbour(1).updateNeighbour(t, p0p2p3p);
		if(t.getNeighbour(2)!=null) t.getNeighbour(2).updateNeighbour(t, p0p1p3p);
		if(t.getNeighbour(3)!=null) t.getNeighbour(3).updateNeighbour(t, p0p1p2p);
		
		if(p.isDegenerate()){
			
			if(p.getDegCase()==CVertex.DegenerateCase.ONFACE){
				CVertex index=p.getDegPointOpposite();
				
				if(index!=t.getPoint(0) && index!=t.getPoint(1) && index!=t.getPoint(2)){
					p0p1p2p.setFlat(true);
				}
				else {
					flips.addFlip(3, p0p1p2p);
				}				
				if(index!=t.getPoint(0) && index!=t.getPoint(1) && index!=t.getPoint(3)){
					p0p1p3p.setFlat(true);
				}
				else {
					flips.addFlip(3, p0p1p3p);
				}				
				if(index!=t.getPoint(0) && index!=t.getPoint(2) && index!=t.getPoint(3)){
					p0p2p3p.setFlat(true);
				}
				else {
					flips.addFlip(3, p0p2p3p);
				}				
				if(index!=t.getPoint(1) && index!=t.getPoint(2) && index!=t.getPoint(3)){
					p1p2p3p.setFlat(true);
				}
				else {
					flips.addFlip(3, p1p2p3p);
				}
				if(p0p1p2p.isFlat())	flips.addFlip(3, p0p1p2p);
				if(p0p1p3p.isFlat())	flips.addFlip(3, p0p1p3p);
				if(p0p2p3p.isFlat())	flips.addFlip(3, p0p2p3p);
				if(p1p2p3p.isFlat())	flips.addFlip(3, p1p2p3p);
			}
			
			if(p.getDegCase()==CVertex.DegenerateCase.ONEDGE){
				
			}
			
		}
		
		else{
			flips.addFlip(3, p0p1p2p);
			flips.addFlip(3, p0p1p3p);
			flips.addFlip(3, p0p2p3p);
			flips.addFlip(3, p1p2p3p);
		}
		
		flips.addTetrahedron(p0p1p2p);
		flips.addTetrahedron(p0p1p3p);
		flips.addTetrahedron(p0p2p3p);
		flips.addTetrahedron(p1p2p3p);

		t.setModified(true);
		
		return p1p2p3p;
		
	}
	
}
