package ProGAL.geom3d.complex.delaunayComplex;

import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.predicates.Predicates;


class Walk {
	private Predicates primitives;
	
	public Walk(Predicates primitives) {
		this.primitives=primitives;
	}


	public CTetrahedron walk(CTetrahedron t, CVertex p){
		int i;
		boolean next = false;
		
		while(true){
			for(i=0;i<4;i++){
				
				if(primitives.inplane(t.getPoint((i+1)%4),t.getPoint((i+2)%4), t.getPoint((i+3)%4), p)){
					p.setDegenerate(true);				
					p.setDegCase(CVertex.DegenerateCase.ONFACE);
					p.setDegPointOpposite(t.getPoint(i));


					if(primitives.inplane(t.getPoint((i+1)%4),t.getPoint((i+2)%4), t.getPoint(i%4), p)){
						p.setDegCase(CVertex.DegenerateCase.ONEDGE);
						p.setDegPointA(t.getPoint((i+1)%4));
						p.setDegPointB(t.getPoint((i+2)%4));
					}

					if(primitives.inplane(t.getPoint((i+2)%4),t.getPoint((i+3)%4), t.getPoint(i%4), p)){
						p.setDegCase(CVertex.DegenerateCase.ONEDGE);
						p.setDegPointA(t.getPoint((i+2)%4));
						p.setDegPointB(t.getPoint((i+3)%4));
					}

					if(primitives.inplane(t.getPoint((i+1)%4),t.getPoint((i+3)%4), t.getPoint(i%4), p)){
						p.setDegCase(CVertex.DegenerateCase.ONEDGE);
						p.setDegPointA(t.getPoint((i+1)%4));
						p.setDegPointB(t.getPoint((i+3)%4));
					}
					return t;				
				}
				Predicates.PlaneConfig pc = primitives.diffsides(t.getPoint((i+1)%4),t.getPoint((i+2)%4),t.getPoint((i+3)%4),p,t.getPoint(i));
//				System.out.println("walk .. "+t);
				if(pc==Predicates.PlaneConfig.DIFF){
					next=true;
					t = t.getNeighbour(i);		
					break;
				}

			}
			if(next) 	next=false;
			else		return t;
		}


	}

}
