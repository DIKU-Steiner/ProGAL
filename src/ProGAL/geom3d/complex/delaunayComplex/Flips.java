package ProGAL.geom3d.complex.delaunayComplex;

import java.util.Stack;

import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CVertex.DegenerateCase;
import ProGAL.geom3d.predicates.Predicates;
import ProGAL.geom3d.predicates.Predicates.PlaneConfig;
import ProGAL.geom3d.predicates.Predicates.SphereConfig;

class Flips {
	enum ApexConfig { CONVEX, CONCAVE, COPLANAR };

	Predicates primitives;
	
	private Stack<CTetrahedron> tetrahedrastack;
	private Stack<Flip> flipstack;
	private final Flip23 f23;
	private final Flip32 f32;
	private final Flip44 f44;

	public Flips(Predicates primitives) {
		this.primitives=primitives;
		tetrahedrastack = new Stack<CTetrahedron>();
		flipstack = new Stack<Flip>();
		f23 = new Flip23(this);
		f32 = new Flip32(f23,this);
		f44 = new Flip44(f23, this);
	}
	
	public CTetrahedron fixDelaunay(){
		CTetrahedron next_t=null;
		
		while(!flipstack.empty()){
			
			Flip f= flipstack.pop();
			CTetrahedron t=f.getT();
			int pid = f.getId();
			CVertex p = t.getPoint(pid);
			CTetrahedron t2= t.getNeighbour(pid);

			if(!t.isModified()){				

				if(t2!=null && t.isFlat()){					

					if(p.getDegCase()==DegenerateCase.ONFACE){
						int did = apex(t,t2);
						next_t = f23.flip23(t,pid,did);
					}else if(p.getDegCase()==DegenerateCase.ONEDGE){

					}					
				}					

				else if(t2!=null){
					int did = apex(t,t2);
					CVertex d = t2.getPoint(did);

					if(primitives.insphere(t,d)==SphereConfig.INSIDE){
						ApexConfig flipcase = apexConfig(t, p, pid, t2, d);

						if(flipcase==ApexConfig.CONVEX){

							next_t = f23.flip23(t, pid, did);
						}
						else if(flipcase==ApexConfig.CONCAVE){

							if(f23.getT3()!=null){

								next_t = f32.flip32(t,t2,f23.getT3(), pid, did);

								f23.setT3(null);

							}						
						}
						else if(flipcase==ApexConfig.COPLANAR){

							if(config44(t,t2,p,d)){
								//HMM?
							}
							else {

							}

						}

					}

				}
			}
		}
		return next_t;
	}

	private boolean config44(CTetrahedron t1, CTetrahedron t2, CVertex p, CVertex d) {
		int cid = f23.getCid();

		CVertex c = t1.getPoint(cid);
		int c2id= t2.findpoint(c);

		CTetrahedron t3 = t2.getNeighbour(c2id);
		CTetrahedron t4 = t1.getNeighbour(cid);

		if(t3!=null && t4!=null){
			int d3id = t3.findpoint(d);

			if(t3.getNeighbour(d3id)==t4){
				f44.flip44(t1, t2, t3, t4, p, d);
				return true;
			}
		}		
		return false;
	}

	public void addFlip(int id, CTetrahedron t){
		flipstack.push(new Flip(id,t));
	}
	
	public void addTetrahedron(CTetrahedron t){
		tetrahedrastack.push(t);
	}

	public int apex(CTetrahedron t1, CTetrahedron t2){
		for(int i=0;i<4;i++){
			if(t2.getNeighbour(i)==t1){
				return i;
			}
		}
		System.out.println("Problemer med apex\n");
		return -1;
	}

	public ApexConfig apexConfig(CTetrahedron t1, CVertex p, int pid, CTetrahedron t2, CVertex d){
		boolean concave = false;
		PlaneConfig case1, case2, case3;


		if((case1=primitives.diffsides(p, t1.getPoint((pid+1)%4),t1.getPoint((pid+2)%4),t1.getPoint((pid+3)%4),d))==PlaneConfig.DIFF){
			f23.setT3(findthird(t1,t2, (pid+3)%4));
			f23.setA1((pid+1)%4);
			f23.setB1((pid+2)%4);
			concave = true;
			if(f23.getT3()!=null) return ApexConfig.CONCAVE;
		}
		if((case2=primitives.diffsides(p, t1.getPoint((pid+1)%4),t1.getPoint((pid+3)%4),t1.getPoint((pid+2)%4),d))==PlaneConfig.DIFF){
			f23.setT3(findthird(t1,t2, (pid+2)%4));
			f23.setA1((pid+1)%4);
			f23.setB1((pid+3)%4);
			concave = true;
			if(f23.getT3()!=null) return ApexConfig.CONCAVE;
		}		

		if((case3=primitives.diffsides(p, t1.getPoint((pid+2)%4),t1.getPoint((pid+3)%4),t1.getPoint((pid+1)%4),d))==PlaneConfig.DIFF){
			f23.setT3(findthird(t1,t2, (pid+1)%4));
			f23.setA1((pid+2)%4);
			f23.setB1((pid+3)%4);
			concave = true;
			if(f23.getT3()!=null) return ApexConfig.CONCAVE;
		}
		if(concave) return ApexConfig.CONCAVE;
		if(case1==PlaneConfig.COPLANAR){
			f23.setCid((pid+3)%4);
			f23.setA1((pid+1)%4);
			f23.setB1((pid+2)%4);

			return ApexConfig.COPLANAR;
		}
		if(case2==PlaneConfig.COPLANAR){			
			f23.setCid((pid+2)%4);
			f23.setA1((pid+1)%4);
			f23.setB1((pid+3)%4);

			return ApexConfig.COPLANAR;
		}
		if(case3==PlaneConfig.COPLANAR){
			f23.setCid((pid+1)%4);
			f23.setA1((pid+2)%4);
			f23.setB1((pid+3)%4);			
			return ApexConfig.COPLANAR;
		}
		else return ApexConfig.CONVEX;
	}

	public CTetrahedron findthird(CTetrahedron t1, CTetrahedron t2, int c1){
		int c2= findpoint(t2, t1.getPoint(c1));

		if(t1.getNeighbour(c1)==t2.getNeighbour(c2)){
			return t1.getNeighbour(c1);
		}

		return null;
	}

	public int findpoint(CTetrahedron t, CVertex p){
		for(int i = 0; i<4; i++){
			if(t.getPoint(i)==p) {
				return i;
			}
		}
		System.out.println("Problemer med findpoint");
		//never happends:
		return -1;
	}

	public Stack<CTetrahedron> getTetrahedrastack() {
		return tetrahedrastack;
	}
	
	public Flip23 getFlip23() { return f23; }
	public Flip32 getFlip32() { return f32; }
	public Stack<Flip> getFlipstack() { return flipstack; }

	public void setTetrahedrastack(Stack<CTetrahedron> tetrahedrastack) {
		this.tetrahedrastack = tetrahedrastack;
	}
	
}


