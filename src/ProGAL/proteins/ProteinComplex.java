package ProGAL.proteins;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.complex.CEdge;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CTriangle;
import ProGAL.geom3d.complex.CVertex;
import ProGAL.geom3d.complex.alphaComplex.AlphaComplex;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.proteins.PDBFile;
import ProGAL.proteins.PDBFile.AtomRecord;

public class ProteinComplex extends AlphaComplex {
	public final PDBFile pdbFile;
	private List<Hole> holes;
	private List<Cavity> cavities;


	public ProteinComplex(PDBFile f) {
		this(f,2.8);
	}

	public ProteinComplex(PDBFile f, double waterRad) {
		super(f.getAtomCoords(), waterRad);
		this.pdbFile = f;
	}

	public List<CEdge> getEdges() {
		List<CEdge> ret = new ArrayList<CEdge>();
		for(CEdge e: super.edges){
			if(getDepth(e)>=0) ret.add(e);
		}
		return ret;
	}

	public List<CTetrahedron> getTetrahedra() {
		List<CTetrahedron> ret = new ArrayList<CTetrahedron>();
		for(CTetrahedron t: super.tetrahedra){
			if(getDepth(t)>=0) ret.add(t);
		}
		return ret;
	}

	public List<CTriangle> getTriangles() {
		List<CTriangle> ret = new ArrayList<CTriangle>();
		for(CTriangle t: super.triangles){
			if(getDepth(t)>=0) ret.add(t);
		}
		return ret;
	}


	public CVertex getVertex(Point p) {
		for(CVertex v: vertices){
			if(v.distance(p)<0.001) return v;
		}
		return null;
	}
	public AtomRecord getAtom(Point p){
		for(AtomRecord ar: pdbFile.getAtomRecords()){
			if(ar.coords.distanceSquared(p)<0.001){
				return ar;
			}
		}
		return null;
	}
	public char getAtomType(Point p){
		return getAtom(p).atomType.charAt(0);
	}


	public boolean isBuried(CVertex v){
		Set<CTetrahedron> hull = getVertexHull(v);
		for(CTetrahedron t: hull){
			if(getDepth(t)<0) return false;
		}
		return true;
	}

	public boolean isBuried(AtomRecord ar){
		return isBuried(getVertex(ar.coords));
	}


	public List<Hole> getHoles(){
		if(holes==null){
			holes = getHoles(alpha);
		}
		return holes;
	}

	public List<Hole> getHoles(double probeRad){
		List<Hole> ret = new LinkedList<Hole>();
		for(CTetrahedron t: tetrahedra){
			if(getDepth(t)<=0) continue;
			if(getInAlpha(t)<probeRad) continue;
			boolean include = true;
			for(Hole h: ret){ if(h.tetrahedra.contains(t)) { include = false; break; } }
			if(!include) continue;

			ret.add(new Hole(t, probeRad));
		}
		return ret;
	}

	public class Hole{
		public final List<CTetrahedron> tetrahedra = new ArrayList<CTetrahedron>();

		public Hole(CTetrahedron representative, double probeRad){
			List<CTetrahedron> candidates = new LinkedList<CTetrahedron>();
			candidates.add(representative);

			while(!candidates.isEmpty()){
				CTetrahedron c = candidates.remove(0);
				if(getDepth(c)<=0) continue;
				if(c.circumRadius()<probeRad) continue;
				if(tetrahedra.contains(c)) continue;
				tetrahedra.add(c);

				for(int n=0;n<4;n++)
					if(c.getTriangle(n).circumradius()>=probeRad) candidates.add(c.getNeighbour(n));
			}
		}
		public Hole(CTetrahedron representative){
			this(representative, alpha);
		}

		public double averageDepth(){
			int sum = 0;
			for(CTetrahedron t: tetrahedra) sum+=getDepth(t);
			return (sum*1.0)/tetrahedra.size();
		}
		public double volume(){
			double sum = 0;
			for(CTetrahedron t: tetrahedra) sum+=t.getVolume();
			return sum;
		}
		public String toString(){
			return String.format("ProteinComplex.Hole[%d tetrahedra, volume=%.2f, avg depth=%.2f]", tetrahedra.size(), volume(), averageDepth());
		}
	}





	public List<Cavity> getCavities(){
		if(cavities==null){
			cavities = new LinkedList<Cavity>();
			for(CTetrahedron t: tetrahedra){
				if(getDepth(t)>=0) continue;
				if(getInAlpha(t)<alpha) continue;
				if(getInAlpha(t)>=Cavity.largeProbeRadius) continue;
				boolean include = true;
				for(Cavity v: cavities){ if(v.tetrahedrons.contains(t)) { include = false; break; } }
				if(!include) continue;

				cavities.add(new Cavity(t));
			}
		}
		return cavities;
	}

	public class Cavity{
		public final List<CTetrahedron> tetrahedrons = new ArrayList<CTetrahedron>();
		private final static double largeProbeRadius = 3000;

		public Cavity(CTetrahedron representative){
			List<CTetrahedron> candidates = new LinkedList<CTetrahedron>();
			candidates.add(representative);

			while(!candidates.isEmpty()){
				CTetrahedron c = candidates.remove(0);
				if(getDepth(c)>=0) continue;
				if(c.circumRadius()<alpha) continue;
				if(c.circumRadius()>largeProbeRadius) continue;
				if(tetrahedrons.contains(c)) continue;
				tetrahedrons.add(c);

				candidates.add(c.getNeighbour(0));
				candidates.add(c.getNeighbour(1));
				candidates.add(c.getNeighbour(2));
				candidates.add(c.getNeighbour(3));
			}
		}
		public double volume(){
			double sum = 0;
			for(CTetrahedron t: tetrahedrons) sum+=t.getVolume();
			return sum;
		}
		public String toString(){
			return String.format("ProteinComplex.Cavity[%d tetrahedra, volume=%.2f]", tetrahedrons.size(), volume());
		}
	}


	public static void main(String[] args){
//		PDBFile f = new PDBFile("/Users/ras/Documents/Datasets/CASP8Training/T0490/T0490.pdb");
		PDBFile f = new PDBFile("/Users/pawel/Downloads/1X0O.pdb");
		ProteinComplex pc = new ProteinComplex(f);
		System.out.println(pc.getTetrahedra().size());
				J3DScene scene = J3DScene.createJ3DSceneInFrame();
				List<AtomRecord> allAtoms = f.getAtomRecords();
				AtomRecord prev = null;
				for(AtomRecord ca: f.getCARecords()){
					if(prev!=null){
						scene.addShape(new LSS(prev.coords, ca.coords, 0.1),java.awt.Color.GRAY,5);
					}
					scene.addText(allAtoms.indexOf(ca)+"", ca.coords, 1.0f);
					
					prev = ca;
				}
				System.out.println("Finding shortest paths:");
		//		
		//		List<CVertex> vertices = pc.getVertices();
		//		System.out.printf("> %d-%d: %d\n",0,1,		pc.getCovalentBondDistance(vertices.get(0), vertices.get(1)));
		//		System.out.printf("> %d-%d: %d\n",0,10,		pc.getCovalentBondDistance(vertices.get(0), vertices.get(10)));
		//		System.out.printf("> %d-%d: %d\n",0,100,	pc.getCovalentBondDistance(vertices.get(0), vertices.get(100)));
		//		System.out.printf("> %d-%d: %d\n",50,100,	pc.getCovalentBondDistance(vertices.get(50), vertices.get(100)));
		//		System.out.printf("> %d-%d: %d\n",50,60,	pc.getCovalentBondDistance(vertices.get(50), vertices.get(60)));
	}


}
