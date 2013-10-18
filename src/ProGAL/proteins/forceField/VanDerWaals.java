package ProGAL.proteins.forceField;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import ProGAL.geom3d.kineticDelaunay.KineticDelaunayTessellation;
import ProGAL.geom3d.kineticDelaunay.KineticDelaunayTessellation.ProblemInstanceType;
import ProGAL.geom3d.kineticDelaunay.Vertex;
import ProGAL.geom3d.kineticDelaunay.Tet;

public class VanDerWaals {

	// van der Waals radii (rA) and van der Waals well depths (eA) from the AMBER force field
	double rN = 1.8240;                 double eN = 0.1700;
	double rC = 1.9080;                 double eC = 0.1094;
	double rO = 1.6612;                 double eO = 0.2100;
	double rH = 0.6000;                 double eH = 0.0157;
	double rP = 2.1000;                 double eP = 0.2000;
	double rS = 2.0000;                 double eS = 0.2500;
	
	double rNN = (rN + rN)/2;             double rNN2 = rNN*rNN;    
	double rCC = (rC + rC)/2; 			  double rCC2 = rCC*rCC;
	double rNC = (rN + rC)/2;			  double rNC2 = rNC*rNC;
	
	double eNN = Math.sqrt(eN*eN); 
	double eCC = Math.sqrt(eC*eC); 
	double eNC = Math.sqrt(eN*eC);
	
	double[][] r = {{rNN2,rNC2,rNC2},{rNC2,rCC2,rCC2},{rNC2,rCC2,rCC2}};
	double[][] e = {{eNN,eNC,eNC},{eNC,eCC,eCC},{eNC,eCC,eCC}};
	Vertex u, v;
	
	public double getVdWPotential(KineticDelaunayTessellation kDT) { return getVdWPotential(kDT, 10000000.0);}

	public double getVdWPotential(KineticDelaunayTessellation kDT, double cutoff) {
		double cutoff2 = cutoff*cutoff;
		double sqDist;
		
		double energy = 0.0;
		if (kDT.getInstanceType() == ProblemInstanceType.pdbNCC) {
			int n = kDT.getNrVertices();
			double rij2, rij6;
			int iType, jType;
			for (int i = 4; i < n; i++) {
				u = kDT.getVertex(i);
				iType = (i-4)%3;
				for (int j = i+4; j < n+4; j++ ) {
					v = kDT.getVertex(j);
					sqDist = u.distanceSquared(v);
					if (sqDist < cutoff2) {
						jType = (j-4)%3;
						rij2 = r[iType][jType]/sqDist;
						rij6 = rij2*rij2*rij2;
						energy += e[iType][jType]*rij6*(rij6 - 2);
					}
				}
			}
		}
		return energy;
	}
	
	public double getVdWPotentialDTVeryFast(KineticDelaunayTessellation kDT, int maxDepth) {
		int n = kDT.getNrVertices();
		double rmj2, rmj6;
		int m, j, mType, jType;
		Tet tet;
		HashSet<Vertex> processed = new HashSet<Vertex>();
		ArrayList<Vertex> reachedVertices = new ArrayList<Vertex>();
		ArrayList<Vertex> tetVertices = new ArrayList<Vertex>();
		for (int i = 0; i < 4; i++) {
			kDT.getVertex(i).setDepth(-1);
			processed.add(kDT.getVertex(i));
		}
		double energy = 0.0;
		
		for (int i = 4; i < n; i++) {
			tet = kDT.getVertex(i).getTet();
			tetVertices.clear();
			for (int k = 0; k < 4; k++) if (!processed.contains(tet.getCorner(k))) tetVertices.add(tet.getCorner(k));
			if (!tetVertices.isEmpty()) {
				reachedVertices = tet.breadthFirstVertices(maxDepth);
				for (Vertex v0 : tetVertices) {
					processed.add(v0);
					m = v0.getId();			
					mType = (m-4)%3;
					for (Vertex v : reachedVertices) {
						j = v.getId();
						if (m+4 <= j) {
							jType = (j-4)%3;
							rmj2 = r[mType][jType]/v0.distanceSquared(v);
							rmj6 = rmj2*rmj2*rmj2;
							energy += e[mType][jType]*rmj6*(rmj6 - 2);
						}
					}
				}
				for (Vertex v : reachedVertices) v.flag = false;
			}
		}
		return energy;
	}
		
	public double getVdWPotentialDTFast(KineticDelaunayTessellation kDT) {
		List<Vertex> level2Vertices = new ArrayList<Vertex>();
		List<Vertex> level3Vertices = new ArrayList<Vertex>();
		double energy = 0.0;
		
		int n = kDT.getNrVertices();
		double alpha = kDT.getAlpha();
		double alpha2 = 4*alpha*alpha;
		double rij2, rij6;
		int iType, jType;
		int i, j;
		
		long startAdjComp = System.nanoTime();
		for (i = 4; i < n; i++) kDT.getVertex(i).computeAdjacentVerticesFast(alpha2);
		System.out.printf(" Adjacency computation time in miliseconds %.2f\n", (System.nanoTime() - startAdjComp)/1000000.0);
		
		for (i = 4; i < n-4; i++) {
			Vertex v0 = kDT.getVertex(i);
			v0.setDepth(0);
			
			for (Vertex v1 : v0.getAdjacentVerticesFast()) v1.setDepth(1);
		
			for (Vertex v1 : v0.getAdjacentVerticesFast()) {
				for (Vertex v2 : v1.getAdjacentVerticesFast()) {
					if (v2.getDepth() == null) {
						level2Vertices.add(v2);				
						v2.setDepth(2);
					}
				}
			}
			
			for (Vertex v2 : level2Vertices) {
				for (Vertex v3 : v2.getAdjacentVerticesFast()) {
					if (v3.getDepth() == null) {
						level3Vertices.add(v3);				
						v3.setDepth(3);
					}
				}
			}

			
			v0.setDepth(null);
			for (Vertex v : v0.getAdjacentVerticesFast()) v.setDepth(null);
			for (Vertex v : level2Vertices) v.setDepth(null);
			for (Vertex v : level3Vertices) v.setDepth(null);
			iType = (i-4)%3;
			for (Vertex b : v0.getAdjacentVerticesFast()) {
				j = b.getId() - 4;
				if (i <= j) {
					jType = j%3;
					rij2 = r[iType][jType]/v0.distanceSquared(b);
					rij6 = rij2*rij2*rij2;
					energy += e[iType][jType]*rij6*(rij6 - 2);
				}
			}
			for (Vertex b : level2Vertices) {
				j = b.getId() - 4;
				if (i <= j) {
					jType = j%3;
					rij2 = r[iType][jType]/v0.distanceSquared(b);
					rij6 = rij2*rij2*rij2;
					energy += e[iType][jType]*rij6*(rij6 - 2);
				}
			}
			for (Vertex b : level3Vertices) {
				j = b.getId()-4;
				if (i <= j) {
					jType = j%3;
					rij2 = r[iType][jType]/v0.distanceSquared(b);
					rij6 = rij2*rij2*rij2;
					energy += e[iType][jType]*rij6*(rij6 - 2);
				}
			}
			level2Vertices.clear();
			level3Vertices.clear();
		}
		return energy;
	}

	public static void main(String[] args){
		long start = System.nanoTime();
		double alpha = 100000;
		KineticDelaunayTessellation kDT = new KineticDelaunayTessellation(alpha, ProblemInstanceType.pdbNCC);
		System.out.printf("Construction time of kDT in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);

		for (int i = 0; i < 4; i++) kDT.getVertex(i).setDepth(-1);
		System.out.println("Number of tetrahedra in DT: " + kDT.getTetrahedra().size());
		VanDerWaals vDW = new VanDerWaals();

		start = System.nanoTime();
		int cutoff = 10000000;
		double energy = vDW.getVdWPotential(kDT, cutoff);
		System.out.printf(" Simple cutoff " + cutoff + ", time in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);
		System.out.println("Energy = " + energy);
		
		start = System.nanoTime();
		double energyFastDT = vDW.getVdWPotentialDTFast(kDT);
		System.out.printf(" alpha complex (alpha = " + alpha + "), time in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);
		System.out.println("EnergyDTFast = " + energyFastDT);
		
		start = System.nanoTime();
		double energyVeryFastDT = vDW.getVdWPotentialDTVeryFast(kDT, 25);
		System.out.printf(" alpha complex (alpha = " + alpha + "), time in miliseconds %.2f\n", (System.nanoTime() - start)/1000000.0);
		System.out.println("EnergyDTVeryFast = " + energyVeryFastDT);

	}
	
	
}
