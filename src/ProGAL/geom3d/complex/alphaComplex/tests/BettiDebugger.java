package ProGAL.geom3d.complex.alphaComplex.tests;

import java.awt.Color;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.List;

import ProGAL.dataStructures.UnionFind;
import ProGAL.geom3d.Simplex;
import ProGAL.geom3d.complex.CTetrahedron;
import ProGAL.geom3d.complex.CTriangle;
import ProGAL.geom3d.complex.alphaComplex.AlphaFiltration;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.proteins.PDBFile;

public class BettiDebugger {
	static int i;
	static J3DScene scene;
	static AlphaFiltration af;
	static UnionFind<CTetrahedron> uf;
	static CTetrahedron bigTet =  new CTetrahedron(null,null,null,null);
	
	static void iterate(){
		List<Simplex> simplices = af.getSimplices();
		Simplex s = simplices.get(i);
		System.out.print("simplex "+i+" ");
		if(s instanceof CTriangle){
			CTriangle t = (CTriangle)s;
			CTetrahedron n0 = t.getAdjacentTetrahedron(0);
			CTetrahedron n1 = t.getAdjacentTetrahedron(1);
			if(n0.containsBigPoint()) n0 = bigTet;
			if(n1.containsBigPoint()) n1 = bigTet;
			CTetrahedron s0 = uf.find(n0);
			CTetrahedron s1 = uf.find(n1);
			if(  s0!=s1  ){
				uf.union(n0, n1);
				if(Math.min(t.getAdjacentTetrahedron(0).circumRadius(), t.getAdjacentTetrahedron(1).circumRadius())<5)
				scene.addShape(new LSS(t.getAdjacentTetrahedron(0).incenter(), t.getAdjacentTetrahedron(1).incenter(), 0.05), Color.GRAY, 7);
				System.out.print("MARKED : ");
			}
			System.out.print("circumrad: "+t.circumradius()+" ");
//			scene.addShape(t, new Color(100,100,100,80));
		}
		System.out.println(s);
		i--;
		
	}
	
	
	public static void main(String[] args) {
		PDBFile f1 = new PDBFile("/Users/ras/Documents/Datasets/NMR_Xray1/Xray/12_1R69.pdb");
	
		af = new AlphaFiltration(f1.getAtomCoords());
		i = af.getSimplices().size()-1; 
		scene = J3DScene.createJ3DSceneInFrame();
		uf = new UnionFind<CTetrahedron>();
		
		System.out.println("Simulating betti2 calculation");
		scene.getCanvas().addKeyListener(new KeyAdapter(){
			public void keyReleased(KeyEvent e) {
				int skip = 0;
				if(e.getKeyCode()==KeyEvent.VK_SPACE) skip = 1;
				if(e.getKeyCode()==KeyEvent.VK_ENTER) skip = 100;

				while(skip-->0){
					iterate();
				}
				
			}
		});
	}
	

}
