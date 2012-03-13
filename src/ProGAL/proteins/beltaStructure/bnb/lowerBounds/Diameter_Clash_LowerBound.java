package ProGAL.proteins.beltaStructure.bnb.lowerBounds;

import java.util.List;

import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.proteins.beltaStructure.bnb.BnBNode;
import ProGAL.proteins.beltaStructure.bnb.BnBSolver;
import ProGAL.proteins.structure.AminoAcidChain;
import ProGAL.proteins.structure.Atom;
import ProGAL.proteins.structure.CBond;

/** 
 * A lower bound that seeks to minimize the diameter of the structure while not 
 * allowing any clashes.
 * 
 * @author R.Fonseca
 */
public class Diameter_Clash_LowerBound implements LowerBound {
	@Override
	public double lowerBound(BnBNode n, BnBSolver solver) {
		AminoAcidChain chain = solver.getChain();
		List<Integer> residues = solver.getDefinedResidues(n);
		
		double clashSq = 2*2;
		
		double maxDist = 0;
		for(int r1: residues){
			for(int r2: residues){
				if(r2>r1-3) break;
				for(Atom a1: chain.aminoAcid(r1).atoms()){
					for(Atom a2: chain.aminoAcid(r2).atoms()){
						double d = a1.distanceSquared(a2);
						
						if(d<clashSq) {
							return Double.POSITIVE_INFINITY;
						}
						if(d>maxDist) maxDist = d;
					}
				}
			}
		} 
//		if(maxDist<3) {
//			J3DScene scene = J3DScene.createJ3DSceneInFrame();
//			for(Atom a: solver.getChain().atoms()){
//				scene.addShape(new Sphere(a, 0.5));
//			}
//			for(CBond b: solver.getChain().covalentBonds()){
//				scene.addShape(new LSS(b.a1,b.a2, 0.2));
//			}
//			System.out.println(residues);
//			throw new Error("This cant be right "+n);
//		}
		
		return maxDist;
	}

}
