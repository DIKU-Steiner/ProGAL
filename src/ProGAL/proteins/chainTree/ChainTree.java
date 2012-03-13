package ProGAL.proteins.chainTree;

import java.util.LinkedList;
import java.util.Queue;

import ProGAL.geom3d.Line;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.math.Constants;
import ProGAL.math.Matrix;
import ProGAL.proteins.belta.PrimaryStructure;

/** 
 * A chain tree representation of the structure of a protein sub-chain. 
 *  
 * @author P.Winter, R.Fonseca
 */
public class ChainTree {
	
	private Matrix firstTransformation;
	private CTNode root;
	private CTNode[] bonds;
	double[] torsions;
	double[] bondLengths;
	double[] bondAngles;
	boolean[] locked;
	
	/** Construct a chain tree representing the structure of the specified sequence.
	 * @param ps Primary structure specifying the sequence
	 */
	public ChainTree(PrimaryStructure ps){
		this(ps,0,ps.sequence.length());
	}
	
	
	/** 
	 * Construct a chain tree representing the structure of the residues in the specified interval.
	 * 
	 * @param ps Primary structure specifying the sequence
	 * @param firstRes The residue-index of the first residue in the loop.
	 * @param lastRes The residue-index of the last residue in the loop.
	 */
	public ChainTree(PrimaryStructure ps, int firstRes, int lastRes){
		this.firstTransformation = Matrix.createIdentityMatrix(4);
		
		int residues = Math.max(firstRes, lastRes)-Math.min(firstRes, lastRes)+1;
		int length = residues*3;

		//Initialize arrays
		torsions = new double[length];
		bondLengths = new double[length];
		bondAngles = new double[length];
		bonds = new CTNode[length];
		locked = new boolean[length];
		for(int i=0;i<bonds.length;i++) {
			torsions[i] = Math.PI;
			bondLengths[i] = 1.7;
			if(i<bonds.length-1) bondAngles[i] = 2*Math.PI/3;
			if(i%3==2) locked[i] = true;
		}
		
		this.root = buildChainTree();
	}
	
	/** Used only once in constructor to build internal nodes */
	private CTNode buildChainTree(){
		Queue<CTNode> queue = new LinkedList<CTNode>();
		//Build leaves
		for(int i=0;i<bonds.length;i++){
			bonds[i] = new CTNode(this, i);
			queue.add(bonds[i]);
		}
		
		//Build balanced tree
		int sz = bonds.length, sz2;
		while (sz > 1) {
			sz2 = sz/2;
			for (int k = 0; k < sz2; k++)
				queue.add(  new CTNode(queue.poll(), queue.poll())  );
			if (2*sz2 != sz) {
				queue.add(queue.poll());
				sz2++; 
			}
			sz = sz2;
		}
		return queue.poll();
	}
	
	private Matrix getTransformation(int i, int j){
		//This is the dumb way
//		Matrix tmp = i==0?firstTransformation.clone():Matrix.createIdentityMatrix(4);
//		for(int c=i;c<j;c++){
//			tmp.multiplyThis(bonds[c].transformation);
//		}
//		if(true) return tmp;

		//Clever tree traversal
		CTNode split = root;
		while( split.left != null ){
			if (j-1 <= split.left.high) split = split.left;
			else {
				if (split.right.low <= i) split = split.right;
				else break;	
			}
		}
		if ((split.low == i) && (split.high == j-1)) {
			if(i==0)	return firstTransformation.multiply(split.transformation);
			else 		return split.transformation;
		}
		CTNode nd = split;
		Matrix rotMatrix = Matrix.createIdentityMatrix(4);
		
		// left subtree
		nd = nd.left;
		while (nd != null) {
			if (nd.low == i) { rotMatrix = nd.transformation.multiply(rotMatrix);nd = null; }
			else {
				if (i <= nd.left.high) { rotMatrix = nd.right.transformation.multiply(rotMatrix);nd = nd.left; }
				else nd = nd.right;
			}
		}
		// right subtree
		nd = split.right;
		while (nd != null) {
			if (nd.high == j-1) { rotMatrix = rotMatrix.multiply(nd.transformation);nd = null; }
			else {
				if (nd.right.low <= j-1) { rotMatrix = rotMatrix.multiply(nd.left.transformation);nd = nd.right; } 
				else nd = nd.left;
			}
		}

		if(i==0) 	return firstTransformation.multiply(rotMatrix);
		else		return rotMatrix;
	}

	private static Point getPos(Matrix m){				return new Point(m.get(0, 3), m.get(1, 3), m.get(2, 3));	}
	private static Vector getVec(Matrix m, int vec){	return new Vector(m.get(0, vec), m.get(1, vec), m.get(2, vec));	}
	
	private void ccdIteration(Point[] target, int i){
		Matrix trans_i = getTransformation(0, i);
		Line l = new Line(	getPos(trans_i),getVec(trans_i, 0));

		Matrix[] trans = new Matrix[target.length];
		trans[0] = getTransformation(0, bonds.length-target.length);
		for(int t=1;t<target.length;t++)
			trans[t] = trans[t-1].multiply(bonds[bonds.length-target.length+t].transformation);
		
		double nomSum = 0, denomSum = 0;
		for(int t=0;t<target.length;t++){
			Point M = getPos(trans[t]);
			Point O = l.orthogonalProjection(M);
			Vector f = O.vectorTo(target[t]); 
			Vector r = O.vectorTo(M);
			double rLen = r.length();
			if(rLen<Constants.EPSILON) continue;
			r.divideThis(rLen);
			Vector s = r.cross(l.getDir()).normalizeThis();
			nomSum += f.dot(s)*rLen;
			denomSum += f.dot(r)*rLen;
		}
        double tmp = Math.sqrt(  denomSum*denomSum+nomSum*nomSum  );
        double cosAlpha = denomSum/tmp;
        double sinAlpha = nomSum/tmp;
        double alpha = -Math.atan2(sinAlpha, cosAlpha);
		torsions[i] = torsions[i]+alpha;
		bonds[i].updateTransformation();
	}
	
	
	
	// ======== Public methods =========
	
	/** 
	 * Return the torsion angle of the specified residue and bond. This is a constant time 
	 * operation.
	 * @param res The residue-index in this chain-tree (0 is the first residue)
	 * @param bond Indication of which bond in the residue (0=phi, 1=psi and 2=omega).
	 * @return The torsion angle of the specified bond
	 */
	public double getTorsionAngle(int res, int bond){
		return torsions[res*3+bond];
	}
	
	/** 
	 * Changes the torsion angle of the specified residue and bond. Since the CT-hierarchy 
	 * must be traversed to the root this method runs O(lgn) in the number of residues. 
	 * @param res The residue-index in this chain-tree (0 is the first residue)
	 * @param bond Indication of which bond in the residue (0=phi, 1=psi and 2=omega).
	 * @param value The new value of the torsion angle
	 */
	public void setTorsionAngle(int res, int bond, double value){
		int bondIdx = res*3+bond;
		torsions[bondIdx] = value;
		bonds[bondIdx].updateTransformation();
	}
	
	public void setLocked(int res, int bond){
		locked[res*3+bond] = true;
	}
	
	public boolean isLocked(int res, int bond){
		return locked[res*3+bond];
	}
	
	/** 
	 * Change the orientation and translation of the entire chain. In practice this is done by changing 
	 * the transformation applied to the first leaf. 
	 * @param m The new first-transformation
	 */
	public void setFirstTransformation(Matrix m){
		if(m.getM()!=4 || m.getN()!=4) throw new RuntimeException("First transformation must be 4x4");
		this.firstTransformation = m;
	}
	
	/** 
	 * Return the position of a specific backbone atom. The internal structure of the chain tree is used 
	 * making the running-time O(lgn).
	 * @param res The residue index in this chain tree (0 is the first residue). 
	 * @param atom The atom index in the residue (0=N, 1=CA and 2=C).
	 * @return The position of the specified atom
	 */
	public Point getBackboneAtom(int res, int atom){
		Matrix trans = this.getTransformation(0, res*3+atom);
		return new Point(trans.get(0, 3), trans.get(1, 3), trans.get(2, 3));
	}
	
	/** 
	 * Collect all backbone atoms using the leaves only. Because all positions are being calculated, 
	 * the internal structure of the tree is not used. This gives a O(n)-time method instead of O(nlgn). 
	 */
	public Point[] getAllBackboneAtoms(){
		Point[] ret = new Point[bonds.length-2];
		Matrix trans = firstTransformation.clone();
		for(int i=0;i<bonds.length-2;i++){
			ret[i] = getPos(trans);
			trans.multiplyThis(bonds[i].transformation);
		}
//		for(int i=0;i<bonds.length-2;i++){
//			ret[i] = getPos(getTransformation(0, i));
//		}
		return ret;
	}
	
	/**
	 * Perform cyclic coordinate descent on this chain so the positions of the end-atoms overlap with 
	 * those specified in the <code>target</code>-array.  
	 * @param target The desired position of the <code>target.length</code> last atoms in the chain.
	 * @param iterations The number of iterations of CCD to perform
	 * @return The root-mean-square deviation of the last atoms from the target-points.
	 */
	public double closeCCD(Point[] target, int iterations){
		for(int it=0;it<iterations;it++){
			for(int i=0;i<bonds.length-1;i++){
				if(locked[i]) continue;
				ccdIteration(target, i);
			}
		}
		
		//Calculate rmsd from target
		double rms = 0;
		Matrix trans = getTransformation(0, bonds.length-target.length);
		rms += target[0].distanceSquared(getPos(trans));
		for(int t=1;t<target.length;t++){
			trans = trans.multiply(bonds[bonds.length-target.length+t].transformation);
			rms += target[t].distanceSquared(getPos(trans));
		}
		rms = Math.sqrt(rms/target.length);
		return rms;
	}
	
	
	
	public static void main(String[] args){
		PrimaryStructure ps = new PrimaryStructure("AAAAA");
		ChainTree loop = new ChainTree(ps,0,4);
		loop.setTorsionAngle(0, 0, 180*Math.PI/180);
		loop.setTorsionAngle(0, 1, 180*Math.PI/180);
		loop.setTorsionAngle(2, 0, 180*Math.PI/180);
		loop.setTorsionAngle(2, 1, 180*Math.PI/180);
		loop.setTorsionAngle(4, 0, 180*Math.PI/180);
		loop.setTorsionAngle(4, 1, 180*Math.PI/180);
//		loop.setTorsionAngle(0, 0, -60*Math.PI/180);
//		loop.setTorsionAngle(0, 1, -30*Math.PI/180);
//		loop.setTorsionAngle(2, 0, -60*Math.PI/180);
//		loop.setTorsionAngle(2, 1, -30*Math.PI/180);
//		loop.setTorsionAngle(4, 0, -60*Math.PI/180);
//		loop.setTorsionAngle(4, 1, -30*Math.PI/180);
		loop.getTransformation(0, 1).toConsole();
		Vector x = new Vector(2,1,0).normalizeThis();
		Vector y = new Vector(3,0,1).cross(x).normalizeThis();
		Vector z = x.cross(y);
		Vector d = new Vector(1,0,2);
		Matrix first = Matrix.create4x4ColumnMatrix(x, y, z, d);
		loop.setFirstTransformation(first);
		loop.getTransformation(0, 1).toConsole();
	}

}
