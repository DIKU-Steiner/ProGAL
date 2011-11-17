package ProGAL.proteins.chainTree;

import ProGAL.geom3d.Vector;
import ProGAL.geom3d.volumes.Volume;
import ProGAL.math.Matrix;

import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.PI;

class CTNode {
	final ChainTree loop;
	Volume boundingVolume;
	CTNode left, right, parent = null;
	Matrix transformation;
	int low, high;
	int height;
	
	CTNode(CTNode left, CTNode right){
		this.left = left;
		this.right = right;
		left.parent = this;
		right.parent = this;
		transformation = left.transformation.multiply(right.transformation);
		this.loop = left.loop;
		this.low = left.low;
		this.high = right.high;
		this.height = Math.max(left.height, right.height)+1;
	}
	
	CTNode(ChainTree loop, int bond){
		this.loop = loop;
		this.low = bond;
		this.high = bond;
		this.height = 0;
		updateTransformation();
	}
	
	public boolean isLeaf() { return height==0; }
	
	public boolean isInternal() { return height!=0; }
	
	/** 
	 * Recursively updates first this nodes transformation matrix and 
	 * then every node up to the root. O(lgn) running time.
	 */
	void updateTransformation(){
		if(isLeaf()){
//			System.out.printf("updateTransformation .. low: %d\n",low);
			double angle = loop.bondAngles[low];
			Vector x = new Vector(    cos(PI-angle),    sin(PI-angle),0);
			Vector y = new Vector(cos(3*PI/2-angle),sin(3*PI/2-angle),0);
			Vector z = new Vector(0,0,1);
			Vector d = new Vector(loop.bondLengths[low], 0, 0);

			Vector v = new Vector(1,0,0);
			v.rotateIn(x, loop.torsions[low]);
			v.rotateIn(y, loop.torsions[low]);
			v.rotateIn(z, loop.torsions[low]);
			
			this.transformation = Matrix.create4x4ColumnMatrix(x, y, z, d);
		}else{
			this.transformation = left.transformation.multiply(right.transformation);
		}
		
		if(parent!=null) parent.updateTransformation();
	}
	
	
}
