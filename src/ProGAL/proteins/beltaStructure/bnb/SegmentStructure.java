package ProGAL.proteins.beltaStructure.bnb;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointList;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.ClickListener;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Matrix;
import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.structure.AminoAcidChain;
import ProGAL.proteins.structure.generators.CABAminoAcidGenerator;

public class SegmentStructure implements Branchable{
	private final AminoAcidChain chain;
	private final SecondaryStructure ss;
	public final SSSegment seg;
	private int directions, rotations;
	private final Point prevAtom;
	private final boolean backwards, last;
	/** First array indicates structure, second indicates atoms */
	private final Point[][] atomStructures;
	
	/** Construct a segment structure for the specified segment with <code>dirs</code> directions, 
	 * <code>rots</code> rotations
	 * @param dirs
	 * @param rots
	 * @param ss Secondary structure of entire chain
	 * @param seg Segment to determine the structure of
	 * @param chain The entire amino acid chain. Kept because the <code>setStructure</code> updates the
	 * chain atoms within the bounds of the segment.
	 */
	public SegmentStructure(int dirs, int rots, SecondaryStructure ss, SSSegment seg, AminoAcidChain chain){
		this.directions = dirs;
		this.rotations = rots;
		this.ss = ss;
		this.seg = seg;
		this.chain = chain;
		
		//Determine backwards
		boolean backwards = true;
		for(int s=seg.segmentIndex-1;s>=0;s--){
			if(ss.segments[s].type==SSType.STRAND) {backwards = false; break;}
		}
		this.backwards = backwards;
		
		//Determine last
		last = ss.segments[ss.segments.length-1]==seg;

		//Determine prevAtom
		if(backwards) 	this.prevAtom = chain.atom(seg.end, 0);
		else			this.prevAtom = chain.atom(seg.start, 0);
		
		//Set up structures
		if(seg.type==SSType.COIL) 
			atomStructures = buildCoilStructures();
		else	
			atomStructures = buildHelixStructures();
	}
	


	public List<Integer> definedResidues() {
		List<Integer> ret = new LinkedList<Integer>();
		for(int i=seg.start;i<seg.end;i++) ret.add(i);
		return ret;
	}

	/** 
	 * Sets the structure of this segment. If <code>backward==false</code> this method places all calphas and 
	 * cbetas from (including) <code>seg.start</code> to (not including) <code>seg.end</code>. If 
	 * <code>backwards</code> 
	 */
	public void setStructure(int s) {
		Point[] points = atomStructures[s];
		Vector v;
		int start, end;
		if(backwards) 	{
			v =  points[points.length-1].vectorTo(prevAtom);
			start = seg.start;
		}else{
			v = points[0].vectorTo(prevAtom);
			start = seg.start+1;
		}
		end = last||backwards?seg.end-1:seg.end;
		int c=0;
//		System.out.printf("setStructure(%d) .. start: %d, end: %d, c: %d, %s\n",s,start,end,c, prevAtom);
		for(int r=start;r<=end;r++){
			chain.atom(r, 0).setCoord(points[c++]).addThis(v);
			chain.atom(r, 1).setCoord(points[c++]).addThis(v);
		}
	}

	/** Return the number of structures in this segment. Typically the number of rotations 
	 * times the number of directions. */
	public int getStructures() {
		return directions*rotations;
	}

	static SegmentStructure segstruc1; static int s1 = 0;
	static SegmentStructure segstruc2; static int s2 = 0;
	static SegmentStructure segstruc3; static int s3 = 0;
	static J3DScene scene;
	public static void main(String[] args){
		AminoAcidChain chain = new AminoAcidChain(     "AAAAAAAAAAAAAAAAAAAAAAAAAAA", new CABAminoAcidGenerator());
		SecondaryStructure ss = new SecondaryStructure(" EE    HHHHHHHHHHHHHH      ");

		segstruc1 = new SegmentStructure(10, 4, ss, ss.segments[2], chain);
		segstruc2 = new SegmentStructure(10, 4, ss, ss.segments[3], chain);
		segstruc3 = new SegmentStructure(10, 4, ss, ss.segments[4], chain);
		segstruc1.setStructure(s1);
		segstruc2.setStructure(s2);
		segstruc3.setStructure(s3);
		scene = J3DScene.createJ3DSceneInFrame();
		scene.addClickListener(new ClickListener() {
			public void shapeClicked(Shape shape, MouseEvent event) {
				if(shape==null){
					s1++;
					if(s3==40) {s3=0; s2++;}
					if(s2==40) {s2=0; s1++;}
					if(s1==40) {s1=0;}
					segstruc1.setStructure(s1);
					segstruc2.setStructure(s2);
					segstruc3.setStructure(s3);
					scene.repaint();
					return;
				}
				
				System.out.println(shape);
			}
		});
		SSSegment seg = ss.segments[2];
		for(int r=seg.start;r<seg.end;r++){
			System.out.println(chain.atom(r,0));
			scene.addShape(new Sphere(chain.atom(r, 0), 0.55), new Color(100,100,100,150));
		}
		System.out.println();
		seg = ss.segments[3];
		for(int r=seg.start;r<seg.end;r++){
			System.out.println(chain.atom(r,0));
			scene.addShape(new Sphere(chain.atom(r, 0), 0.55), new Color(200,00,00,150));
		}
		System.out.println();
		seg = ss.segments[4];
		for(int r=seg.start;r<seg.end;r++){
			System.out.println(chain.atom(r,0));
			scene.addShape(new Sphere(chain.atom(r, 0), 0.55), new Color(100,100,100,150));
		}
	}
	
	protected Point[][] buildHelixStructures(){
//		System.out.println("buildHelixStructures() "+(seg.length*2+2));
		Point[][] ret = new Point[getStructures()][seg.length*2+2];
		
		//Create an initial base helix that all others are (rotated) duplicates of.
		Point[] baseHelix = new Point[seg.length*2+2];
		double rad = 2.7, frequency = 2*Math.PI/3.6;
		for(int r=0;r<=seg.length;r++)
			baseHelix[2*r] = new Point(rad*Math.cos(frequency*r),rad*Math.sin(frequency*r),1.5*r);
		rad = 3.7;
		for(int r=0;r<=seg.length;r++)
			baseHelix[2*r+1] = new Point(rad*Math.cos(frequency*r),rad*Math.sin(frequency*r),1.5*r-1.0);
		
		//Put first atom at 0,0,0
		Vector v = baseHelix[0].toVector();
		for(int a=0;a<baseHelix.length;a++)
			baseHelix[a].subtractThis(v);
		
		//Rotate so direction from first to last CA is (0,0,1)
		Vector z = baseHelix[0].vectorTo(baseHelix[baseHelix.length-2]).normalizeThis();
		Vector x = new Vector(0,1,0).crossThis(z);
		Vector y = z.cross(x);
		Matrix rot = Matrix.createRowMatrix(x, y, z);
		for(int a=0;a<baseHelix.length;a++)
			rot.multiplyIn(baseHelix[a]);
//		System.out.printf("Base: First: %s, last: %s\n",baseHelix[0], baseHelix[baseHelix.length-2]);
		
//		J3DScene scene = J3DScene.createJ3DSceneInFrame();
//		scene.setAxisEnabled(true);
//		for(int a=2;a<baseHelix.length;a+=2){
//			scene.addShape(new LSS(baseHelix[a-2], baseHelix[a], 0.1));
//		}
//		Color col;
		
		//Now create all the other structures
		List<Point> dirs = PointList.generatePointsOnSphere(directions);
		for(int d=0;d<directions;d++){
			Vector dir = dirs.get(d).toVector();
			z = dir.normalize();
			x = new Vector(1,1,0).crossThis(z).normalizeThis();
			y = z.cross(x);
			rot = Matrix.createRowMatrix(x, y, z);
			for(int r=0;r<rotations;r++){
				double theta = (2*Math.PI*r)/rotations;
				Matrix rot1 = new Matrix(new double[][]{
						new double[]{Math.cos(theta), -Math.sin(theta),0},
						new double[]{Math.sin(theta),  Math.cos(theta),0},
						new double[]{0,0,1},
				});
				rot1 = rot.multiply(rot1);
//				System.out.printf("%d %d %d x:%s y:%s z:%s\n",d,r,d*rotations+r, x,y,z);
//				col = new Color(Randomization.randBetween(0, 250),Randomization.randBetween(0, 250),Randomization.randBetween(0, 250));
				for(int a=0;a<seg.length*2+2;a++){
					ret[d*rotations+r][a] = rot1.multiply(baseHelix[a]);
//					System.out.printf("> [%d][%d]: %s %s\n",d*rotations+r, a, ret[d*rotations+r][a], baseHelix[a]);
//					if(d==0 && a%2==0 && a>0){
//						scene.addShape(new LSS(rot1.multiply(baseHelix[a-2]), rot1.multiply(baseHelix[a]),0.1), col);
//					}
				}
//				if(d==0){
//					System.out.printf("structure: First: %s, last: %s\n",ret[d*rotations+r][0], ret[d*rotations+r][baseHelix.length-2]);
//					System.out.println("rot1 * baseLast = retLast");
//					rot1.toConsole();
//					baseHelix[baseHelix.length-2].toConsole();
//					rot1.multiply(baseHelix[baseHelix.length-2]).toConsole();
//				}
			}
		}
		
		return ret;
	}

	protected Point[][] buildCoilStructures(){
//		Point[][] ret = new Point[getStructures()][seg.length*2];
//		return ret;
		return buildHelixStructures();
	}


	
}
