package ProGAL.proteins.beltaStructure;

import java.awt.Color;
import java.util.List;
import java.util.Locale;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Matrix;
import ProGAL.proteins.PDBFile;
import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SheetTopology;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;

public class CanonicalSheet {
	public final SheetTopology sheetTopology;
	
	public CanonicalSheet(SheetTopology st){
		this.sheetTopology=st;
		strandCenters = new Point[sheetTopology.strands.size()];
		strandOrientations = new Vector[sheetTopology.strands.size()];
	}
	
	/** 
	 * Points representing the centers of strands in the topology. The order of center-points 
	 * follows the order of <code>sheetTopology.strands</code>.
	 */
	private Point[] strandCenters;

	/** 
	 * Vectors representing the orientation of strands in the topology. The order of orientation-
	 * vectors follows the order of <code>sheetTopology.strands</code>. 
	 */
	private Vector[] strandOrientations;
	
	private Vector sheetCenter;
	private Matrix sheetOrientation;

	public Point strandCenter(int strand){
		return sheetOrientation.multiply(strandCenters[strand]).addThis(sheetCenter);
	}
	public Vector strandOrientation(int strand){
		return sheetOrientation.multiply(strandOrientations[strand]);
	}
	
	
	public static CanonicalSheet createFromPDB(PDBFile pdbFile, SheetTopology st){
		CanonicalSheet ret = new CanonicalSheet(st);
		List<Point> caCoords = pdbFile.getCACoords();
		int N = st.strands.size();

		//Read strand center and directions from pdb
		for(int i=0;i<N;i++){
			SSSegment strand = ret.sheetTopology.secondaryStructure.getStrands()[st.strands.get(i)];
			Point startCoord = caCoords.get(strand.start);
			Point endCoord = caCoords.get(strand.end-1);
			ret.strandCenters[i] = Point.getMidpoint(startCoord, endCoord);
			ret.strandOrientations[i] = startCoord.vectorTo(endCoord).normalizeThis();
		}
		
		//Set the sheet center and orientation
		ret.sheetCenter = ret.strandCenters[N/2].toVector();
		Vector[] orientations = new Vector[3];
		orientations[1] = ret.strandOrientations[N/2].clone();
		if(N==2) 	orientations[0] = ret.strandCenters[  0  ].vectorTo(ret.strandCenters[  1  ]).normalizeThis();
		else		orientations[0] = ret.strandCenters[N/2-1].vectorTo(ret.strandCenters[N/2+1]).normalizeThis();
		
		//Set sheet orientation
		orientations[2] = orientations[0].crossThis(orientations[1]).normalizeThis();
		orientations[0] = orientations[1].cross(orientations[2]);
		ret.sheetOrientation = Matrix.createColumnMatrix(orientations[0],orientations[1],orientations[2]);
		Matrix orientationInv = ret.sheetOrientation.invert();
		
		
		//Normalize strand centers and orientations
		for(int i=0;i<N;i++) {
			ret.strandCenters[i].subtractThis(ret.sheetCenter);
			orientationInv.multiplyIn(ret.strandCenters[i]);
			orientationInv.multiplyIn(ret.strandOrientations[i]);
		}
		
		//Make uni-directed strand orientations
		for(int i=N/2+1;i<N;i++){
			if(ret.strandOrientations[i].dot(ret.strandOrientations[i-1])<0) 
				ret.strandOrientations[i].multiplyThis(-1);
		}
		for(int i=N/2-1;i>=0;i--){
			if(ret.strandOrientations[i].dot(ret.strandOrientations[i+1])<0) 
				ret.strandOrientations[i].multiplyThis(-1);
		}

		return ret;
	}
	public static void main(String[] args){
		Locale.setDefault(Locale.ENGLISH);
		SecondaryStructure ss = new SecondaryStructure(" EEEEEEE GGGHHHHHHHHHHHH   HHHHHHHHHT SEEEEEEE HHHHHHHHHHHHHHT EEEE ");
		BetaTopology bTop = new BetaTopology(ss);
		bTop.setPaired(0, 1);
		bTop.setPaired(0, 2);
		SheetTopology sTop = bTop.getSheets().get(0);
		PDBFile pdbFile = new PDBFile("/Users/ras/Downloads/1CTF.pdb", true);
		CanonicalSheet cs = CanonicalSheet.createFromPDB(pdbFile, sTop);
		System.out.println("Done");
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		Point prev = null;
		for(Point p: pdbFile.getCACoords()){
			if(prev!=null){
				Color col = Color.gray;
				if(ss.getSegmentContainingResidue(pdbFile.getCACoords().indexOf(p)).type==SSType.STRAND){
					if(ss.getSegmentContainingResidue(pdbFile.getCACoords().indexOf(prev)).type==SSType.STRAND)
						col = Color.yellow;
					else 
						col = new Color(150,150,70);
				}else{
					if(ss.getSegmentContainingResidue(pdbFile.getCACoords().indexOf(prev)).type==SSType.STRAND)
						col = new Color(150,150,70);
				}
				scene.addShape(new LSS(prev,p,0.2), col);
			}
			prev = p;
		}
		for(int s=0;s<cs.sheetTopology.strands.size();s++){
			scene.addShape(new Sphere(cs.strandCenter(s),0.8),Color.YELLOW);
			scene.addShape(new LSS(cs.strandCenter(s), cs.strandCenter(s).addThis(cs.strandOrientation(s).multiplyThis(3)),0.3), Color.YELLOW.darker());
		}
		scene.setAxisEnabled(true);
	}
}
