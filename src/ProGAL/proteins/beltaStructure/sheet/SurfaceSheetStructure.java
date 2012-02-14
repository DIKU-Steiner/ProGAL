package ProGAL.proteins.beltaStructure.sheet;

import java.util.Arrays;

import ProGAL.geom3d.Vector;
import ProGAL.geom3d.surface.ParametricParaboloid;
import ProGAL.geom3d.surface.ParametricSurface;
import ProGAL.math.Matrix;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.beltaStructure.sheetLoop.PartialStructure;
import ProGAL.proteins.structure.AminoAcid;
import ProGAL.proteins.structure.AminoAcidChain;

/**
 * A representation of the backbone-structure of a beta-sheet based on an arbitrary 
 * parametric surface. Given only a sheet-alignment (which also specifies sheet-topology 
 * and secondary structure) as input a surface is created and positions of N, CA, C, O 
 * and CB atoms can be retrieved using the <code>getAtomPosition</code>-method. The  
 * surface specifying the structure can be retrieved and modified. Subsequent calls to 
 * <code>getAtomPosition</code> will reflect this change. Changes to the sheet-alignment
 * or secondary structure will also be reflected in the structure.
 * 
 * @author R.Fonseca
 */
public class SurfaceSheetStructure implements PartialStructure{
	/** The sheet alignment of this structure. Also holds reference to sheet-topology and secondary 
	 * structure. Any change to the sheet alignment will be reflected in the structure of this sheet.*/
	public final SheetAlignment sheetAlignment;

	private ParametricSurface surface;
	private ProGAL.geom2d.Point centerPos;
	private ProGAL.geom2d.Vector rowDis;
	private ProGAL.geom2d.Vector colDis;
	private ProGAL.geom3d.Vector[] atomPos = new ProGAL.geom3d.Vector[5];
	private boolean flip = false;

	/**
	 * Constructs a sheet structure using the specified sheet alignment and a hyperbolic 
	 * paraboloid as the surface. 
	 */
	public SurfaceSheetStructure(SheetAlignment sa){
		this.sheetAlignment = sa;

		this.surface = new ParametricParaboloid(0.6,-0.6,0.06);

		this.centerPos = new ProGAL.geom2d.Point(0,0);
		this.rowDis = new ProGAL.geom2d.Vector(-0.4, 3.7);
		this.colDis = new ProGAL.geom2d.Vector(5.3, 0.8);

		this.atomPos[0] = new ProGAL.geom3d.Vector(-0.47,-0.97,-0.35);//N
		this.atomPos[1] = new ProGAL.geom3d.Vector( 0.10,0.43,-0.99); //CA
		this.atomPos[2] = new ProGAL.geom3d.Vector(-0.51,1.74,-0.26); //C
		this.atomPos[3] = new ProGAL.geom3d.Vector(-1.74,1.85,-0.07); //O
		this.atomPos[4] = new ProGAL.geom3d.Vector(-0.34,0.47,-2.46); //CB
	}
	
	/**
	 * Return a reference to the surface that specifies this sheets structure. Currently, only 
	 * a hyperbolic paraboloid is possible, so the object returned by this method can be cast to 
	 * a <code>ParametricParaboloid</code>. The surface can be changed an
	 * @return
	 */
	public ParametricSurface getSurface(){
		return surface;
	}

	/* Interprets the residue immediately after a strand as belonging to the strand. 
	 * this is very useful for the getLoopTransform method. */
	private static boolean strandContainsRes(SSSegment seg, int res){
		return res>=seg.start && res<=seg.end;
	}

	/**
	 * Return the specified atom-position. The result is affected by the state of the surface 
	 * (can be altered via the <code>getSurface()</code>-object) and by the state of the 
	 * <code>sheetAlignment</code>-object. 
	 * 
	 * Note that significant improvements can still be made to the speed of this method.  
	 * @param res The residue index, referring to the residue in the primary sequence. The result
	 * is not specified if the residue is not part of the sheet. 
	 * @param atom The backbone-atom within the specified residue. 0 corresponds to N, 1 to CA, 
	 * 2 to C, 3 to O and 4 to CB. 
	 * @return A 3D point of the atom center
	 */
	public ProGAL.geom3d.Point getAtomPosition(int res, int atom){
		//Vertical strands. 
		int[] grid = getGridPoint(res);
		ProGAL.geom2d.Point residueCenter = centerPos.add(rowDis.multiply(grid[0]).addThis(colDis.multiply(grid[1])));

		SSSegment[] strands = sheetAlignment.sTop.secondaryStructure.getStrands();
		int[] strandOrder = sheetAlignment.sTop.getStrandOrder();
		int[] strandOrientation = sheetAlignment.sTop.getStrandOrientation();
		int column = -1;
		for(int i=0;i<strandOrder.length;i++){
			if(strandContainsRes(strands[strandOrder[i]],res)) {
				column = i;
				break;
			}
		}
		
		boolean up = (strandOrientation[column]==1);
		boolean out = (grid[0]%2==0)^up;
		Vector surNor = surface.getNormal(residueCenter);

		ProGAL.geom2d.Point atomPoint = residueCenter;
		atomPoint.addThis(	colDis.normalize().multiplyThis(atomPos[atom].x()*(out?1:-1)*(flip?-1:1))	);
		atomPoint.addThis(	rowDis.normalize().multiplyThis(atomPos[atom].y()*( up?1:-1))	);
		
		ProGAL.geom3d.Point ret = surface.getPoint(atomPoint);
		ret.addThis(  surNor.multiplyThis(atomPos[atom].z()*(out==up?1:-1)) );

		return ret;
	}

	private int[] getGridPoint(int res){
		SSSegment[] strands = sheetAlignment.sTop.secondaryStructure.getStrands();
		int[] strandOrder = sheetAlignment.sTop.getStrandOrder();
		int[] strandOrientation = sheetAlignment.sTop.getStrandOrientation();

		int column = -1;
		for(int i=0;i<strandOrder.length;i++){
			if(strandContainsRes(strands[strandOrder[i]],res)) {
				column = i;
				break;
			}
		}
		if(column<0)
			throw new RuntimeException(String.format("Failed to find residue %d in strands: %s\n",res, Arrays.toString(strandOrder)));

		int row = res;
		for(int c=column;c!=strandOrder.length/2;){
			//Translate current row index into next
			if(c<strandOrder.length/2){
				if(strandOrientation[c]==strandOrientation[c+1]){//Parallel
					row = (row-strands[strandOrder[c]].start)+sheetAlignment.alignmentPairs[c];
				}else{ //Antiparallel
					row = sheetAlignment.alignmentPairs[c]-(row-strands[strandOrder[c]].start);
				}
				c++;
			} else {
				if(strandOrientation[c-1]==strandOrientation[c]){//Parallel
					row = (row-sheetAlignment.alignmentPairs[c-1])+strands[strandOrder[c-1]].start;
				}else{ //Antiparallel
					row = (sheetAlignment.alignmentPairs[c-1]-row)+strands[strandOrder[c-1]].start;
				}
				c--;
			}
		}
		row = row-strands[strandOrder[strandOrder.length/2]].midpoint();
		if(strandOrientation[strandOrder.length/2]==0)
			row = -row;
		column = column-strandOrder.length/2;
		return new int[]{row,column};
	}

	
	/** 
	 * TODO 
	 * @param seg The first segment of a loop following a strand in this sheet structure
	 * @return
	 */
	public Matrix getLoopTransform(SSSegment seg){
		ProGAL.geom3d.Point lastC = getAtomPosition(seg.start-1, 2);
		ProGAL.geom3d.Point nextN = getAtomPosition(seg.start, 0);
		ProGAL.geom3d.Point nextCA = getAtomPosition(seg.start, 1);
		
		ProGAL.geom3d.Vector x = nextN.vectorTo(nextCA).normalizeThis();
		ProGAL.geom3d.Vector z = lastC.vectorTo(nextN).crossThis(x).normalizeThis();
		ProGAL.geom3d.Vector y = z.cross(x);
		
		Matrix m = Matrix.create4x4ColumnMatrix(x, y, z, nextN.toVector());
		return m;
	}
	
	
	public String toString(){
		return String.format("SurfaceSheetStructure[strands: %s]",this.sheetAlignment.sTop.getStrandOrderString());
	}


	@Override
	public void updateAtoms(AminoAcidChain chain) {
		SSSegment[] strands = sheetAlignment.sTop.secondaryStructure.getStrands();
		for(int sIdx: sheetAlignment.sTop.strands){
			SSSegment strand = strands[sIdx];
			for(int r=strand.start;r<=strand.end;r++){
				int atoms = 5;
				AminoAcid aa = chain.aminoAcid(r);
				if(aa.type()=='G') atoms = 4;//Glycine do not have a CB-atom
				
				for(int a=0;a<atoms;a++){
					aa.atom(a).set(getAtomPosition(r, a));
				}
			}
		}
	}
}
