package ProGAL.proteins.beltaStructure.sheet;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import ProGAL.geom2d.Point;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.viewer.ClickListener;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Matrix;
import ProGAL.math.Randomization;
import static ProGAL.math.Randomization.randBetween;
import ProGAL.proteins.PDBFile.AtomRecord;
import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.PDBFile;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SheetTopology;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.beltaStructure.sheetLoop.PartialStructure;
import ProGAL.proteins.structure.AminoAcidChain;

public class SurfaceSheetStructure implements PartialStructure{
	public final SheetAlignment sheetAlignment;

	private ParametricSurface surface;
	private ProGAL.geom2d.Point centerPos;
	private ProGAL.geom2d.Vector rowDis;
	private ProGAL.geom2d.Vector colDis;
	private ProGAL.geom3d.Vector[] atomPos = new ProGAL.geom3d.Vector[4];
	private boolean flip = false;


	public SurfaceSheetStructure(SheetAlignment sa){
		this.sheetAlignment = sa;

		this.surface = new ParametricParaboloid(0.61,-0.7,0.06);

		this.centerPos = new ProGAL.geom2d.Point(0,0);
		this.rowDis = new ProGAL.geom2d.Vector(-0.4, 3.7);
		this.colDis = new ProGAL.geom2d.Vector(5.3, 1.1);

		this.atomPos[0] = new ProGAL.geom3d.Vector(-0.53,-1.08,-0.01); //N
		this.atomPos[1] = new ProGAL.geom3d.Vector( 0.03,0.28,-0.24); //CA
		this.atomPos[2] = new ProGAL.geom3d.Vector(-0.48,1.59,-0.09); //C
		this.atomPos[3] = new ProGAL.geom3d.Vector(-1.82,1.71,-0.06); //O
//		this.atomPos[0] = new ProGAL.geom3d.Vector( 0.3, -1.3, 0); //N
//		this.atomPos[1] = new ProGAL.geom3d.Vector( -0.7,    0, 0.4); //CA
//		this.atomPos[2] = new ProGAL.geom3d.Vector( 0.3,  1.3, 0); //C
//		this.atomPos[3] = new ProGAL.geom3d.Vector(1.5,  1.3, 0); //O
//		this.atomPos[0] = new ProGAL.geom3d.Vector( -0.3, -1.3, 0); //N
//		this.atomPos[1] = new ProGAL.geom3d.Vector( 0.7,    0, 1.4); //CA
//		this.atomPos[2] = new ProGAL.geom3d.Vector( -0.3,  1.3, 0); //C
//		this.atomPos[3] = new ProGAL.geom3d.Vector(-1.5,  1.3, 0); //O
		
	}
	
	public ParametricSurface getSurface(){
		return surface;
	}

	/** Also interprets the residue immediately after a strand as belonging to the strand. 
	 * this is very useful for the getLoopTransform method. */
	private static boolean strandContainsRes(SSSegment seg, int res){
		return res>=seg.start && res<=seg.end;
	}

	public ProGAL.geom3d.Point getAtomPosition(int res, int atom){
		int[] grid = getGridPoint(res);
		ProGAL.geom2d.Point residueCenter = centerPos.add(rowDis.multiply(grid[0]).addThis(colDis.multiply(grid[1])));

		SSSegment[] strands = sheetAlignment.sTop.secondaryStructure.getStrands();
		int[] strandOrder = sheetAlignment.sTop.getStrandOrder();
		int[] strandOrientation = sheetAlignment.sTop.getStrandOrientation();
		int column = -1;
//		SSSegment strand = null;
		for(int i=0;i<strandOrder.length;i++){
			if(strandContainsRes(strands[strandOrder[i]],res)) {
				column = i;
//				strand = strands[strandOrder[i]];
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
		ret.addThis(  surNor.multiplyThis(atomPos[atom].z()*(out?1:-1)) );

		return ret;
	}

	private int[] getGridPoint(int res){
		SSSegment[] strands = sheetAlignment.sTop.secondaryStructure.getStrands();
		int[] strandOrder = sheetAlignment.sTop.getStrandOrder();
		int[] strandOrientation = sheetAlignment.sTop.getStrandOrientation();

		int column = -1;
//		SSSegment strand = null;
		for(int i=0;i<strandOrder.length;i++){
			if(strandContainsRes(strands[strandOrder[i]],res)) {
				column = i;
//				strand = strands[strandOrder[i]];
				break;
			}
		}
		if(column<0)
			throw new RuntimeException(String.format("Failed to find residue %d in strands: %s\n",res, Arrays.toString(strandOrder)));

		int row = res;
		boolean up = true;
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
//	public String toString(){
//		return String.format("SurfaceSheetStructure[\n %s\n rowDis:%s, colDis:%s\n centerPos:%s\n atomPos:%s\n]",
//				surface.toString(),rowDis.toString(),colDis.toString(),centerPos.toString(), Arrays.toString(atomPos));
//	}

	public static void main(String[] args) {
		Locale.setDefault(Locale.ENGLISH);
		optimizeSheetStructure();
		//		test1();
	}
	public static void test1(){
		//		ProGAL.geom3d.Point p1 = new ProGAL.geom3d.Point(22.019, 36.131, 15.009);
		//		ProGAL.geom3d.Point p2 = new ProGAL.geom3d.Point(26.494, 31.098, 13.057);
		//		System.out.println(p1.distance(p2));
		//		
		//		p1 = new ProGAL.geom3d.Point(23.283, 26.878, 16.769);
		//		p2 = new ProGAL.geom3d.Point(20.494, 26.163, 20.158);
		//		System.out.println(p1.distance(p2));

		SecondaryStructure ss = new SecondaryStructure(" EEEEE EEEE EEEE EEEEEEEE EEEE ");
		BetaTopology bt = new BetaTopology(ss);
		bt.setPaired(4, 0);
		bt.setPaired(3, 4);
		bt.setPaired(2, 3);
		bt.setPaired(1, 2);
		SheetTopology st = bt.getSheets().get(0);
		SheetAlignment sa = new SheetAlignment(st);
		//		sa.setAligned(1,27);
		//		sa.setAligned(29,23);
		//		sa.setAligned(19, 12);
		//		sa.setAligned(12, 10);
		System.out.println(Arrays.toString(sa.alignmentPairs));
		SurfaceSheetStructure sheetStruc = new SurfaceSheetStructure(sa);
		int res = 0;
		//		System.out.println((res=4)+": "+Arrays.toString(sheetStruc.getGridPoint(res)));
		//		System.out.println((res=1)+": "+Arrays.toString(sheetStruc.getGridPoint(res)));
		//		System.out.println((res=26)+": "+Arrays.toString(sheetStruc.getGridPoint(res)));
		//		System.out.println((res=8)+": "+Arrays.toString(sheetStruc.getGridPoint(res)));
		//		System.out.println((res=12)+": "+Arrays.toString(sheetStruc.getGridPoint(res)));

		System.out.print("Creating scene .. ");
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		System.out.println("done");
		for(SSSegment strand: ss.getStrands()){
			for(int r=strand.start;r<strand.end;r++){
				for(int a=0;a<4;a++){

					ProGAL.geom3d.Point pos = sheetStruc.getAtomPosition(r, a);
					System.out.println(r+" "+a );
					Color col = Color.gray;
					if(a==0) col = Color.BLUE;
					if(a==3) col = Color.RED;
					scene.addShape(new Sphere(pos, 0.7), col);
					System.out.println("Added atom at "+pos);
				}
			}
		}
	}


	public static void optimizeSheetStructure(){
		//		PDBFile f = new PDBFile(PDBWebReader.downloadPDBFile("1TTA"));
		PDBFile f = new PDBFile("src/edu/belta/sheetStructure/regularStrand.pdb");
		BetaTopology bt = f.getBetaTopology();
		SheetTopology st = bt.getSheets().get(0);
		SheetAlignment sa = new SheetAlignment(st, f);
		sa.setAligned(47, 16);
		SurfaceSheetStructure sheetStruc = new SurfaceSheetStructure(sa);
		sheetStructure = sheetStruc;
		((ParametricParaboloid)sheetStruc.surface).setDisplacement(new Vector(30.717,  28.358,   4.553));

		J3DScene scene = J3DScene.createJ3DSceneInFrame();

		updater = new SheetUpdater(scene);
		for(SSSegment strand: bt.secondaryStructure.getStrands()){
			for(int r=strand.start;r<strand.end;r++){
				addAtom(f.getAtom(r,  "N"), scene);
				addAtom(f.getAtom(r,  "CA"), scene);
				addAtom(f.getAtom(r,  "C"), scene);
				addAtom(f.getAtom(r,  "O"), scene);
				realPoints.add(f.getAtom(r,  "N").coords);
				realPoints.add(f.getAtom(r, "CA").coords);
				realPoints.add(f.getAtom(r,  "C").coords);
				realPoints.add(f.getAtom(r,  "O").coords);
			}
		}
		
		JPanel p = new JPanel(new GridLayout(2,6));
		scene.frame.add(p,BorderLayout.SOUTH);
		p.add(optimizePlacementButton());
		p.add(optimizeSurfaceButton());
		p.add(optimizeGridButton());
		p.add(optimizeStructureButton());
		p.add(flipButton());
		p.add(printButton());
		
		p.add(surfacePanel());
		scene.frame.validate();
		scene.centerCamera();
		scene.autoZoom();
	}
	private static SurfaceSheetStructure sheetStructure;
	private static List<ProGAL.geom3d.Point> sheetPoints = new ArrayList<ProGAL.geom3d.Point>();
	private static List<ProGAL.geom3d.Point> realPoints = new ArrayList<ProGAL.geom3d.Point>();

	private static double getRMSD(){
		double sqSum = 0;
		for(int i=0;i<sheetPoints.size();i++){
			sqSum += sheetPoints.get(i).distanceSquared(realPoints.get(i));
		}

		return Math.sqrt(sqSum/sheetPoints.size());
	}

	private static SheetUpdater updater;

	private static class SheetUpdater{
		private final List<Sphere> atomShapes = new ArrayList<Sphere>();
		private J3DScene scene;
		private void updatePointsFromSheetStructure(){
			int c= 0;;
			for(SSSegment strand: sheetStructure.sheetAlignment.sTop.secondaryStructure.getStrands()){
				for(int r=strand.start;r<strand.end;r++){
					for(int a=0;a<4;a++){
						sheetPoints.get(c++).set(sheetStructure.getAtomPosition(r, a));
					}
				}
			}
			scene.repaint();
		}
		SheetUpdater(J3DScene scene){
			this.scene = scene;
			for(SSSegment strand: sheetStructure.sheetAlignment.sTop.secondaryStructure.getStrands()){
				for(int r=strand.start;r<strand.end;r++){
					for(int a=0;a<4;a++){
						ProGAL.geom3d.Point pos = sheetStructure.getAtomPosition(r, a).clone();
						sheetPoints.add(pos);
						Sphere sphere = new Sphere(pos,0.6);
						atomShapes.add(sphere);

						Color col = Color.gray;
						if(a==0) col = Color.BLUE;
						if(a==3) col = Color.RED;
						scene.addShape(sphere, col);

						class AtomClick implements ClickListener{
							int res; String atom;ProGAL.geom3d.Point p;
							AtomClick(ProGAL.geom3d.Point p, int res, int atom){
								this.res = res; this.p = p;
								switch(atom){ 
								case 0: this.atom = "N";break;
								case 1: this.atom = "CA";break;
								case 2: this.atom = "C";break;
								case 3: this.atom = "O";break;
								}
							}
							public void shapeClicked(Shape shape, MouseEvent event) {
								if(shape instanceof Sphere && ((Sphere)shape).getCenter()==p)
									System.out.printf("Residue %d, atom %s\n",res,atom);
							}
						}
						scene.addClickListener(new AtomClick(pos, r, a));
					}


				}
			}
		}
	}

	private static void addAtom(AtomRecord ar, J3DScene scene){
		Sphere s = new Sphere(ar.coords, 0.7);
		Color col = new Color(100,100,100,100);
		if(ar.element.equalsIgnoreCase("N")) col = new Color(0,0,250,100);
		if(ar.element.equalsIgnoreCase("O")) col = new Color(250,0,0,100);
		scene.addShape(s, col);

		class AtomClick implements ClickListener{
			AtomRecord ar;
			AtomClick(AtomRecord ar){
				this.ar = ar;
			}
			public void shapeClicked(Shape shape, MouseEvent event) {
				if(shape instanceof Sphere && ((Sphere)shape).getCenter()==ar.coords)
					System.out.println(ar);
			}

		}
		scene.addClickListener(new AtomClick(ar));
		
	}

	private static final int maxIterations = 300;

	private static JButton optimizePlacementButton(){
		JButton ret = new JButton("Placement");
		ret.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				for(int it=0;it<maxIterations;it++){
					ParametricSurface oldSurface = sheetStructure.surface.clone();
					double oldRMSD = getRMSD();
					double rand = Math.random();

					if(rand<0.4){//Permute displacement
						double delta = 1; 
						if(sheetStructure.surface instanceof ParametricParaboloid){
							ParametricParaboloid pp = (ParametricParaboloid)sheetStructure.surface;
							pp.setDisplacement(pp.getDisplacement().addThis(new Vector(
									Randomization.randBetween(-delta,delta),
									Randomization.randBetween(-delta,delta),
									Randomization.randBetween(-delta,delta)
									)));
						}
					}else{ 

						if(sheetStructure.surface instanceof ParametricParaboloid){
							ParametricParaboloid pp = (ParametricParaboloid)sheetStructure.surface;
							Matrix m = pp.getRotation();
							int i = Randomization.randBetween(0, 3);
							Vector v = new Vector(m.getColumn(i));
							double delta = 0.4;
							v.addThis(new Vector(
									Randomization.randBetween(-delta,delta),
									Randomization.randBetween(-delta,delta),
									Randomization.randBetween(-delta,delta)
									));
							v.normalizeThis();
							Vector v1=null,v2=null,v3=null;
							switch(i){
							case 0: 
								v1 = v;
								v2 = new Vector(m.getColumn(2)).crossThis(v1).normalizeThis();
								v3 = v1.cross(v2);break;
							case 1: 
								v1 = v.cross(new Vector(m.getColumn(2))).normalizeThis();
								v2 = v;
								v3 = v1.cross(v2);break;
							case 2: 
								v1 = new Vector(m.getColumn(1)).crossThis(v);
								v2 = v.cross(v1);
								v3 = v;break;
							}
							pp.setRotation(Matrix.createColumnMatrix(v1, v2, v3));
						}
					}
					updater.updatePointsFromSheetStructure();
					double newRMSD = getRMSD();
					if(oldRMSD<newRMSD){
						sheetStructure.surface = oldSurface;
						updater.updatePointsFromSheetStructure();
						double finalRMSD = getRMSD();
						if(finalRMSD!=oldRMSD) 
							throw new RuntimeException(oldRMSD+" "+newRMSD+" "+finalRMSD);
					}
				}
			}});


		return ret;
	}
	private static JButton optimizeSurfaceButton(){
		JButton ret = new JButton("Surface");
		ret.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				for(int it=0;it<maxIterations;it++){
					ParametricSurface oldSurface = sheetStructure.surface.clone();
					double oldRMSD = getRMSD();
					double rand = Math.random();

					double delta = 0.01; 
					if(sheetStructure.surface instanceof ParametricParaboloid){
						ParametricParaboloid pp = (ParametricParaboloid)sheetStructure.surface;
						if(rand<0.33)		pp.setA(pp.getA()+Randomization.randBetween(-delta, delta));
						else if(rand<0.66)	pp.setB(pp.getB()+Randomization.randBetween(-delta, delta));
						else				pp.setC(pp.getC()+Randomization.randBetween(-delta, delta));
					}
					updater.updatePointsFromSheetStructure();
					double newRMSD = getRMSD();
					if(oldRMSD<newRMSD){
						sheetStructure.surface = oldSurface;
						updater.updatePointsFromSheetStructure();
						double finalRMSD = getRMSD();
						if(finalRMSD!=oldRMSD) 
							throw new RuntimeException(oldRMSD+" "+newRMSD+" "+finalRMSD);
					}
				}
			}});


		return ret;
	}

	private static JButton optimizeGridButton(){
		JButton ret = new JButton("Grid");
		ret.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				for(int it=0;it<maxIterations;it++){
					ProGAL.geom2d.Vector oldRowDis = sheetStructure.rowDis.clone();
					ProGAL.geom2d.Vector oldColDis = sheetStructure.colDis.clone();
					double oldRMSD = getRMSD();
					double rand = Math.random();

					double delta = 0.1; 
					if(rand<0.5)		sheetStructure.rowDis.addThis(new ProGAL.geom2d.Vector(randBetween(-delta, delta), randBetween(-delta, delta)));
					else				sheetStructure.colDis.addThis(new ProGAL.geom2d.Vector(randBetween(-delta, delta), randBetween(-delta, delta)));
					
					updater.updatePointsFromSheetStructure();
					double newRMSD = getRMSD();
					if(oldRMSD<newRMSD){
						sheetStructure.rowDis = oldRowDis;
						sheetStructure.colDis = oldColDis;
						updater.updatePointsFromSheetStructure();
						double finalRMSD = getRMSD();
						if(finalRMSD!=oldRMSD) 
							throw new RuntimeException(oldRMSD+" "+newRMSD+" "+finalRMSD);
					}
				}
			}});

		return ret;
	}

	private static JButton optimizeStructureButton(){
		JButton ret = new JButton("Structure");
		ret.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				for(int it=0;it<maxIterations;it++){
					Vector[] oldPos = {
							sheetStructure.atomPos[0].clone(), 
							sheetStructure.atomPos[1].clone(), 
							sheetStructure.atomPos[2].clone(), 
							sheetStructure.atomPos[3].clone() };
					Point oldCenter = sheetStructure.centerPos.clone();
					
					double oldRMSD = getRMSD();
					double rand = Math.random();

					double delta = 0.1; 
					Vector randVec = new Vector(randBetween(-delta,delta), randBetween(-delta,delta), randBetween(-delta,delta));
					if(rand<0.2)		sheetStructure.atomPos[0].addThis(randVec);
					else if(rand<0.4)	sheetStructure.atomPos[1].addThis(randVec);
					else if(rand<0.6)	sheetStructure.atomPos[2].addThis(randVec);
					else if(rand<0.7)	sheetStructure.atomPos[3].addThis(randVec);
					else				sheetStructure.centerPos.addThis(new ProGAL.geom2d.Vector(randBetween(-delta,delta),randBetween(-delta,delta)));
					
					updater.updatePointsFromSheetStructure();
					double newRMSD = getRMSD();
					if(oldRMSD<newRMSD){
						sheetStructure.atomPos = oldPos;
						sheetStructure.centerPos = oldCenter;
						updater.updatePointsFromSheetStructure();
						double finalRMSD = getRMSD();
						if(finalRMSD!=oldRMSD) 
							throw new RuntimeException(oldRMSD+" "+newRMSD+" "+finalRMSD);
					}
				}
			}});

		return ret;
	}
	private static JButton flipButton(){
		JButton ret = new JButton("Flip");
		ret.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				
				sheetStructure.flip = !sheetStructure.flip;
				updater.updatePointsFromSheetStructure();
			}});

		return ret;
	}
	private static JButton printButton(){
		JButton ret = new JButton("Print");
		ret.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent arg0) {
				System.out.println(sheetStructure);
				double conv=180/Math.PI;
				System.out.printf("N-CA-C %.2f\n", Point.getAngle(sheetPoints.get(0), sheetPoints.get(1), sheetPoints.get(2))*conv);
				System.out.printf("CA-C-N %.2f\n", Point.getAngle(sheetPoints.get(1), sheetPoints.get(2), sheetPoints.get(4))*conv);
				System.out.printf("C-N-CA %.2f\n", Point.getAngle(sheetPoints.get(2), sheetPoints.get(4), sheetPoints.get(5))*conv);
				System.out.println();
				System.out.printf(" N-CA-C-N   (Psi) %.2f\n", ProGAL.geom3d.Point.getDihedralAngle(sheetPoints.get(0), sheetPoints.get(1), sheetPoints.get(2), sheetPoints.get(4))*conv);
				System.out.printf("CA-C-N-CA (Omega) %.2f\n", ProGAL.geom3d.Point.getDihedralAngle(sheetPoints.get(1), sheetPoints.get(2), sheetPoints.get(4), sheetPoints.get(5))*conv);
				System.out.printf(" C-N-CA-C   (Phi) %.2f\n", ProGAL.geom3d.Point.getDihedralAngle(sheetPoints.get(2), sheetPoints.get(4), sheetPoints.get(5), sheetPoints.get(6))*conv);
				
				
			}});

		return ret;
	}
	private static JPanel surfacePanel(){
		JSlider aSlider = new JSlider(-100, 100);
		aSlider.addChangeListener(new ChangeListener(){
			public void stateChanged(ChangeEvent arg0) {
				int val = ((JSlider)arg0.getSource()).getValue();
				double min = -1;
				double max = 1;
				double newVal = ((val+100)/200.0)*(max-min)+min;
				((ParametricParaboloid)sheetStructure.surface).setC(newVal);
				System.out.printf("a = %.3f\n",newVal);
				updater.updatePointsFromSheetStructure();
			}
		});	
		JSlider bSlider = new JSlider(-100, 100);
		bSlider.addChangeListener(new ChangeListener(){
			public void stateChanged(ChangeEvent arg0) {
				int val = ((JSlider)arg0.getSource()).getValue();
				double min = -1;
				double max = 1;
				double newVal = ((val+100)/200.0)*(max-min)+min;
				((ParametricParaboloid)sheetStructure.surface).setC(newVal);
				System.out.printf("b = %.3f\n",newVal);
				updater.updatePointsFromSheetStructure();
			}
		});
		JSlider cSlider = new JSlider(-100, 100);
		cSlider.addChangeListener(new ChangeListener(){
			public void stateChanged(ChangeEvent arg0) {
				int val = ((JSlider)arg0.getSource()).getValue();
				double min = -0.2;
				double max = 0.2;
				double newVal = ((val+100)/200.0)*(max-min)+min;
				((ParametricParaboloid)sheetStructure.surface).setC(newVal);
				System.out.printf("c = %.3f\n",newVal);
				updater.updatePointsFromSheetStructure();
			}
		});
		
		JPanel ret = new JPanel(new GridLayout(3,2));
		ret.add(new JLabel("a "));
		ret.add(aSlider);
		ret.add(new JLabel("b "));
		ret.add(bSlider);
		ret.add(new JLabel("c "));
		ret.add(cSlider);
		return ret;
	}


	@Override
	public void updateAtoms(AminoAcidChain chain) {
		SSSegment[] strands = sheetAlignment.sTop.secondaryStructure.getStrands();
		for(int sIdx: sheetAlignment.sTop.strands){
			SSSegment strand = strands[sIdx];
			for(int r=strand.start;r<strand.end;r++){
				for(int a=0;a<4;a++){
					chain.atom(r, a).set(getAtomPosition(r, a));
				}
			}
		}
	}
}
