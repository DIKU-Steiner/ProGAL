package ProGAL.proteins.beltaStructure.bnb;

import java.awt.event.MouseEvent;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.surface.ParametricParaboloid;
import ProGAL.geom3d.viewer.ClickListener;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SheetTopology;
import ProGAL.proteins.beltaStructure.sheet.SheetAlignment;
import ProGAL.proteins.beltaStructure.sheet.SurfaceSheetStructure;
import ProGAL.proteins.structure.AminoAcidChain;
import ProGAL.proteins.structure.generators.CABAminoAcidGenerator;

public class SheetStructure implements Branchable{
	private final AminoAcidChain chain;
	final SheetTopology sheetTop;
	private final SurfaceSheetStructure sheetStruc;
	private int structures;
	
	public SheetStructure(int structures, SheetTopology st, AminoAcidChain chain){
		this.structures = structures;
		this.sheetTop = st;
		this.chain = chain;
		
		SheetAlignment sa = new SheetAlignment(st);
		this.sheetStruc = new SurfaceSheetStructure(sa);
	}

	static SheetStructure sstruc;static int s1 = 0;
	static J3DScene scene;
	public static void main(String[] args){
		
		
		AminoAcidChain chain = new AminoAcidChain(     "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", new CABAminoAcidGenerator());
		SecondaryStructure ss = new SecondaryStructure(" EEEEE   HHHHHHHH  EEEE    EEEE  ");
		BetaTopology bt = new BetaTopology(ss,new boolean[][]{
				new boolean[]{false, true,false},
				new boolean[]{false,false,false},
				new boolean[]{false, true,false}
				});
		SheetTopology st = bt.getSheets().get(0);
		sstruc = new SheetStructure(3, bt.getSheets().get(0), chain);
		
		sstruc.setStructure(0);
		
		scene = J3DScene.createJ3DSceneInFrame();
		for(int strandIdx: st.strands){
			SSSegment strand = ss.getStrands()[strandIdx];
			for(int r=strand.start;r<strand.end;r++){
				scene.addShape(new Sphere(chain.atom(r, 0),0.3));
				scene.addShape(new Sphere(chain.atom(r, 1),0.2));
			}
		}
		scene.addClickListener(new ClickListener() {
			public void shapeClicked(Shape shape, MouseEvent event) {
				if(shape==null){
					s1++;
					if(s1%3==0) s1 = 0;
					sstruc.setStructure(s1);
					scene.repaint();
					return;
				}
			}
		});
	}
	

	public List<Integer> definedResidues() {
		List<Integer> ret = new LinkedList<Integer>();
		for(int s: sheetTop.strands){
			SSSegment seg = sheetTop.secondaryStructure.getStrands()[s];
			for(int i=seg.start;i<seg.end;i++) ret.add(i);
		}
		return ret;
	}
	
	@Override
	public void setStructure(int s) {
//		System.out.println("SheetStructure.setStructure("+s);
		double minC = 0.02, maxC = 0.1;
		double dC = (maxC-minC)/structures;
		double c = minC+s*dC;
		if(structures%2!=0) c+=dC/2;
		
		((ParametricParaboloid)sheetStruc.getSurface()).setC(c);
		for(int strandIdx: sheetTop.strands){
			SSSegment strand = sheetTop.secondaryStructure.getStrands()[strandIdx];
			for(int r=strand.start;r<=strand.end;r++){
				
				Point n = sheetStruc.getAtomPosition(r, 0);
				Point ca = sheetStruc.getAtomPosition(r, 1);
				Point cp = sheetStruc.getAtomPosition(r, 2);
				Vector can = ca.vectorTo(n);
				Vector cac = ca.vectorTo(cp);
				Vector tmp = can.cross(cac).multiplyThis(0.54);
				tmp.addThis( can.multiplyThis(-0.65) );
				tmp.addThis( cac.multiplyThis(-0.53) );
				Point cb = ca.add( tmp  );
				chain.atom(r, 0).setCoord(ca);
				chain.atom(r, 1).setCoord(cb);
			}
			
		}
	}

	@Override
	public int getStructures() {
		return structures;
	}
}
