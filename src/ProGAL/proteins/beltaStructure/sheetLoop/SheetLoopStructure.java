package ProGAL.proteins.beltaStructure.sheetLoop;

import java.util.ArrayList;
import java.util.List;

import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.PrimaryStructure;
import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.belta.SheetTopology;
import ProGAL.proteins.beltaStructure.loop.LoopStructure;
import ProGAL.proteins.beltaStructure.sheet.SheetAlignment;
import ProGAL.proteins.beltaStructure.sheet.SurfaceSheetStructure;
import ProGAL.proteins.structure.AminoAcidChain;
import ProGAL.proteins.structure.Atom;
import ProGAL.proteins.structure.generators.HeavyAtomAminoAcidGenerator;

public class SheetLoopStructure extends AminoAcidChain implements PartialStructure {

	public final PrimaryStructure primaryStructure;
	public final SecondaryStructure secondaryStructure;
	public final BetaTopology betaTopology;

	final List<SurfaceSheetStructure> sheetStructures = new ArrayList<SurfaceSheetStructure>();
	final List<LoopStructure> loopStructures = new ArrayList<LoopStructure>();
	
	
	public SheetLoopStructure(BetaTopology bTop){
		super(bTop.secondaryStructure.primaryStructure.sequence, new HeavyAtomAminoAcidGenerator());
		
		this.betaTopology = bTop;
		this.secondaryStructure = bTop.secondaryStructure;
		this.primaryStructure = bTop.secondaryStructure.primaryStructure;
		
		for(SheetTopology st: bTop.getSheets()){
			sheetStructures.add(new SurfaceSheetStructure(new SheetAlignment(st)));
		}
		
		SSSegment[] segments = secondaryStructure.segments;
		for(int s1 = 0;s1<segments.length;s1++){
			SSSegment seg1 = segments[s1];
			if(seg1.type==SSType.STRAND){
				for(int s2=s1+1;s2<segments.length;s2++){
					if(segments[s2].type==SSType.STRAND){
						Atom[] loopTarget = {
								super.atom(segments[s2].start,"N"),
								super.atom(segments[s2].start,"CA"),
								super.atom(segments[s2].start,"C"),
						};
						loopStructures.add(new LoopStructure(secondaryStructure, s1+1,s2-1,loopTarget));
						
						s1 = s2-1;
						break;
					}
				}
			}
		}
	}
	
	@Override
	public void updateAtoms(AminoAcidChain chain){
		for(SurfaceSheetStructure sss: sheetStructures) sss.updateAtoms(this);
		
		//Set firsttransforms
		SurfaceSheetStructure sss = sheetStructures.get(0);//TODO: Modify to support more sheets
		for(LoopStructure ls: loopStructures) {
			ls.setFirstTransform(sss.getLoopTransform(ls.segment1));
			ls.updateAtoms(this);
		}
	}
	
}
