package ProGAL.proteins.beltaStructure.bnb;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.PriorityQueue;
import java.util.SortedSet;
import java.util.TreeSet;

import ProGAL.proteins.belta.BetaTopology;
import ProGAL.proteins.belta.PrimaryStructure;
import ProGAL.proteins.belta.SSType;
import ProGAL.proteins.belta.SecondaryStructure;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;
import ProGAL.proteins.belta.SheetTopology;
import ProGAL.proteins.beltaStructure.bnb.lowerBounds.LowerBound;
import ProGAL.proteins.structure.AminoAcidChain;
import ProGAL.proteins.structure.generators.CABAminoAcidGenerator;

public class BnBSolver {
	private final PrimaryStructure primaryStructure;
	private final SecondaryStructure secondaryStructure;
	@SuppressWarnings("unused")
	private final BetaTopology betaTopology;
	private final AminoAcidChain chain;
	private final SortedSet<BnBNode> best = new TreeSet<BnBNode>();
	private final Branchable[] parts;
	private final List<Integer>[] partsDefined; 
	private LowerBound lowerBound;
	
	@SuppressWarnings("unchecked")
	public BnBSolver(BetaTopology bt, LowerBound lowerBound){
		this.betaTopology = bt;
		this.secondaryStructure = bt.secondaryStructure;
		this.primaryStructure = secondaryStructure.primaryStructure;
		this.chain = new AminoAcidChain(primaryStructure.sequence, new CABAminoAcidGenerator());
		this.lowerBound = lowerBound;
		
		int loopSegments = 0;
		for(SSSegment seg: secondaryStructure.segments){
			if( seg.type!=SSType.STRAND ) loopSegments++;
		}
		List<SheetTopology> sheets = bt.getSheets();
		parts = new Branchable[sheets.size()+loopSegments];
		partsDefined = new List[sheets.size()+loopSegments];
		int i;
		for(i=0;i<sheets.size();i++){
			parts[i] = new SheetStructure(1, bt.getSheets().get(i), chain);
			partsDefined[i] = parts[i].definedResidues();
			if(i>0) partsDefined[i].addAll(partsDefined[i-1]);
		}
		for(SSSegment seg: secondaryStructure.segments){
			if( seg.type!=SSType.STRAND ) {
				parts[i] = new SegmentStructure(8, 4, secondaryStructure, seg, chain);
				partsDefined[i] = parts[i].definedResidues();
				if(i>0) partsDefined[i].addAll(partsDefined[i-1]);
				Collections.sort(partsDefined[i]);
				i++;
			}
		}
	}
	
	public void run(){
		best.clear();
		
		double upperBound = Double.POSITIVE_INFINITY;
		PriorityQueue<BnBNode> nodePool = new PriorityQueue<BnBNode>(1000000);
		nodePool.add(new BnBNode(null, -1));
		
		while(!nodePool.isEmpty()){
			BnBNode n = nodePool.poll();
//			System.out.println(n);
			if(n.part == parts.length-1){
				if(n.lowerBound<upperBound)
					upperBound = n.lowerBound;
				
				if(best.size()<1000) {
					System.out.println("Added: "+n);
					best.add(n);
				}
				else{
					if(n.lowerBound<best.last().lowerBound){
						best.remove(best.last());
						best.add(n);
						System.out.println("Added: "+n);
					}
				}
			}else{
				if(n.lowerBound<upperBound){
//					branch(n, nodePool);
					nodePool.addAll(branch(n));
				}
			}
		}
	}

	private List<BnBNode> branch(BnBNode n){
		LinkedList<BnBNode> ret = new LinkedList<BnBNode>();
		for(int s=0;s<parts[n.part+1].getStructures();s++){
			BnBNode newNode = new BnBNode(n, s);
			updateChain(newNode);
			newNode.lowerBound = lowerBound.lowerBound(newNode, this);
			ret.add(newNode);
		}
		return ret;
	}
//	private void branch(BnBNode n, Queue<BnBNode> nodePool){
//		for(int s=0;s<parts[n.part+1].getStructures();s++){
//			BnBNode newNode = new BnBNode(n, s);
//			updateChain(newNode);
//			newNode.lowerBound = lowerBound.lowerBound(newNode, this);
//			nodePool.add(newNode);
//		}
//	}
	
	public void updateChain(BnBNode n){
		while(n.parent!=null){
			parts[n.part].setStructure(n.structure);
			n=n.parent;
		}
	}
	
	public SortedSet<BnBNode> getBest(){ return best; }

	public AminoAcidChain getChain() {
		return chain;
	}
	public List<Integer> getDefinedResidues(BnBNode n){
		return partsDefined[n.part];
	}
}
