package ProGAL.proteins.belta;

import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Vector;
import ProGAL.proteins.belta.SecondaryStructure.SSSegment;

/**
 * Quick note: One important feature of this extension of the normal PDBFile is that all 
 * indices are shifted so they are consecutive, which is required for betatopologies and 
 * secondary structures to make sense and not crash. 
 * TODO: 
 * <ul>
 *   <li>Fetching secondary structures from other chains than the standard</li>
 *   <li>Comment and test</li>
 * </ul>
 * @author R.Fonseca
 */
public class PDBFile extends ProGAL.proteins.PDBFile {
	private static final long serialVersionUID = -8750143450211694975L;
	
	private List<HelixRecord> helixRecords = new ArrayList<HelixRecord>();
	private List<SheetRecord> sheetRecords = new ArrayList<SheetRecord>();
	private Map<Integer, Integer> shiftMap = new HashMap<Integer,Integer>();
	
	public PDBFile(String path) {
		super(path, true);
		
		for(PDBRecord rec: super.records){
			if(rec instanceof HelixRecord){ helixRecords.add((HelixRecord)rec); }
			if(rec instanceof SheetRecord){ sheetRecords.add((SheetRecord)rec); }
		}
		
		genShiftMap();
		
	}
	
	private void genShiftMap(){
		int index = 0;
		for(AtomRecord a: super.getAtomRecords()){
			if(a.atomType.equalsIgnoreCase("CA")){
				shiftMap.put(a.residueNumber, index);
				index++;
			}
		}
	}

	public int shift(int residueNumber){
		return shiftMap.get(residueNumber);
	}
	
	public AtomRecord getAtom(int res, String atomType){
		for(AtomRecord ar: getAtomRecords()){
			if(shift(ar.residueNumber)==res && ar.atomType.equalsIgnoreCase(atomType))
				return ar;
		}
		return null;
	}
	
	public SecondaryStructure getSecondaryStructure(){
		//Find length
		PrimaryStructure ps = new PrimaryStructure(super.getSequence());
		int length = ps.sequence.length();
		
		char[] types = new char[length];
		for(int i=0;i<length;i++) types[i] = 'C';
		
		for(HelixRecord h: helixRecords)	{
			for(int r=h.initSeqNum;r<=h.endSeqNum;r++){
				try{ types[shift(r)] = 'H'; }catch(NullPointerException exc){}
			}
		}
		for(SheetRecord s: sheetRecords)	{
			for(int r=s.initSeqNum;r<=s.endSeqNum;r++){
				try{ types[shift(r)] = 'E'; }catch(NullPointerException exc){}
			}
			try{//Ensure that two completely adjacent strands have a coil between them
				int start = shift(s.initSeqNum);
				if(start>0 && types[start-1]=='E') types[start-1] = 'C';
			}catch(NullPointerException exc){}
			try{
				int end = shift(s.endSeqNum);
				if(end>types.length-1 && types[end+1]=='E') types[end+1] = 'C';
			}catch(NullPointerException exc){}
		}
		return new SecondaryStructure(ps, new String(types));
	}

	public BetaTopology getBetaTopology(){
		return genBetaTopology(getSecondaryStructure());
	}
	
	private BetaTopology genBetaTopology(SecondaryStructure ssass){
		List<Point> caCoords = getCACoords();

		SSSegment[] strands = ssass.getStrands();
		boolean[][] pairs = new boolean[strands.length][strands.length];
		for(int s1 = 0;s1<strands.length;s1++){
			SSSegment sc1 = strands[s1];
			for(int s2=s1+1;s2<strands.length;s2++){
				SSSegment sc2 = strands[s2];
				double maxMin = maxMinDist(sc1, sc2, caCoords);
				//System.out.printf("Strand %d to %d has maxMin %f\n",s1,s2,maxMin);
				if(maxMin<6){//They are paired
					if(parallel(sc1,sc2,caCoords))	pairs[Math.max(s1, s2)][Math.min(s1, s2)] = true;
					else							pairs[Math.min(s1, s2)][Math.max(s1, s2)] = true;	
				}
			}
		}
		
		BetaTopology nlAss = fixupForLooseStrands(ssass, pairs);//new BetaTopology(ssass, pairs);
		return nlAss;
	}
	
	private BetaTopology fixupForLooseStrands(SecondaryStructure ssass, boolean[][] pairs){
		
		for(int i=0;i<pairs.length;i++){
			boolean hasPartner = false;
			for(int j=0;j<pairs.length;j++){
				if(pairs[i][j] || pairs[j][i]){ hasPartner = true; break; }
			}
			if(!hasPartner){
				SSSegment sc = ssass.getStrands()[i];
				char[] types = new char[ssass.primaryStructure.sequence.length()];
				for(int r=0;r<types.length;r++) 	types[r] = ssass.getType(r).toChar();
				for(int r=sc.start;r<sc.end;r++)	types[r] = 'C';
				
				return genBetaTopology(new SecondaryStructure(ssass.primaryStructure, new String(types)));
			}
		}
		return new BetaTopology(ssass, pairs);
	}
	

	/** Returns true iff strand1 and strand2 are parallel */
	private static boolean parallel(SSSegment strand1, SSSegment strand2, List<Point> caCoords){
		Vector v1 = caCoords.get(strand1.start).vectorTo(caCoords.get(strand1.end-1));
		Vector v2 = caCoords.get(strand2.start).vectorTo(caCoords.get(strand2.end-1));
		return v1.dot(v2)>0;

	}
	/** Returns the largest minimal distance among the <em>m</em> smallest distance pairs between the strands, 
	 * where <em>m</em> is the length of the shortest strand. */
	private static double maxMinDist(SSSegment strand1, SSSegment strand2, List<Point> caCoords){
		SSSegment sc1 = strand1;
		SSSegment sc2 = strand2;
		if(strand1.length>strand2.length){
			sc1 = strand2;
			sc2 = strand1;
		}
		TreeSet<Double> minDists = new TreeSet<Double>();
		for(int r1=sc1.start;r1<sc1.end;r1++){
			double minDist = Float.POSITIVE_INFINITY;
			for(int r2=sc2.start;r2<sc2.end;r2++){
				double dist = caCoords.get(r1).distanceSquared(caCoords.get(r2));
				if(dist<minDist) minDist = dist;
			}
			minDists.add(minDist);
		}
		double maxMin = 0;
		for(int i=0;i<Math.min(sc1.length,3);i++) {
			double min = minDists.pollFirst();
			if(maxMin<min) maxMin = min;
		}
		return Math.sqrt(maxMin);
	}
	
	public static void main(String[] args){
		PDBFile f = new PDBFile("/Users/ras/Downloads/2KJX.pdb");
		BetaTopology bt = f.getBetaTopology();
		System.out.println(bt.secondaryStructure.primaryStructure);
		System.out.println(bt.secondaryStructure);
		System.out.println(bt);
	}
}
