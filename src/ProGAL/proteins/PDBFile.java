package ProGAL.proteins;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.io.IOToolbox;
import ProGAL.io.WebIOToolbox;

/**
 * A class for reading and representing a PDB-file. Currently, all ATOM records are 
 * read as well as HELIX and SHEET records. In some files there are multiple suggestions 
 * for atom-placements. This is resolved by only reading the first suggested atom-placement.
 * 
 * The four-letter PDB id and chain id are guessed based on the filename.
 * 
 * TODO: Represent more than one chain as well as ligands and read all types of records. 
 * @author R. Fonseca
 */
public class PDBFile implements Serializable{
	private static final long serialVersionUID = 1791445213018199901L;
	
	private List<PDBHelix> helixRecords = new ArrayList<PDBHelix>();
	private List<PDBSheet> sheetRecords = new ArrayList<PDBSheet>();
	private List<PDBAtom> atomRecords = new ArrayList<PDBAtom>();
	public String pdbId;
	public char chain;

	/** Construct a PDB-file from the specified file-path */
	public PDBFile(String path){
		this(new File(path));
	}
	/** Construct a PDB-file from the specified file */
	public PDBFile(File f){
		List<PDBRecord> records = readPDBFile(f,false);
		for(PDBRecord r: records){
			if(r instanceof PDBAtom) atomRecords.add((PDBAtom)r);
			else if(r instanceof PDBSheet) sheetRecords.add((PDBSheet)r);
			else if(r instanceof PDBHelix) helixRecords.add((PDBHelix)r);
		}
		pdbId = f.getName().substring(0,4);
		chain = f.getName().charAt(f.getName().length()-5);
	}

	/** Read a PDB-file from www.pdb.org.*/
	public static PDBFile downloadPDBFile(String pdbId){
		return new PDBFile(WebIOToolbox.downloadFile("http://www.pdb.org/pdb/files/"+pdbId.toUpperCase()+".pdb"));
	}
	
	public List<PDBAtom> getAtomRecords(){
		return atomRecords;
	}
	
	public List<Point> getCACoords(){
		List<Point> caCoords = new ArrayList<Point>();
		for(PDBAtom a: atomRecords) 
			if(a.atomType.equalsIgnoreCase("CA"))
				caCoords.add(a.coords);
		return caCoords;
	}
	
	public List<Point> getAtomCoords(){
		List<Point> coords = new ArrayList<Point>();
		for(PDBAtom a: atomRecords) 
			coords.add(a.coords);
		return coords;
	}

	public static List<PDBRecord> readPDBFile(File f, boolean onlySS){
		List<PDBRecord> ret = new ArrayList<PDBRecord>();
		for(String l: IOToolbox.readFromFile(f.getAbsolutePath()).split("\n")){
			try{
				if(!onlySS && l.startsWith("ATOM")) ret.add(new PDBAtom(l));
				if(l.startsWith("SHEET")) ret.add(new PDBSheet(l));
				if(l.startsWith("HELIX")) ret.add(new PDBHelix(l));
			}catch(Exception exc){}
		}
		
		//Filter out duplicate ATOM records
		for(int i=1;i<ret.size();i++){
			if(ret.get(i-1) instanceof PDBAtom && ret.get(i) instanceof PDBAtom){
				PDBAtom atom1 = (PDBAtom)ret.get(i-1);
				PDBAtom atom2 = (PDBAtom)ret.get(i);
				if(atom1.aaType.equals(atom2.aaType) && atom1.atomType.equals(atom2.atomType)){
					ret.remove(i);
					i--;
				}
			}
		}
		return ret;
	}

	public String getSequence(){
		StringBuilder sb = new StringBuilder();
		for(PDBAtom a: atomRecords){
			if(a.atomType.equalsIgnoreCase("CA")){
				sb.append(a.getSingleCharAAType());
			}
		}
		return sb.toString().toString();
	}
	
	public static interface PDBRecord{}

	public static class PDBHelix implements PDBRecord{
		public char initChainID, endChainID;
		public int initSeqNum, endSeqNum;

		public PDBHelix(String pdbLine){
			initChainID = pdbLine.charAt(19);
			initSeqNum = Integer.parseInt(pdbLine.substring(21,25).trim());
			endChainID = pdbLine.charAt(31);
			endSeqNum = Integer.parseInt(pdbLine.substring(33,37).trim());
		}
	}

	public static class PDBSheet implements PDBRecord{
		public int strand;
		public String sheetID;
		public int numStrands, sense;
		//public String initResName, endResName, curResName, prevResName;
		public int initSeqNum, endSeqNum, curResSeq, prevResSeq;
		public char initChainID, endChainID, curChainID, prevChainID;

		public PDBSheet(String pdbLine) throws Exception{
			//System.out.print("PDBSheet("+pdbLine+") .. ");
			try{
				strand = Integer.parseInt(pdbLine.substring(7,10).trim());
				sheetID = pdbLine.substring(11,14).trim();
				numStrands = Integer.parseInt(pdbLine.substring(14,16).trim());
				initChainID = pdbLine.charAt(21);
				initSeqNum = Integer.parseInt(pdbLine.substring(22,26).trim());
				endChainID = pdbLine.charAt(32);
				endSeqNum = Integer.parseInt(pdbLine.substring(33,37).trim());
				sense = Integer.parseInt(pdbLine.substring(38,40).trim());
				//initResName = pdbLine.substring(17,20).trim();
				//endResName = pdbLine.substring(28,31).trim();
				//if(sense!=0 && pdbLine.length()<50) System.err.println("Warning: '"+pdbLine+"' should contain more info");

				if(pdbLine.trim().length()>50){
					curChainID = pdbLine.charAt(49);
					curResSeq = Integer.parseInt(pdbLine.substring(50,54).trim());
					prevChainID = pdbLine.charAt(64);
					prevResSeq = Integer.parseInt(pdbLine.substring(65,69).trim());
				}
			}catch(Exception exc){
				//System.err.println("Warning: Somethings jiffy with: '"+pdbLine+"' (typically sense is not specified)");
				throw new RuntimeException("Warning: Somethings jiffy with: '"+pdbLine+"' (typically sense is not specified)");
			}
		}
	}

	public static class PDBAtom implements PDBRecord{
		public String atomType, aaType;
		public char chain;
		public int atomNumber, residueNumber;
		public Point coords;

		public PDBAtom(String pdbLine){
			atomNumber = Integer.parseInt(pdbLine.substring(5,11).trim());
			atomType = pdbLine.substring(13,16).trim();
			aaType = pdbLine.substring(17, 20);
			chain = pdbLine.charAt(21);
			residueNumber = Integer.parseInt(pdbLine.substring(22,26).trim());
			double x = Double.parseDouble(pdbLine.substring(28,38).trim());
			double y = Double.parseDouble(pdbLine.substring(38,46).trim());
			double z = Double.parseDouble(pdbLine.substring(46,54).trim());
			coords = new Point(x,y,z);
		}
		public String toString(){
			return String.format("PDBAtom[#:%d,res#:%d,type:%4s,coords:%s]",atomNumber, residueNumber, atomType, coords.toString());
		}
		
		public char getSingleCharAAType(){
			if(aaType.equalsIgnoreCase("ALA")) return 'A';  
			if(aaType.equalsIgnoreCase("ARG")) return 'R';  
			if(aaType.equalsIgnoreCase("ASN")) return 'N';  
			if(aaType.equalsIgnoreCase("ASP")) return 'D';  
			if(aaType.equalsIgnoreCase("CYS")) return 'C';  
			if(aaType.equalsIgnoreCase("GLU")) return 'E';  
			if(aaType.equalsIgnoreCase("GLN")) return 'Q';  
			if(aaType.equalsIgnoreCase("GLY")) return 'G';  
			if(aaType.equalsIgnoreCase("HIS")) return 'H';  
			if(aaType.equalsIgnoreCase("ILE")) return 'I';  
			if(aaType.equalsIgnoreCase("LEU")) return 'L'; 
			if(aaType.equalsIgnoreCase("LYS")) return 'K'; 
			if(aaType.equalsIgnoreCase("MET")) return 'M'; 
			if(aaType.equalsIgnoreCase("PHE")) return 'F'; 
			if(aaType.equalsIgnoreCase("PRO")) return 'P'; 
			if(aaType.equalsIgnoreCase("SER")) return 'S'; 
			if(aaType.equalsIgnoreCase("THR")) return 'T'; 
			if(aaType.equalsIgnoreCase("TRP")) return 'W'; 
			if(aaType.equalsIgnoreCase("TYR")) return 'Y'; 
			if(aaType.equalsIgnoreCase("VAL")) return 'V'; 
			return '?';
//			throw new Error("Unknown amino acid type: "+aaType);
		}
	}

}
