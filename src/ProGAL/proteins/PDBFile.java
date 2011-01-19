package ProGAL.proteins;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.io.IOToolbox;
import ProGAL.io.WebIOToolbox;

/**
 * A class for reading and representing a PDB-file. All lines are read but not all are parsed.
 * All the records can be retrieved as they are shown in the PDB-file. The ATOM-records, however, 
 * can also be accessed according to model and chain id. If a model-number and chain id are not 
 * specified the method will assume that the first chain of the first model is meant.  
 *  
 * In some files there are multiple suggestions for atom-placements. This is resolved by only 
 * reading the first suggested atom-placement.
 * 
 * The four-letter PDB id and chain id are guessed based on the filename.
 * 
 * TODO: Represent ligands, heterogens and secondary structure.
 * TODO: Create method that extracts resolution from REMARK-records 
 * TODO: Check the integrity of different models, e.g. that the sequences are the same.
 * @author R. Fonseca
 */
public class PDBFile implements Serializable{
	private static final long serialVersionUID = 1791445213018199901L;
	
	private boolean includeHydrogens = false;
	private final List<PDBRecord> records;
	private final List<PDBModel> models = new ArrayList<PDBModel>();
	
	/** The name of this PDB-file. Typically the PDB-id. */
	public String name;

	
	/** Construct a PDB-file from the specified file-path */
	public PDBFile(String path){
		this(new File(path));
	}
	
	/** Construct a PDB-file from the specified file */
	public PDBFile(File f){
		name = f.getName().replaceAll(".pdb", "");
		records = readPDBFile(f,false);
		buildPDBStructure();
	}

	private void buildPDBStructure(){
		
		//Build models and chains and fill them with ATOM-records
		for(int r=0;r<records.size();r++){
			if(records.get(r) instanceof AtomRecord) r = buildModel(r);
			else if(records.get(r) instanceof ModelRecord) r = buildModel(r);
			
		}
		
		//TODO: Build secondary structure
	}
	
	/**
	 * Build models and chains from the records.
	 * @return The index of the last line corresponding to this record. Ideally this is an ENDMDL record, but it can be 
	 * a TER or ATOM if the PDB-file is bad. 
	 */
	private int buildModel(int r){
		PDBModel model;
		if(records.get(r) instanceof ModelRecord) {
			model = new PDBModel( ((ModelRecord)records.get(r)).number );
			++r;
		}else{
			model = new PDBModel( models.size()>1?models.get(models.size()-1).number+1:0 );
		}
		
		//Now start reading the model
		for(;r<records.size();r++){
			PDBRecord record = records.get(r);
			if(record instanceof EndModelRecord) break;
			if(record instanceof AtomRecord) {
				r = buildChain(r, model);
			}else if(record instanceof HetatmRecord){
				//TODO: Do something about HETATM's
				//TODO: Do something about ANISOU's
			}else if(record instanceof ParentRecord){
				//TODO: Do something about PARENT's (ignore them perhaps?)
			}else{ break; }
		}
		if(model==null) throw new RuntimeException("Expected a model at "+name+":"+(r+1));
		models.add(model);

		return r;
	}
	
	/** Builds a chain consisting of ATOM-records. Assumes that 
	 * <code>(records.get(r) instanceof AtomRecord)</code>. 
	 * @return The last position of a TER or ATOM-record within this chain;  
	 */
	private int buildChain(int r, PDBModel m){
		PDBChain chain = null;
		for(;r<records.size();r++){
			PDBRecord record = records.get(r);
			if(record instanceof TerRecord) break;
			if(record instanceof AtomRecord){
				AtomRecord aRecord = (AtomRecord)record;
				if(chain==null) chain = new PDBChain(aRecord.chain);
				if(chain.chainId!=aRecord.chain) break;
				
				chain.atomRecords.add(aRecord);
			}else if(record instanceof HetatmRecord){
				//TODO: Do something about HETATM's
			}else if(record instanceof AnisouRecord){
				//TODO: Do something about ANISOU's
			}else{ break; }//throw new RuntimeException("Unexpected record at "+name+":"+(r+1)+"\n> "+record); }
		}
		if(chain==null) throw new RuntimeException("Expected a chain at "+name+":"+(r+1)+". Is model empty?");
		m.chains.add(chain);
		return r;
	}
	
	/** Read a PDB-file from www.pdb.org.*/
	public static PDBFile downloadPDBFile(String pdbId){
		return new PDBFile(WebIOToolbox.downloadFile("http://www.pdb.org/pdb/files/"+pdbId.toUpperCase()+".pdb"));
	}


	public List<PDBModel> getModels(){ return models; }

	/** Returns the ATOM-records. Only the records in the first chain of the first model are returned. */
	public List<AtomRecord> getAtomRecords(){
		return getAtomRecords(0,0);
	}

	/** Returns the ATOM-records of CA-atoms. Only the records in the first chain of the first model are returned. */
	public List<AtomRecord> getCARecords(){
		return getCARecords(0,0);
	}
	
	/** Returns the CA-coordinates of atoms in the first chain in the first model. */
	public List<Point> getCACoords(){
		return getCACoords(0,0);
	}

	/** Returns all atom-coordinates in the first chain in the first model. */
	public List<Point> getAtomCoords(){
		return getAtomCoords(0,0);
	}

	/** Returns the ATOM-records of the specified model and chain. */
	public List<AtomRecord> getAtomRecords(int modelNum, int chainNum){
		List<AtomRecord> ret = new ArrayList<AtomRecord>();
		for(AtomRecord ar: models.get(modelNum).chains.get(chainNum).atomRecords){
			if(!includeHydrogens && ar.isHydrogen()) continue;
			
			ret.add(ar);
		}
		
		return ret;
	}

	/** Returns the CA ATOM-records of the specified model and chain. */
	public List<AtomRecord> getCARecords(int modelNum, int chainNum){
		List<AtomRecord> ret = new ArrayList<AtomRecord>();
		for(AtomRecord ar: models.get(modelNum).chains.get(chainNum).atomRecords){
			if(ar.atomName.equalsIgnoreCase("CA")) 
				ret.add(ar);
		}
		
		return ret;
	}
	
	/** Returns the CA-coordinates of the specified model and chain. */
	public List<Point> getCACoords(int modelNum, int chainNum){
		List<Point> caCoords = new ArrayList<Point>();
		for(AtomRecord a: getAtomRecords(modelNum,chainNum)) 
			if(a.atomName.equalsIgnoreCase("CA"))
				caCoords.add(a.coords);
		return caCoords;
	}

	/** Returns all atom-coordinates of the specified model and chain. */
	public List<Point> getAtomCoords(int modelNum, int chainNum){
		List<Point> coords = new ArrayList<Point>();
		for(AtomRecord a: getAtomRecords(modelNum,chainNum)) {
			if(!includeHydrogens && a.isHydrogen()) continue;
			coords.add(a.coords);
		}
		return coords;
	}

	
	
	
	private static List<PDBRecord> readPDBFile(File f, boolean onlySS){
		List<PDBRecord> ret = new ArrayList<PDBRecord>();
		for(String l: IOToolbox.readFromFile(f.getAbsolutePath()).split("\n")){
			try{
				if(l.startsWith("ATOM")) ret.add(new AtomRecord(l));
				else if(l.startsWith("HETATM")) ret.add(new HetatmRecord(l));
				else if(l.startsWith("SHEET")) ret.add(new SheetRecord(l));
				else if(l.startsWith("HELIX")) ret.add(new HelixRecord(l));
				else if(l.startsWith("MODEL")) ret.add(new ModelRecord(l));
				else if(l.startsWith("ENDMDL")) ret.add(new EndModelRecord(l));
				else if(l.startsWith("TER")) ret.add(new TerRecord(l));
				else if(l.startsWith("REMARK")) ret.add(new RemarkRecord(l));
				else if(l.startsWith("PARENT")) ret.add(new ParentRecord(l));
				else if(l.startsWith("ANISOU")) ret.add(new AnisouRecord(l));
				else ret.add(new OtherRecord(l));
			}catch(Exception exc){
				System.err.println("Error parsing "+f.getName()+"\n> "+l);
				exc.printStackTrace();
			}
		}
		
		//Filter out duplicate ATOM records
		for(int i=1;i<ret.size();i++){
			if(ret.get(i-1) instanceof AtomRecord && ret.get(i) instanceof AtomRecord){
				AtomRecord atom1 = (AtomRecord)ret.get(i-1);
				AtomRecord atom2 = (AtomRecord)ret.get(i);
				if(atom1.aaType.equals(atom2.aaType) && atom1.atomName.equals(atom2.atomName)){
					ret.remove(i);
					i--;
				}
			}
		}
		return ret;
	}

	public String getSequence(){
		StringBuilder sb = new StringBuilder();
		for(AtomRecord a: getAtomRecords()){
			if(a.atomName.equalsIgnoreCase("CA")){
				sb.append(a.getSingleCharAAType());
			}
		}
		return sb.toString().toString();
	}
	
	public static interface PDBRecord{}

	public static class HelixRecord implements PDBRecord{
		public char initChainID, endChainID;
		public int initSeqNum, endSeqNum;

		HelixRecord(String pdbLine){
			initChainID = pdbLine.charAt(19);
			initSeqNum = Integer.parseInt(pdbLine.substring(21,25).trim());
			endChainID = pdbLine.charAt(31);
			endSeqNum = Integer.parseInt(pdbLine.substring(33,37).trim());
		}
	}

	public static class SheetRecord implements PDBRecord{
		public int strand;
		public String sheetID;
		public int numStrands, sense;
		//public String initResName, endResName, curResName, prevResName;
		public int initSeqNum, endSeqNum, curResSeq, prevResSeq;
		public char initChainID, endChainID, curChainID, prevChainID;

		SheetRecord(String pdbLine) throws Exception{
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

	public static class AtomRecord implements PDBRecord{
		public String atomName, aaType, element;
		public char chain;
		public int atomNumber, residueNumber;
		public Point coords;

		AtomRecord(String pdbLine){
			try{
			atomNumber = Integer.parseInt(pdbLine.substring(5,11).trim());
			atomName = pdbLine.substring(12,16).trim();
			aaType = pdbLine.substring(17, 20).trim();
			chain = pdbLine.charAt(21);
			residueNumber = Integer.parseInt(pdbLine.substring(22,26).trim());
			double x = Double.parseDouble(pdbLine.substring(28,38).trim());
			double y = Double.parseDouble(pdbLine.substring(38,46).trim());
			double z = Double.parseDouble(pdbLine.substring(46,54).trim());
			coords = new Point(x,y,z);
			element = pdbLine.substring(76,78).trim();
			//TODO: Read occupancy,tempFactor, element and charge
			}catch(StringIndexOutOfBoundsException exc){}//Ignore missing elements
		}
		public boolean isHydrogen() {
			if(element==null || element.isEmpty()){
				return atomName.startsWith("H");
			}else{
				return element.equalsIgnoreCase("H");
			}
		}
		/**
		 * Returns a string representation of this ATOM-record that follows the PDB-file format. 
		 * Following is from http://www.wwpdb.org/documentation/format32/sect9.html
		 * COLUMNS        DATA  TYPE    FIELD        DEFINITION
         * -------------------------------------------------------------------------------------
         *  1 -  6        Record name   "ATOM  "
         *  7 - 11        Integer       serial       Atom  serial number.
         * 13 - 16        Atom          name         Atom name.
         * 17             Character     altLoc       Alternate location indicator.
         * 18 - 20        Residue name  resName      Residue name.
         * 22             Character     chainID      Chain identifier.
         * 23 - 26        Integer       resSeq       Residue sequence number.
         * 27             AChar         iCode        Code for insertion of residues.
         * 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
         * 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
         * 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
         * 55 - 60        Real(6.2)     occupancy    Occupancy.
         * 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
         * 77 - 78        LString(2)    element      Element symbol, right-justified.
         * 79 - 80        LString(2)    charge       Charge  on the atom.
		 */
		public String toString(){
			return String.format("%6s%5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
					"ATOM",
					atomNumber, 
					atomName,
					' ',
					aaType,
					chain,
					residueNumber,
					' ',
					coords.getX(),
					coords.getY(),
					coords.getZ(),
					1.0,
					0.0,
					"",
					""
			);
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
	public static class TerRecord implements PDBRecord{
		public String aaType;
		public char chain;
		public int atomNumber, residueNumber;

		TerRecord(String pdbLine){
			try{
				aaType = pdbLine.substring(17, 20);
				chain = pdbLine.charAt(21);
				residueNumber = Integer.parseInt(pdbLine.substring(22,26).trim());
			}catch(Exception exc){}//Very few people follow the specifications for TER
		}
		public String toString(){
			return String.format("Ter[#:%d,res#:%d]",atomNumber, residueNumber);
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
	public static class HetatmRecord implements PDBRecord{
		public String atomType, haType;
		public char chain;
		public int atomNumber, residueNumber;
		public Point coords;

		HetatmRecord(String pdbLine){
			atomNumber = Integer.parseInt(pdbLine.substring(6,11).trim());
			atomType = pdbLine.substring(13,16).trim();
			haType = pdbLine.substring(17, 20);
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
			if(haType.equalsIgnoreCase("ALA")) return 'A';  
			if(haType.equalsIgnoreCase("ARG")) return 'R';  
			if(haType.equalsIgnoreCase("ASN")) return 'N';  
			if(haType.equalsIgnoreCase("ASP")) return 'D';  
			if(haType.equalsIgnoreCase("CYS")) return 'C';  
			if(haType.equalsIgnoreCase("GLU")) return 'E';  
			if(haType.equalsIgnoreCase("GLN")) return 'Q';  
			if(haType.equalsIgnoreCase("GLY")) return 'G';  
			if(haType.equalsIgnoreCase("HIS")) return 'H';  
			if(haType.equalsIgnoreCase("ILE")) return 'I';  
			if(haType.equalsIgnoreCase("LEU")) return 'L'; 
			if(haType.equalsIgnoreCase("LYS")) return 'K'; 
			if(haType.equalsIgnoreCase("MET")) return 'M'; 
			if(haType.equalsIgnoreCase("PHE")) return 'F'; 
			if(haType.equalsIgnoreCase("PRO")) return 'P'; 
			if(haType.equalsIgnoreCase("SER")) return 'S'; 
			if(haType.equalsIgnoreCase("THR")) return 'T'; 
			if(haType.equalsIgnoreCase("TRP")) return 'W'; 
			if(haType.equalsIgnoreCase("TYR")) return 'Y'; 
			if(haType.equalsIgnoreCase("VAL")) return 'V'; 
			return '?';
//			throw new Error("Unknown amino acid type: "+aaType);
		}
	}

	public static class ModelRecord implements PDBRecord{
		public int number;
		
		private ModelRecord(String pdbLine){
			try{
				number = Integer.parseInt(pdbLine.substring(6).trim());
			}catch(NumberFormatException exc){}//Some write REFINED after model number
		}
	}

	public static class EndModelRecord implements PDBRecord{
		private EndModelRecord(String pdbLine){	}
	}
	
	public static class RemarkRecord implements PDBRecord{
		public String remark;
		
		private RemarkRecord(String pdbLine){
			remark = pdbLine.substring(6);
		}
	}
	public static class ParentRecord implements PDBRecord{
		public String recordString;
		
		private ParentRecord(String pdbLine){
			recordString = pdbLine;
		}
		public String toString(){
			return recordString;
		}
	}

	public static class AnisouRecord implements PDBRecord{
		public String recordString;
		
		private AnisouRecord(String pdbLine){
			recordString = pdbLine;
		}
		public String toString(){
			return recordString;
		}
	}

	public static class OtherRecord implements PDBRecord{
		public String recordString;
		
		private OtherRecord(String pdbLine){
			recordString = pdbLine;
		}
		public String toString(){
			return recordString;
		}
	}
	
	
	public static class PDBModel{
		public int number;
		private List<PDBChain> chains = new ArrayList<PDBChain>();
		
		private PDBModel(int n){
			number = n;
		}

		public List<PDBChain> getChains(){ return chains; }
	}
	
	public static class PDBChain{
		char chainId;
		private List<AtomRecord> atomRecords = new ArrayList<AtomRecord>();

		PDBChain(char id){
			chainId = id;
		}
	}
}
