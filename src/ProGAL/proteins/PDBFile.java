package ProGAL.proteins;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.complex.alphaComplex.AlphaComplex;
import ProGAL.geom3d.complex.CTriangle;
import ProGAL.geom3d.superposition.RMSD;
import ProGAL.geom3d.superposition.Transform;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.io.IOToolbox;
import ProGAL.io.WebIOToolbox;
import ProGAL.proteins.structure.AminoAcidChain;
import ProGAL.proteins.structure.generators.HeavyAtomAminoAcidGenerator;

/**
 * A class for reading and representing a PDB-file. All lines are read but not all are parsed.
 * All the records can be retrieved as they are shown in the PDB-file. The ATOM-records, however, 
 * can also be accessed according to model and chain id. If a model-number and chain id are not 
 * specified the method will assume that the first chain of the first model is meant.  
 *  
 * In some files there are multiple suggestions for atom-placements. This is resolved by only 
 * reading the first suggested atom-placement. The other atom-placements are added to the 
 * alternativeCoords LinkedList in the first atom-record. 
 * 
 * The four-letter PDB id and chain id are guessed based on the filename.
 * 
 * TODO: Represent ligands, heterogens and secondary structure.
 * TODO: Check the integrity of different models, e.g. that the sequences are the same.
 * @author R. Fonseca
 */
public class PDBFile extends File{
	private static final long serialVersionUID = 1791445213018199901L;
	
	private boolean includeHydrogens = true;
	private boolean includeHetAtms = false;
	private int standardModel = 0;
	private int standardChain = 0;
	protected final List<PDBRecord> records;
	protected final List<PDBModel> models = new ArrayList<PDBModel>();
	
	/** The name of this PDB-file. Typically the PDB-id. */
	public String name;

	/** Construct a PDB-file from the specified file-path */
	public PDBFile(String path, boolean hydro){
		this(new File(path), hydro);
	}
	/** Construct a PDB-file from the specified file-path */
	public PDBFile(String path){
		this(new File(path), false);
	}
	
	/** Construct a PDB-file from the specified file */
	public PDBFile(File f, boolean hydro){
		super(f.getAbsolutePath());
		if(f.getName().contains("."))	name = f.getName().substring(0,f.getName().indexOf('.'));
		else							name = f.getName();
		records = readPDBFile(f,false);
		buildPDBStructure();
		if(name.contains("_")) {
			setStandardChain(name.charAt(name.indexOf('_')+1));
		}
		this.includeHydrogens = hydro;
	}

	private void buildPDBStructure(){
		
		//Build models and chains and fill them with ATOM-records
		for(int r=0;r<records.size();r++){
			if(records.get(r) instanceof AtomRecord) 		r = buildModel(r);
			else if(records.get(r) instanceof ModelRecord) 	r = buildModel(r);
			
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
//		if(model==null) throw new RuntimeException("Expected a model at "+name+":"+(r+1));
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
			if(record instanceof TerRecord) continue;
			else if(record instanceof HetatmRecord){
				HetatmRecord aRecord = (HetatmRecord)record;
				if(chain==null) chain = new PDBChain(aRecord.chain);
				if(chain.chainId!=aRecord.chain) break;
				
				//Ignore alternative coordinates
				AtomRecord prevARecord = locateOccupancyPartner(aRecord, chain);
				if(prevARecord!=null){
					if(prevARecord.occupancy>=aRecord.occupancy)
						prevARecord.alternativeCoords.add(aRecord.coords);
					else{
						prevARecord.alternativeCoords.add(prevARecord.coords);
						prevARecord.coords = aRecord.coords;
						prevARecord.occupancy = aRecord.occupancy;
					}
				}else{
					chain.atomRecords.add(aRecord);
				}
			}else if(record instanceof AtomRecord){
				AtomRecord aRecord = (AtomRecord)record;
				if(chain==null) chain = new PDBChain(aRecord.chain);
				if(chain.chainId!=aRecord.chain) break;
				
				//Ignore alternative coordinates
				AtomRecord prevARecord = locateOccupancyPartner(aRecord, chain);
				if(prevARecord!=null){
					if(prevARecord.occupancy>=aRecord.occupancy)
						prevARecord.alternativeCoords.add(aRecord.coords);
					else{
						prevARecord.alternativeCoords.add(prevARecord.coords);
						prevARecord.coords = aRecord.coords;
						prevARecord.occupancy = aRecord.occupancy;
					}
				}else{
					chain.atomRecords.add(aRecord);
				}
			}else if(record instanceof AnisouRecord){
				//TODO: Do something about ANISOU's
			}else{ break; }//throw new RuntimeException("Unexpected record at "+name+":"+(r+1)+"\n> "+record); }
		}
		if(chain==null) throw new RuntimeException("Expected a chain at "+name+":"+(r+1)+". Is model empty?");
		m.chains.add(chain);
		return r;
	}
	
	/** 
	 * Searches through all AtomRecords in chain for a multiple occupancy to r.
	 * Returns null if none is found 
	 */
	private AtomRecord locateOccupancyPartner(AtomRecord r, PDBChain chain){
		for(AtomRecord ar: chain.atomRecords){
			if(ar==r) break;
			if(ar.atomType.equalsIgnoreCase(r.atomType) && ar.residueNumber==r.residueNumber) return ar;
		}
		return null;
	}
	
	/** Read a PDB-file from www.pdb.org.*/
	public static PDBFile downloadPDBFile(String pdbId, boolean hydro){
		return new PDBFile(WebIOToolbox.downloadFile("http://www.pdb.org/pdb/files/"+pdbId.toUpperCase()+".pdb"), hydro);
	}

	public int  getStandardModel()      { return this.standardModel; }
	public void setStandardModel(int m) { this.standardModel = m; }
	public int  getStandardChain()      { return this.standardChain; }
	public void setStandardChain(int c) { this.standardChain = c; }
	public void setIncludeHydrogens(boolean b){ this.includeHydrogens=b; }
	public void setIncludeHetatms(boolean b){ this.includeHetAtms=b; }
	public void setStandardChain(char c){
		PDBModel m = models.get(standardModel);
		for(PDBChain chain: m.chains){
			if(Character.toUpperCase(chain.chainId)==Character.toUpperCase(c)){
				setStandardChain(m.chains.indexOf(chain));
				break;
			}
		}
	}

	public List<PDBModel> getModels(){ return models; }

	/** Returns the ATOM-records. Only the records in the first chain of the first model are returned. */
	public List<AtomRecord> getAtomRecords(){
		return getAtomRecords(standardModel, standardChain);
	}

	/** Returns the ATOM-records. Only the records with the types specified in the comma-separated string
	 * are returned.
	 * @param atomTypes A comma separated list of atom types. Could, for example, be "CA,C,N,O" to specify
	 * backbone atoms only.*/
	public List<AtomRecord> getAtomRecords(String atomTypes){
		String[] atomTypeArr = atomTypes.split(",");
		List<AtomRecord> records = getAtomRecords();
		List<AtomRecord> ret = new ArrayList<AtomRecord>(records.size());
		for(AtomRecord rec: records) {
			for(String atomType: atomTypeArr)
				if(rec.atomType.equalsIgnoreCase(atomType)){
					ret.add(rec);
				}
		}
		return ret;
	}
	
	/** Returns the ATOM-records of CA-atoms. Only the records in the first chain of the first model are returned. */
	public List<AtomRecord> getCARecords(){
		return getCARecords(standardModel, standardChain);
	}
	
	/** Returns the CA-coordinates of atoms in the first chain in the first model. */
	public List<Point> getCACoords(){
		return getCACoords(standardModel, standardChain);
	}

	/** Returns all atom-coordinates in the first chain in the first model. */
	public List<Point> getAtomCoords(){
		return getAtomCoords(standardModel, standardChain);
	}

	/** Returns the coordinates of the specified atom types.
	 * @param atomTypes A comma separated list of atom types. Could, for example, be "CA,C,N,O" to specify
	 * backbone atoms only.*/
	public List<Point> getAtomCoords(String atomTypes){
		List<AtomRecord> records = getAtomRecords(atomTypes);
		List<Point> ret = new ArrayList<Point>(records.size());
		for(AtomRecord rec: records) 
			ret.add(rec.coords);
		return ret;
	}

	public void writeAtomCoords(String fileName, List<Point> atoms) throws FileNotFoundException, UnsupportedEncodingException {
		PrintWriter writer = new PrintWriter(fileName, "UTF-8");
		for (Point p: atoms) writer.println(p.x() + " " + p.y() + " " + p.z());
		writer.close();
	}
	
	/** Returns the ATOM-records of the specified model and chain. */
	public List<AtomRecord> getAtomRecords(int modelNum, int chainNum){
		List<AtomRecord> ret = new ArrayList<AtomRecord>();
		for(AtomRecord ar: models.get(modelNum).chains.get(chainNum).atomRecords){
//			if(!ar.isOnBackbone()) continue;       // to be removed if all atoms are to be included
			if(!includeHydrogens && ar.isHydrogen()) continue;
			if(!includeHetAtms && ar instanceof HetatmRecord) continue;
			
			ret.add(ar);
		}
		return ret;
	}

	public List<HetatmRecord> getHetatmRecords(){
		return getHetatmRecords(standardModel, standardChain);
	}
	
	/** Returns the HETATM-records of the specified model and chain. */
	public List<HetatmRecord> getHetatmRecords(int modelNum, int chainNum){
		List<HetatmRecord> ret = new ArrayList<HetatmRecord>();
//		for(AtomRecord ar: models.get(modelNum).chains.get(chainNum).atomRecords){
//			if(!includeHydrogens && ar.isHydrogen()) continue;
//			if(!(ar instanceof HetatmRecord)) continue;
		for(PDBRecord ar: records){
			if (ar instanceof HetatmRecord){
				if(!includeHydrogens && ((HetatmRecord)ar).isHydrogen()) continue;
				ret.add((HetatmRecord)ar);
			}
		}
		
		return ret;
	}

	/** Returns all records of the specified model and chain. */
	public List<PDBRecord> getRecords(){
		return new ArrayList<PDBRecord>(records);
	}
	/** Returns the CA ATOM-records of the specified model and chain. */
	public List<AtomRecord> getCARecords(int modelNum, int chainNum){
		List<AtomRecord> ret = new ArrayList<AtomRecord>();
		for(AtomRecord ar: models.get(modelNum).chains.get(chainNum).atomRecords){
			if(ar.atomType.equalsIgnoreCase("CA")) 
				ret.add(ar);
		}
		
		return ret;
	}
	
	/** Returns the CA-coordinates of the specified model and chain. */
	public List<Point> getCACoords(int modelNum, int chainNum){
		List<Point> caCoords = new ArrayList<Point>();
		for(AtomRecord a: getAtomRecords(modelNum,chainNum)) 
			if(a.atomType.equalsIgnoreCase("CA"))
				caCoords.add(a.coords);
		return caCoords;
	}

	/** Returns all atom-coordinates of the specified model and chain. */
	public List<Point> getAtomCoords(int modelNum, int chainNum){
		List<Point> coords = new ArrayList<Point>();
		for(AtomRecord a: getAtomRecords(modelNum,chainNum)) {
			if(!includeHydrogens && a.isHydrogen()) continue;
			if(!includeHetAtms && a instanceof HetatmRecord) continue;
			coords.add(a.coords);
		}
//		for(HetatmRecord a: getHetatmRecords(modelNum, chainNum)){
//			if(!includeHydrogens && a.isHydrogen()) continue;
//			coords.add(a.coords);
//		}
		return coords;
	}
	
	public AminoAcidChain getChain(int modelNum, int chainNum){
		AminoAcidChain chain = new AminoAcidChain(this.getSequence(), new HeavyAtomAminoAcidGenerator());
		int c=0;
		int prevRes = getAtomRecords(modelNum,chainNum).get(0).residueNumber;
		for(AtomRecord ar: getAtomRecords(modelNum, chainNum)){
			if(prevRes!=ar.residueNumber) c++;
			prevRes = ar.residueNumber;
			try{
				chain.atom(c, ar.atomType).set(ar.coords);
			}catch(RuntimeException exc){}
		}
//		Point O = new Point(0,0,0);
//		AminoAcid[] aas = chain.aminoAcids();
//		for(int aa=0;aa<aas.length;aa++){
//			for(Atom a: aas[aa].atoms()){
//				if(a.equals(O)) 
//					System.err.printf("Warning: %s%d_%s not set\n",aas[aa].typeThreeLetter(),aa,a.name());
//			}
//		}
		
		return chain;
	}
	
	/** Returns all REMARK records. */
	public List<RemarkRecord> getRemarkRecords(){
		List<RemarkRecord> ret = new ArrayList<RemarkRecord>();
		for(PDBRecord r: records) {
			if(r instanceof RemarkRecord)
				ret.add((RemarkRecord)r);
		}
		return ret;
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
//		for(int i=1;i<ret.size();i++){
//			if(ret.get(i-1) instanceof AtomRecord && ret.get(i) instanceof AtomRecord){
//				AtomRecord atom1 = (AtomRecord)ret.get(i-1);
//				AtomRecord atom2 = (AtomRecord)ret.get(i);
//				if(atom1.aaType.equals(atom2.aaType) && atom1.atomName.equals(atom2.atomName)){
//					ret.remove(i);
//					i--;
//				}
//			}
//		}
		return ret;
	}

	public String getSequence(){
		StringBuilder sb = new StringBuilder();
		for(AtomRecord a: getAtomRecords()){
			if(a.atomType.equalsIgnoreCase("CA")){
				sb.append(a.getSingleCharAAType());
			}
		}
		return sb.toString().toString();
	}
	
	public double getResolution(){
		for(PDBRecord record: records){
			if(record instanceof RemarkRecord){
				RemarkRecord rr = (RemarkRecord)record;
				if(rr.remark.contains("2 RESOLUTION.")){
					int start = rr.remark.indexOf("ANGSTROMS")-5;
					if(start<0) return Double.NaN;
					return Double.parseDouble(rr.remark.substring(start, start+4));
				}
			}
		}
		throw new RuntimeException(this.name+" - Resolution not specified");
	}
	

	public double superposeOnto(PDBFile f){
		List<Point> thisAtoms = new ArrayList<Point>();
		for(AtomRecord ar: getAtomRecords()) thisAtoms.add(ar.coords);
		List<Point> fAtoms = new ArrayList<Point>();
		for(AtomRecord ar: f.getAtomRecords()) fAtoms.add(ar.coords);

		while(fAtoms.size()>thisAtoms.size()) fAtoms.remove(fAtoms.size()-1);
		while(fAtoms.size()<thisAtoms.size()) thisAtoms.remove(thisAtoms.size()-1);
		
		Transform t = RMSD.optimalSuperposition(thisAtoms, fAtoms);
		t.transformIn(thisAtoms);
		return RMSD.getRMSD(thisAtoms, fAtoms);
	}
	
	
	
	
	
	
	
	public static interface PDBRecord{
	}

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
//				sense = Integer.parseInt(pdbLine.substring(38,40).trim());
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
		/** The type of atom. Backbone types are N, CA, C and O. Side-chain types can be CB, CG1 etc. */
		public String atomType;
		/** The type of amino acid in a three-letter code, for example: TYR, ALA or ARG */
		public String aaType;
		/** The atom element. Typically either of C, N, O, S or H. */
		public String element;
		/** A character identifying the chain */
		public char chain;
		public int atomNumber, residueNumber;
		public Point coords;
		public double occupancy;
		public LinkedList<Point> alternativeCoords = new LinkedList<Point>();

		AtomRecord(String pdbLine){
			try{
			atomNumber = Integer.parseInt(pdbLine.substring(6,11).trim());
			atomType = pdbLine.substring(12,16).trim();
			aaType = pdbLine.substring(17, 20).trim();
			chain = pdbLine.charAt(21);
			residueNumber = Integer.parseInt(pdbLine.substring(22,26).trim());
			double x = Double.parseDouble(pdbLine.substring(28,38).trim());
			double y = Double.parseDouble(pdbLine.substring(38,46).trim());
			double z = Double.parseDouble(pdbLine.substring(46,54).trim());
			coords = new Point(x,y,z);
			occupancy = Double.parseDouble(pdbLine.substring(56,60).trim());
			element = pdbLine.substring(76,78).trim();
			if(element.isEmpty()) element = pdbLine.charAt(13)+"";
			//TODO: Read occupancy,tempFactor, and charge
			}catch(StringIndexOutOfBoundsException exc){
				element = pdbLine.charAt(13)+""; 
			}
		}
		public boolean isOnBackbone() {
			return atomType.equals("C") || atomType.equals("CA") || atomType.equals("N");
		}
		public boolean isHydrogen() {
			if(element==null || element.isEmpty()){
				return atomType.startsWith("H");
			}else{
				return element.equalsIgnoreCase("H");
			}
		}
		
		/** Returns true if the atom can act as a hydrogen donor in a hydrogen bond */
		public boolean isHydrogenDonor() {
			if (element.equalsIgnoreCase("N")) {
				if (atomType.equalsIgnoreCase("N")) return true;
				if (aaType.equalsIgnoreCase("HIS")) {
					if (atomType.equalsIgnoreCase("NE2")) return true;
					if (atomType.equalsIgnoreCase("ND1")) return true;
					return false;
				}
				if (aaType.equalsIgnoreCase("LYS") && atomType.equalsIgnoreCase("NZ")) return true;
				if (aaType.equalsIgnoreCase("ASN") && atomType.equalsIgnoreCase("ND2")) return true;
				if (aaType.equalsIgnoreCase("GLN") && atomType.equalsIgnoreCase("NE2")) return true;
				if (aaType.equalsIgnoreCase("ARG")) {
					if (atomType.equalsIgnoreCase("NE")) return true;
					if (atomType.equalsIgnoreCase("NH1")) return true;
					if (atomType.equalsIgnoreCase("NH2")) return true;
					return false;
				}
				if (aaType.equalsIgnoreCase("TRP") && atomType.equalsIgnoreCase("NE1")) return true;
				return false;
			}
			if (element.equalsIgnoreCase("O")) {
				if (aaType.equalsIgnoreCase("SER") && atomType.equalsIgnoreCase("OG")) return true;
				if (aaType.equalsIgnoreCase("THR") && atomType.equalsIgnoreCase("OG1")) return true;
				if (aaType.equalsIgnoreCase("TYR") && atomType.equalsIgnoreCase("OH")) return true;
				if (aaType.equalsIgnoreCase("ASN") && atomType.equalsIgnoreCase("OD1")) return true;
				return false;
			}
			return false;
		}
		
		/** Returns true if the atom can act as a hydrogen acceptor in a hydrogen bond */
		public boolean isHydrogenAcceptor() {
			if (element.equalsIgnoreCase("O")) {
				if (atomType.equalsIgnoreCase("O")) return true;
				if (aaType.equalsIgnoreCase("ASP")) {
					if (atomType.equalsIgnoreCase("OD1")) return true;
					if (atomType.equalsIgnoreCase("OD2")) return true;
					return false;
				}
				if (aaType.equalsIgnoreCase("GLU")) {
					if (atomType.equalsIgnoreCase("OE1")) return true;
					if (atomType.equalsIgnoreCase("OE2")) return true;
					return false;
				}
				if (aaType.equalsIgnoreCase("ASN") && atomType.equalsIgnoreCase("OD1")) return true;
				if (aaType.equalsIgnoreCase("GLN") && atomType.equalsIgnoreCase("OE1")) return true;
				if (aaType.equalsIgnoreCase("SER") && atomType.equalsIgnoreCase("OG")) return true;
				if (aaType.equalsIgnoreCase("THR") && atomType.equalsIgnoreCase("OG1")) return true;
				if (aaType.equalsIgnoreCase("TYR") && atomType.equalsIgnoreCase("OH")) return true;
				return false;
			}
			if (element.equalsIgnoreCase("N")) {
				if (aaType.equalsIgnoreCase("HIS") && atomType.equalsIgnoreCase("ND1")) return true;
				if (aaType.equalsIgnoreCase("GLN") && atomType.equalsIgnoreCase("NE2")) return true;
				if (aaType.equalsIgnoreCase("ASN") && atomType.equalsIgnoreCase("ND2")) return true;
				return false;
			}
			return false;
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
			return String.format("%6s%5d %-4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s",
					(this instanceof HetatmRecord)?"HETATM":"ATOM  ",
					atomNumber, 
					atomType,
					' ',
					aaType,
					chain,
					residueNumber,
					' ',
					coords.x(),
					coords.y(),
					coords.z(),
					occupancy,
					0.0,
					element==null?"":element,
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
	public static class HetatmRecord extends AtomRecord{
		HetatmRecord(String pdbLine){
			super(pdbLine);
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
		
		private RemarkRecord(String pdbLine){	remark = pdbLine.substring(6);		}
		public String toString(){	return "REMARK"+remark;	}
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
		public List<AtomRecord> getAtoms(){
			return new ArrayList<AtomRecord>(atomRecords);
		}
	}

	public static void main(String[] args) throws FileNotFoundException, UnsupportedEncodingException {
		
		PDBFile f = new PDBFile("/home/daisy/Downloads/2oed_cs_244_samples/sample_000267200000_2_91.576430.pdb", true);//folded
//		PDBFile f = new PDBFile("/home/daisy/Downloads/2oed_cs_244_samples/sample_000000200000_2_284.178215.pdb", true);
		
//		PDBFile f = new PDBFile("/home/daisy/Downloads/1X0O.pdb", true);
		J3DScene scene = J3DScene.createJ3DSceneInFrame();
		List<Point> points = f.getAtomCoords();
		AlphaComplex ac = new AlphaComplex(points, 2.8);
		System.out.println("Point set size : "+points.size());
		/*List<Point> pointList = new ArrayList<Point>(points.size()/2);
		Vector translate = new Vector(-points.get(0).getCoord(0), -points.get(0).getCoord(1), -points.get(0).getCoord(2)); 
		for (int i=0 ; i<343 ; i++) {
			pointList.add(points.get(i).add(translate));
		}*/
/*		List<AtomRecord> AR = f.getAtomRecords();
//		points.addAll(f.getAtomCoords(0, 1));
		int maxAA = AR.get(AR.size()-1).residueNumber;
		int indexAA = (int)Math.ceil(maxAA*0.5);
		System.out.println("IndexAA = "+indexAA);
		
		f.writeAtomCoords("/Users/pawel/Downloads/1X0O_Coordinates.tex", points);
		
		//int j = 0;
		//for (j=0 ; j<AR.size() ; j++) {
		//	if (AR.get(j).residueNumber==indexAA ) break;
		//}
		//Line l = new Line(points.get(80), points.get(83));
		System.out.println(points.size());
//		AlphaComplex ac = new AlphaComplex(points, 2.0);
		
		l.toScene(scene, 0.1, Color.CYAN);
		List<Integer> rotIndices = new ArrayList<Integer>(Arrays.asList(84, 85, 86, 89, 90, 91, 92, 93, 94, 95, 96, 97));
		j = 0;
		for (int i=0 ; i<points.size() ; i++) {
			Sphere s = new Sphere(points.get(i), 0.9);
			Color clr= Color.black;
			//if (AR.get(i).residueNumber>= indexAA) {
			if (rotIndices.contains(i)) {
				clr = Color.RED;
			}
			scene.addShape(s, clr);
//			scene.addShape(new TextShape(Integer.toString(i)+" - "+AR.get(i-1).aaType, pointList.get(i-1), 0.3));
		}
		for (CEdge edg : ac.getEdges()) {
			System.out.println("new Edge");
			LSS lss = new LSS(edg.getA(), edg.getB(), 0.02);
			scene.addShape(lss, new Color(0, 0, 150, 100));
			try { Thread.sleep(30); } catch (InterruptedException e) {}
		}*/
		
/*		for(CTetrahedron tetr : ac.getTetrahedra()) {
			scene.addShape(tetr, new Color(0, 255, 255, 100));
		}*/
		for (CTriangle tri : ac.getAlphaShape(2.8)) {
				Triangle t = new Triangle(tri.getP1(), tri.getP2(), tri.getP3());
				t.toScene(scene, new Color(50,205,50, 200));
		}
		
		
	}

}
