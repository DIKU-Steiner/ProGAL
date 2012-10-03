package ProGAL.proteins.structure.generators;

import java.util.LinkedList;
import java.util.List;

import ProGAL.proteins.structure.AminoAcid;
import ProGAL.proteins.structure.Atom;
import ProGAL.proteins.structure.CBond;


public class AllAtomAminoAcidGenerator implements AtomGenerator{
	public CBond[] generateBonds(AminoAcid[] aminoAcids){
		List<CBond> cBonds = new LinkedList<CBond>();
		String[][] rules = {//Rules for all bonds except across amino acids and within pro
				{"N","CA"},
				{"C","O"},
				{"CA","C"},
				{"CA","CB"},
				{"CB","CG"},
				{"CB","CG1"},
				{"CB","CG2"},
				{"CB","OG"},
				{"CB","OG1"},
				{"CD","CE"},
				{"CD","NE"},
				{"CD","OE1"},
				{"CD","OE2"},
				{"CD1","CE1"},
				{"CD1","NE1"},
				{"CD2","CE2"},
				{"CD2","CE3"},
				{"CE","NZ"},
				{"CE1","CZ"},
				{"CE2","CZ"},
				{"CE2","CZ2"},
				{"CE3","CZ3"},
				{"CG","CD"},
				{"CG","CD1"},
				{"CG","CD2"},
				{"CG","ND2"},
				{"CG","OD1"},
				{"CG","OD2"},
				{"CG1","CD1"},
				{"CZ","NH1"},
				{"CZ","NH2"},
				{"CZ","OH"},
				{"CZ2","CH2"},
				{"CZ3","CH2"},
				{"NE","CZ"},
				{"NE1","CE2"},
				
				//Hydrogens
				{"N","H1"},
				{"N","H"},
				{"CA","HA"},
				{"","H"},
				{"","H"},
				{"","H"},
				{"","H"},
				{"","H"},
				{"","H"},
				
		};
		AminoAcid prev = null;
		for(AminoAcid aa: aminoAcids){
			for(String[] rule: rules) {
				try{
					Atom a1 = aa.atom(rule[0]);
					Atom a2 = aa.atom(rule[1]);
					cBonds.add(createBond(a1,a2));
				}catch(RuntimeException exc){}
			}
			if(aa.type()=='P')//Weird PRO
				cBonds.add(  createBond(aa.atom("CD"), aa.atom("N"))  );
			if(prev!=null){
				cBonds.add(  createBond(prev.atom("C"), aa.atom("N"))  );
			}
			prev = aa;
		}

		CBond[] allBonds = new CBond[cBonds.size()];
		int c=0;
		for(CBond b: cBonds) allBonds[c++] = b;
		return allBonds;
	}
	static CBond createBond(Atom a1, Atom a2){
		CBond bond = new CBond(a1,a2);
		if(a1.covalentBonds()==null) a1.setCovalentBonds(new CBond[0]);
		if(a2.covalentBonds()==null) a2.setCovalentBonds(new CBond[0]);
		Atom a = a1;
		CBond[] newBonds = new CBond[a.covalentBonds().length+1];
		for(int i=0;i<newBonds.length-1;i++) newBonds[i] = a.covalentBonds()[i];
		newBonds[newBonds.length-1] = bond;
		a.setCovalentBonds(newBonds);

		a = a2;
		newBonds = new CBond[a.covalentBonds().length+1];
		for(int i=0;i<newBonds.length-1;i++) newBonds[i] = a.covalentBonds()[i];
		newBonds[newBonds.length-1] = bond;
		a.setCovalentBonds(newBonds);

		return bond;
	}



	public Atom[] generateAtoms(char type){
		switch(type){
		case 'A': return generateAtomsA();
		case 'R': return generateAtomsR();	
		case 'N': return generateAtomsN();
		case 'D': return generateAtomsD();
		case 'C': return generateAtomsC();
		case 'E': return generateAtomsE();
		case 'Q': return generateAtomsQ();
		case 'G': return generateAtomsG();
		case 'H': return generateAtomsH();
		case 'I': return generateAtomsI();
		case 'L': return generateAtomsL();
		case 'K': return generateAtomsK();
		case 'M': return generateAtomsM();
		case 'F': return generateAtomsF();
		case 'P': return generateAtomsP();
		case 'S': return generateAtomsS();
		case 'T': return generateAtomsT();
		case 'W': return generateAtomsW();
		case 'Y': return generateAtomsY();
		case 'V': return generateAtomsV();
		}
		return null;
	}

	static Atom[] generateAtomsA(){ 
		Atom[] ret = {
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB1"),
				new Atom("HB2"),
				new Atom("HB3"),
		};
		return ret;
	}
	static Atom[] generateAtomsR(){ 
		Atom[] ret = {
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("CD"),
				new Atom("NE"),
				new Atom("CZ"),
				new Atom("NH1"),
				new Atom("NH2"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
				new Atom("HG2"),
				new Atom("HG3"),
				new Atom("HD2"),
				new Atom("HD3"),
				new Atom("HE"),
				new Atom("HH11"),
				new Atom("HH12"),
				new Atom("HH21"),
				new Atom("HH22"),
		};

		return ret;
	}
	static Atom[] generateAtomsN(){
		Atom[] ret = {
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("OD1"),
				new Atom("ND2"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
				new Atom("HD21"),
				new Atom("HD22"),
		};
		return ret;
	}
	static Atom[] generateAtomsD(){ 
		Atom[] ret = {
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("OD1"),
				new Atom("OD2"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
		};

		return ret;
	}
	static Atom[] generateAtomsC(){ 
		Atom[] ret = {
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("SG"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
		};

		return ret;
	}
	static Atom[] generateAtomsE(){  
		Atom[] ret = {
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("CD"),
				new Atom("OE1"),
				new Atom("OE2"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
				new Atom("HG2"),
				new Atom("HG3"),
		};

		return ret;
	}
	static Atom[] generateAtomsQ(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("CD"),
				new Atom("OE1"),
				new Atom("NE2"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
				new Atom("HG2"),
				new Atom("HG3"),
				new Atom("HE21"),
				new Atom("HE22"),
		};
	}
	static Atom[] generateAtomsG(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("H"),
				new Atom("HA2"),
				new Atom("HA3"),
		};
	}
	static Atom[] generateAtomsH(){
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("ND1"),
				new Atom("CD2"),
				new Atom("CE1"),
				new Atom("NE2"),
				new Atom("H"), 
				new Atom("HA"), 
				new Atom("HB2"), 
				new Atom("HB3"), 
				new Atom("HD1"),
				new Atom("HD2"),
				new Atom("HE1"),
		};
	}
	static Atom[] generateAtomsI(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG1"),
				new Atom("CG2"),
				new Atom("CD1"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB"),
				new Atom("HG12"),
				new Atom("HG13"),
				new Atom("HG21"),
				new Atom("HG22"),
				new Atom("HG23"),
				new Atom("HD11"),
				new Atom("HD12"),
				new Atom("HD13"),
		};
	}
	static Atom[] generateAtomsL(){
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("CD1"),
				new Atom("CD2"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
				new Atom("HG"),
				new Atom("HD11"),
				new Atom("HD12"),
				new Atom("HD13"),
				new Atom("HD21"),
				new Atom("HD22"),
				new Atom("HD23"),
		};
	}
	static Atom[] generateAtomsK(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("CD"),
				new Atom("CE"),
				new Atom("NZ"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
				new Atom("HG2"),
				new Atom("HG3"),
				new Atom("HD2"),
				new Atom("HD3"),
				new Atom("HE2"),
				new Atom("HE3"),
				new Atom("HZ1"),
				new Atom("HZ2"),
				new Atom("HZ3"),
		};
	}
	static Atom[] generateAtomsM(){
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("SD"),
				new Atom("CE"),
				new Atom("H"),  
				new Atom("HA"),  
				new Atom("HB2"),  
				new Atom("HB3"),  
				new Atom("HG2"),  
				new Atom("HG3"),  
				new Atom("HE1"),  
				new Atom("HE2"),  
				new Atom("HE3"),
		};
	}
	static Atom[] generateAtomsF(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("CD1"),
				new Atom("CD2"),
				new Atom("CE1"),
				new Atom("CE2"),
				new Atom("CZ"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
				new Atom("HD1"),
				new Atom("HD2"),
				new Atom("HE1"),
				new Atom("HE2"),
				new Atom("HZ"),
		};
	}
	static Atom[] generateAtomsP(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("CD"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
				new Atom("HG2"),
				new Atom("HG3"),
				new Atom("HD2"),
				new Atom("HD3"),
		};
	}
	static Atom[] generateAtomsS(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("OG"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB2"),
				new Atom("HB3"),
				new Atom("HG"),
		};
	}
	static Atom[] generateAtomsT(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("OG1"),
				new Atom("CG2"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB"),
				new Atom("HG1"),
				new Atom("HG21"),
				new Atom("HG22"),
				new Atom("HG23"),
		};
	}
	static Atom[] generateAtomsW(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("CD1"),
				new Atom("CD2"),
				new Atom("NE1"),
				new Atom("CE2"),
				new Atom("CE3"),
				new Atom("CZ2"),
				new Atom("CZ3"),
				new Atom("CH2"),
				new Atom("H"),  
				new Atom("HA"),  
				new Atom("HB2"),  
				new Atom("HB3"),  
				new Atom("HD1"),  
				new Atom("HE1"),  
				new Atom("HE3"),  
				new Atom("HZ2"),  
				new Atom("HZ3"),  
				new Atom("HH2"),
		};
	}
	static Atom[] generateAtomsY(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG"),
				new Atom("CD1"),
				new Atom("CD2"),
				new Atom("CE1"),
				new Atom("CE2"),
				new Atom("CZ"),
				new Atom("OH"),
				new Atom("H"),  
				new Atom("HA"),  
				new Atom("HB2"),  
				new Atom("HB3"),  
				new Atom("HD1"),  
				new Atom("HD2"),  
				new Atom("HE1"),  
				new Atom("HE2"),  
				new Atom("HH"), 
		};
	}
	static Atom[] generateAtomsV(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O"),
				new Atom("CB"),
				new Atom("CG1"),
				new Atom("CG2"),
				new Atom("H"),
				new Atom("HA"),
				new Atom("HB"),
				new Atom("HG11"),
				new Atom("HG12"),
				new Atom("HG13"),
				new Atom("HG21"),
				new Atom("HG22"),
				new Atom("HG23"),
		};
	}
}
