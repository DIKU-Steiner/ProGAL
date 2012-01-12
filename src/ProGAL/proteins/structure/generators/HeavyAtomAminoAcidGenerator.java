package ProGAL.proteins.structure.generators;

import java.util.LinkedList;
import java.util.List;

import ProGAL.proteins.structure.AminoAcid;
import ProGAL.proteins.structure.Atom;
import ProGAL.proteins.structure.CBond;


public class HeavyAtomAminoAcidGenerator implements AtomGenerator{
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
				{"NE1","CE2"}
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
				new Atom("CB")
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
				new Atom("ND2")
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
				new Atom("OD2")
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
				new Atom("SG")
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
				new Atom("OE2")
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
				new Atom("NE2")
		};
	}
	static Atom[] generateAtomsG(){ 
		return new Atom[]{
				new Atom("N"),
				new Atom("CA"),
				new Atom("C"),
				new Atom("O")
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
				new Atom("NE2")
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
		};
	}
}
