package ProGAL.proteins.structure.generators;

import java.util.LinkedList;
import java.util.List;

import ProGAL.proteins.structure.AminoAcid;
import ProGAL.proteins.structure.Atom;
import ProGAL.proteins.structure.CBond;


public class CABAminoAcidGenerator implements AtomGenerator{
	public CBond[] generateBonds(AminoAcid[] aminoAcids){
		List<CBond> cBonds = new LinkedList<CBond>();
		AminoAcid prev = null;
		for(AminoAcid aa: aminoAcids){
			
			Atom ca = aa.atom("CA");
			Atom cb = aa.atom("CB");
			cBonds.add(createBond(ca, cb));
			
			if(prev!=null){
				cBonds.add(  createBond(prev.atom("CA"), ca)  );
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
		Atom[] ret = {
				new Atom("CA"),
				new Atom("CB")
		};
		return ret;
	}
}
