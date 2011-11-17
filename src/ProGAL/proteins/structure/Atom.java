package ProGAL.proteins.structure;

import ProGAL.geom3d.Point;

public class Atom extends Point{
	private static final long serialVersionUID = -6262882698432071090L;
	
	private final double radius;
	private final char element;
	private final String name;
	final int id;
//	private final Point center;
	private CBond[] covalentBonds;
	
	protected AminoAcid aminoAcid;
	protected int index; 

	public Atom(String name){
		super(0,0,0);
		this.name = name;
		this.element = name.charAt(0);
		this.id = name.hashCode();
		switch(element){
		case 'C': radius = 1.7;break;
		case 'N': radius = 1.55;break;
		case 'O': radius = 1.52;break;
		case 'S': radius = 1.8;break;
		case 'H': radius = 1.2;break;
		default: radius = 0;
		}
//		this.center = new Point(0,0,0);
	}

	public Atom(String name, char element, double radius){
		super(0,0,0);
		this.name = name;
		this.element = element;
		this.id = name.hashCode();
		this.radius = radius;
//		center = new Point(0,0,0);
	}
	
	public void setCovalentBonds(CBond[] bonds){
		this.covalentBonds = bonds;
	}
	
	public double radius(){ return radius; }
	public char element() { return element; }
//	public Point center() { return center; }
	public String name()  { return name; }
	public AminoAcid aminoAcid(){ return aminoAcid; }
	public int index()	  { return index; }
	public CBond[] covalentBonds(){ return covalentBonds; } 
	
	public boolean isBB(){
		return 
		name.equalsIgnoreCase("CA") ||
		name.equalsIgnoreCase("C") ||
		name.equalsIgnoreCase("N") ||
		name.equalsIgnoreCase("O");
	}
	
	public String toString(){
		return String.format("Atom[%s_%d,%s,%.2f,%s]",
				aminoAcid.typeThreeLetter(),
				aminoAcid.index, 
				name,
				radius,
				super.toString());
	}
	
	public Atom clone(){
		Atom ret = new Atom(name);
		ret.set(this);
		return ret;
	}
}
