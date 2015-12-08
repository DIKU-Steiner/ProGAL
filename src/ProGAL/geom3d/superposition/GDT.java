package ProGAL.geom3d.superposition;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import ProGAL.geom3d.Point;
import ProGAL.proteins.PDBFile;

public class GDT {
	public static final int REMOVE_ONE = 0;
	public static final int REMOVE_TEN_PERCENT = 1;

	public static void main(String[] args){
		PDBFile f1 = new PDBFile("/Users/ras/Documents/Datasets/CASP8Training/T0490/decoys/servers/schenk-torda-server_TS1");
		PDBFile f2 = new PDBFile("/Users/ras/Documents/Datasets/CASP8Training/T0490/T0490.pdb");
		List<Point> coords1 = f1.getCACoords();
		List<Point> coords2 = f2.getCACoords();
		coords1.remove(0);coords1.remove(0);coords1.remove(0);

		double cutoff = 8;
		System.out.printf("%.1f %.1f\n",cutoff, getGDT(coords1, coords2, cutoff));
	}

	public static double getGDT(List<Point> points1, List<Point> points2, double cutoff){
		return getGDT(points1, points2,cutoff, REMOVE_ONE);
	}
	public static double getGDT(List<Point> points1, List<Point> points2, double cutoff, int method){
		if(points1.size()!=points2.size()) throw new RuntimeException("Sizes of point-lists must match");
		
		int n = points1.size();
		
		//points1 remains static .. a clone of points2 is moved around
		List<Point> moving = new ArrayList<Point>(n);
		for(Point p: points2) moving.add(p.clone());

		switch(method){
		case REMOVE_ONE: return getGDT_REMOVE_ONE(points1, moving, cutoff);
//		case REMOVE_TEN_PERCENT: return getGDT_REMOVE_TEN_PERCENT(points1, moving);
		}
		throw new RuntimeException("Unknown method "+method);
	}
	
	private static double getGDT_REMOVE_ONE(List<Point> points1, List<Point> moving, double cutoff){
		List<Integer> included = new LinkedList<Integer>();
		for(int i=0;i<points1.size();i++) included.add(i);
		Transform t = RMSD.optimalSuperposition(moving, points1, included);
		t.transformIn(moving);
		double[] distances = getDistances(points1,moving);
		List<Integer>[] within = withinCutoff(distances, cutoff, included, null);
		
		while(!within[1].isEmpty()){
			System.out.println("Within:  "+within[0]);
			System.out.println("Without: "+within[1]);
			int max = maxDist(distances, included);
			included.remove(included.indexOf(max));
			
			RMSD.optimalSuperposition(moving, points1, included);
			t.transformIn(moving);
			distances = getDistances(points1,moving);
			within = withinCutoff(distances, cutoff, included, within);
		}
		
		return (within[0].size()*100.0)/points1.size();
	}
	
	private static int maxDist(double[] distances, List<Integer> included){
		int maxIdx = included.get(0);
		for(int i: included){
			if(distances[i]>distances[maxIdx]) maxIdx = i;
		}
		return maxIdx;
	}
	
	private static double[] getDistances(List<Point> points1, List<Point> points2){
		double[] ret = new double[points1.size()];
		Iterator<Point> i1 = points1.iterator();
		Iterator<Point> i2 = points2.iterator();
		int i=0;
		while(i1.hasNext()){
			ret[i++] = i1.next().distance(i2.next());
		}
		return ret;
	}
	
	/**
	 * Returns two lists specifying which indices of point-pairs are within cutoff distance (first list) 
	 * and which are not (second list). If the last argument is null a new pair of lists are allocated, 
	 * otherwise those in ret are cleared and reused. 
	 */
	@SuppressWarnings("unchecked")
	private static List<Integer>[] withinCutoff(double[] distances, double cutoff, List<Integer> included, List<Integer>[] ret){
		if(ret==null) {
			ret = new List[2];
			ret[0] = new LinkedList<Integer>();
			ret[1] = new LinkedList<Integer>();
		}
		else { ret[0].clear(); ret[1].clear(); }
		
		
		for(int i: included){
			if(distances[i]<=cutoff)
				ret[0].add(i);
			else
				ret[1].add(i);
		}
		
		return ret;
	}
}