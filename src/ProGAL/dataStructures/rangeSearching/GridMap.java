package ProGAL.dataStructures.rangeSearching;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import ProGAL.geomNd.Point;

public class GridMap implements RangeSearchDataStructure{
	private final HashMap<Integer, List<Point>> map = new HashMap<Integer, List<Point>>();
	private final double cellSize;
	private final int dim;
	
	
	public GridMap(double cellSize, int dimensions){
		this.cellSize = cellSize;
		this.dim = dimensions;
	}
	

	public GridMap(double cellSize, List<Point> points){
		this.cellSize = cellSize;
		this.dim = points.get(0).getDimensions();
		for(Point p: points) addPoint(p);
	}
	
	public int[] cellFromPoint(Point p){
		int[] ret = new int[dim];
		for(int i=0;i<dim;i++) 
			ret[i] = (int)(p.get(i)/cellSize);
		return ret;
	}
	
	private List<Point> query(int[] lowCell, int[] highCell, int d, List<Point> ret){
		int origLow = lowCell[d];
		for(int i=lowCell[d];i<=highCell[d];i++){
			lowCell[d] = i;
			if(d<dim-1) query(lowCell,highCell,d+1,ret);
			else{
				List<Point> points = map.get(Arrays.hashCode(lowCell));
				if(points!= null) ret.addAll(points);
			}
		}
		lowCell[d] = origLow;
		return ret;
	}
	
	public List<Point> query(Point low, Point high){
		List<Point> ret = new LinkedList<Point>();
		int[] lowCell = cellFromPoint(low);
		int[] highCell = cellFromPoint(high);
		query(lowCell,highCell, 0, ret);
		
		//Some points might be in the queried cells but not in the window. Remove these. 
		for(Iterator<Point> it = ret.iterator(); it.hasNext(); ){
			Point p = it.next();
			for(int d=0;d<dim;d++){
				if(p.get(d)<low.get(d) || p.get(d)>high.get(d)){
					it.remove();
					break;
				}
			}
		}
		return ret;
	}
	
	public void addPoint(Point p){
		int cell = Arrays.hashCode(cellFromPoint(p));
		List<Point> pl = map.get(cell);
		if(pl!=null) pl.add(p);
		else {
			pl = new LinkedList<Point>();
			pl.add(p);
			map.put(cell, pl);
		}
	}

	
	
	public static void main(String[] args) {
		GridMap gm = new GridMap(1, 3);
		gm.addPoint(new Point(new double[]{0.7,0.7,0.7}));
		gm.addPoint(new Point(new double[]{0.3,0.3,0.3}));
		gm.addPoint(new Point(new double[]{-0.3,-0.3,-0.3}));
		gm.addPoint(new Point(new double[]{2.2,2.2,3.2}));
		List<Point> points = gm.query(new Point(new double[]{0.5,0.5,0.5}), new Point(new double[]{2.5,2.5,3.5}));
		System.out.println(points);

	}
}
