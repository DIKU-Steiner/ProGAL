package ProGAL.geom3d.tessellation.BowyerWatson;

import java.awt.Color;
import java.awt.event.MouseEvent;
import java.math.BigDecimal;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import ProGAL.geom3d.Point;
import ProGAL.geom3d.PointWeighted;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.viewer.ClickListener;
import ProGAL.geom3d.viewer.J3DScene;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.math.Matrix;
import ProGAL.math.Randomization;

public class RegularTessellation {
	private final Tetr bigTetr;
	private final LinkedList<Tetr> tetras = new LinkedList<Tetr>();
	private final LinkedList<PointWeighted> points = new LinkedList<PointWeighted>();
	private int lastSize = 1;
	private final Queue<Tetr> lastQueue = new LinkedList<Tetr>();

	public static void main(String[] args){
		Tetr t = new Tetr(
				new PointWeighted(9.80,9.19,-4.53, 0),
				new PointWeighted(8.73,9.54,-0.53, 0),
				new PointWeighted(9.44,9.78, 5.35, 0),
				new PointWeighted(9.97,8.47, 2.45, 0));

		Randomization.seed(0);
		int insertedPoints = 1000000;
		LinkedList<PointWeighted> points= new LinkedList<PointWeighted>();
		for(int i=0;i<insertedPoints;i++){
			PointWeighted pw = new PointWeighted( 
					Randomization.randBetween(-10.0, 10.0),
					Randomization.randBetween(-10.0, 10.0),
					Randomization.randBetween(-10.0, 10.0) ,0);
			points.add(pw);
		}
		long start = System.nanoTime();
		//		DelaunayComplex rc = new DelaunayComplex(points);
		RegularTessellation rc = new RegularTessellation(points);
		long end = System.nanoTime();
		System.out.printf("%d points took %.3fms\n",insertedPoints,(end-start)/1000000.0);


	}

	public RegularTessellation(){
		bigTetr = createBigTetr();
	}


	public RegularTessellation(List<PointWeighted> points){
		this();
		int c=0;
		for(PointWeighted pw: points)	{
			insertPoint(pw);
			if(++c%10000==0){
				System.out.println(c);
			}
		}
	}

	private void insertPoint(PointWeighted p){
		Tetr tet = walk(p);

		LinkedList<Tetr> star = new LinkedList<Tetr>();//TODO: Replace with HashSet for faster contains method
		collectStar(p, tet, star);
		LinkedList<Tetr> newStar = new LinkedList<Tetr>();

		//Create new tets
		for(Tetr t: star){
			for(int i=0;i<4;i++){
				if(t.neighbors[i]==null || !star.contains(t.neighbors[i])) {
					Tetr newTet = new Tetr(p, t.corners[(i+1)&3], t.corners[(i+2)&3], t.corners[(i+3)&3]);
					newTet.neighbors[0] = t.neighbors[i];
					if(t.neighbors[i]!=null) {
						oneConnect(t.neighbors[i],newTet);
					}
					newStar.add(newTet);
				}
			}
		}
		//Connect tets
		for(Tetr t1: newStar){
			for(Tetr t2: newStar){
				if(t1==t2) break;
				connect(t1,t2);
			}
		}

		points.add(p);

		lastQueue.add(newStar.get(0));
		lastSize = Math.max(1,(int)Math.pow(points.size(), 0.25));
		while(lastQueue.size()>lastSize) lastQueue.poll();
	}

	private final static int[] shared1 = new int[3];
	private final static int[] shared2 = new int[3];

	/** Assumes that t1.corner[0]==t2.corner[0] */
	private static void connect(Tetr t1, Tetr t2){

		int s = 1;
		shared1[0] = 0;
		shared2[0] = 0;

		for(int i=1;i<4;i++){
			int partner = -1;
			for(int c=1;c<4;c++) if(t2.corners[c]==t1.corners[i]) {partner=c; break;}
			if(partner>0) {
				shared1[s]=i;
				shared2[s]=partner;
				s++;
				if(s==3 || (i>s+1)) break;
			}
		}
		if(s==3) {
			t1.neighbors[excluded(shared1)] = t2;
			t2.neighbors[excluded(shared2)] = t1;
		}
	}

	/** arr is assumed to hold three distinct values between 0 and 3. This method returns 
	 * the value not present in the array. */
	private static int excluded(int[] arr){
		search: for(int i=1;i<4;i++){
			for(int v: arr) if(v==i) continue search;
			return i;
		}
	return -1;
	}

	/** Only connects t1 to t2 and assumes that t2.corners[0] is not on the shared face. */
	private static void oneConnect(Tetr t1, Tetr t2){
		for(int i=0;i<4;i++){
			int idx = t2.cornerIdx(t1.corners[i]);
			if(idx<0) {
				t1.neighbors[i] = t2;
				return;
			}
		}
	}


	private void collectStar(PointWeighted p, Tetr t, Collection<Tetr> star){
		//		System.out.printf("collectStar(%s, %s ...)\n",p,t);
		if(regular(t.corners[0], t.corners[1], t.corners[2], t.corners[3], p)>0) {
			star.add(t);
			if(t.neighbors[0]!=null && !star.contains(t.neighbors[0])) collectStar(p, t.neighbors[0], star);
			if(t.neighbors[1]!=null && !star.contains(t.neighbors[1])) collectStar(p, t.neighbors[1], star);
			if(t.neighbors[2]!=null && !star.contains(t.neighbors[2])) collectStar(p, t.neighbors[2], star);
			if(t.neighbors[3]!=null && !star.contains(t.neighbors[3])) collectStar(p, t.neighbors[3], star);
		}
	}

	public Tetr walk(PointWeighted p){
		Tetr t = bigTetr;
		double minDistSq = Double.POSITIVE_INFINITY;
		for(Tetr last: lastQueue){
			double dSq = last.corners[1].distanceSquared(p);
			if(dSq<minDistSq){
				minDistSq = dSq;
				t = last;
			}
		}

		bigloop: while(true){
			//			System.out.println(t);
			for(int i=0;i<4;i++){
				int orient = (int)Math.signum(orient(t.corners[(i+1)&3],t.corners[(i+2)&3],t.corners[(i+3)&3],p));
				if(orient!=t.cornerSides[i] && orient!=0){
					t = t.neighbors[i];
					continue bigloop;
				}
			}
			return t;
		}
	}


	private Tetr createBigTetr(){
		double max = 100;
		PointWeighted c0 = new PointWeighted(-max,-max, 0, 0);
		PointWeighted c1 = new PointWeighted( max,-max, 0, 0);
		PointWeighted c2 = new PointWeighted( 0, max,-max, 0);
		PointWeighted c3 = new PointWeighted( 0, max, max, 0);
		Tetr ret = new Tetr(c0,c1,c2,c3);
		return ret;
	}

	public boolean bigTet(Tetr t){
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				if(t.corners[i]==bigTetr.corners[j])
					return true;
			}
		}
		return false;
	}

	private final ShewchuckPredicates pred = ShewchuckPredicates.getInstance();
	/** If the result of insphere is positive the point p is inside the circumsphere of t */
	protected double regular(PointWeighted c0, PointWeighted c1, PointWeighted c2, PointWeighted c3, PointWeighted p){
		double insphere = pred.insphere(c0.getCoords(), c1.getCoords(), c2.getCoords(), c3.getCoords(), p.getCoords());
		if(insphere==0) return 0;
		else {
			double orient = pred.orient(c0.getCoords(),c1.getCoords(),c2.getCoords(),c3.getCoords());
			return insphere*orient;
		}
	}

	/** Returns a positive number if q is 'below' the oriented plane spanned by p0, p1 and p2. When looking at the plane 
	 * from 'above' p0, p1 and p2 occur in counterclockwise order. */
	protected double orient(PointWeighted c0, PointWeighted c1, PointWeighted c2, PointWeighted q){
		return pred.orient(c0.getCoords(), c1.getCoords(), c2.getCoords(), q.getCoords());
	}
	
	
//	private final BigDecimalPredicates pred = BigDecimalPredicates.getInstance();
//	private final BigDecimal zero = new BigDecimal(0);
//	
//	/** If the result of insphere is positive the point p is inside the circumsphere of t */
//	protected double regular(PointWeighted c0, PointWeighted c1, PointWeighted c2, PointWeighted c3, PointWeighted p){
//		BigDecimal insphere = pred.insphere(c0.getCoords(), c1.getCoords(), c2.getCoords(), c3.getCoords(), p.getCoords());
//		if(insphere.compareTo(zero)==0) return 0;
//		else {
//			BigDecimal orient = pred.orient(c0.getCoords(),c1.getCoords(),c2.getCoords(),c3.getCoords());
//			return insphere.multiply(orient).compareTo(zero);
//		}
//	}
//
//	/** Returns a positive number if q is 'below' the oriented plane spanned by p0, p1 and p2. When looking at the plane 
//	 * from 'above' p0, p1 and p2 occur in counterclockwise order. */
//	protected double orient(PointWeighted c0, PointWeighted c1, PointWeighted c2, PointWeighted q){
//		return pred.orient(c0.getCoords(), c1.getCoords(), c2.getCoords(), q.getCoords()).compareTo(zero);
//	}
}
