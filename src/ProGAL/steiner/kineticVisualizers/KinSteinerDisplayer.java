package ProGAL.steiner.kineticVisualizers;

import java.awt.Color;
import java.awt.event.MouseEvent;

import ProGAL.steiner.bnb.MinimumSpanningTree;
import ProGAL.steiner.bnb.PointPlacementWS;
import ProGAL.steiner.bnb.SteinerBnB;
import ProGAL.steiner.bnb.Topology;
import ProGAL.geom2d.Circle;
import ProGAL.geom2d.LineSegment;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.Shape;
import ProGAL.geom2d.Vector;
import ProGAL.geom2d.viewer.ClickListener;
import ProGAL.geom2d.viewer.J2DScene;

public class KinSteinerDisplayer implements ClickListener{
	private final J2DScene scene;
	private final SteinerBnB bnb;
	private final PointPlacementWS pp;
	private final int N;

	private volatile Point[] scSites;
	private volatile Point[] scPoints;
	private final double siteRad;

	public KinSteinerDisplayer(ProGAL.geomNd.Point[] sites){
		this.scene = J2DScene.createJ2DSceneInFrame();
		this.bnb = new SteinerBnB(sites);
		this.pp = new PointPlacementWS(sites);
		this.N = sites.length;

		Point max = new Point(-1000,-1000);
		Point min = new Point(1000,1000);
		for(int i=0;i<N;i++){
			if(sites[i].get(0)>max.get(0)) max.set(0,sites[i].get(0));
			if(sites[i].get(1)>max.get(1)) max.set(1,sites[i].get(1));
			if(sites[i].get(0)<min.get(0)) min.set(0,sites[i].get(0));
			if(sites[i].get(1)<min.get(1)) min.set(1,sites[i].get(1));
		}
		siteRad = Math.max(max.get(0)-min.get(0), max.get(1)-min.get(1))/100.0;

		scSites = new Point[N];
		scPoints = new Point[2*N-2];
		for(int i=0;i<N;i++){
			scSites[i] = new Point(sites[i].getCoords());
			scPoints[i] = scSites[i];
		}

		for(int i=N;i<N+N-2;i++)
			scPoints[i] = new Point();

		updateScene();
		scene.centerCamera();
		scene.autoZoom();
		scene.addClickListener(this);
	}


	private java.awt.Point clickedPoint = null;
	private volatile Circle clickedCircle = null;

	@Override
	public void shapeClicked(Shape shape, final MouseEvent event) {
		if(clickedCircle!=null){ 
			clickedCircle = null;
			return;
		}

		if(shape instanceof Circle)
			clickedCircle = (Circle)shape;
		else return;

		Runnable run = new Runnable(){
			public void run(){
				while(clickedCircle!=null){
					java.awt.Point mouseP = scene.getCanvas().getMousePosition();
					if(mouseP==null){ clickedCircle = null; break; }
					Point mouseLoc = scene.transformPoint(mouseP);
					Point siteLoc = clickedCircle.center();
					double dX = (mouseLoc.x()-siteLoc.x())*0.05;
					double dY = (mouseLoc.y()-siteLoc.y())*0.05;
					clickedCircle.getCenter().addThis(new Vector(dX,dY));
					updateScene();
					try {
						Thread.sleep(30);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
		};
		new Thread(run).start();


		/*
		if(shape!=null){
			if(shape instanceof Circle){
				clickedPoint = event.getPoint();
				clickedCircle = (Circle)shape;
			}else
				clickedPoint = null;
		}
		if(shape==null && clickedPoint!=null){

			final double dX = event.getPoint().x - clickedPoint.x;
			final double dY = event.getPoint().y - clickedPoint.y;

			Runnable run = new Runnable(){
				public void run(){
					//Start moving
					Vector v = new Vector(dX*0.000033/2,-dY*0.000033/2);
					clickedPoint = null;

					for(int i=0;i<200;i++){
						clickedCircle.getCenter().addThis(v);
						updateScene();
						try {
							Thread.sleep(20);
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
					}
				}
			};
			new Thread(run).start();
		}
/**/
	}

	void updateScene(){

		//			for(int i=0;i<N;i++)
		//				scSites[i].set(pp.points[i]);

		MinimumSpanningTree mst = new MinimumSpanningTree(scSites);
		Topology optTop = bnb.solve(N, mst.getLength());
		pp.updateSteinerPoints(optTop, 0.001);

		for(int i=N;i<N+N-2;i++)
			scPoints[i].set(pp.points[i]);

		scene.removeAllShapes();


		//Draw SMT
		for(int e=0;e<optTop.edges.length;e++){
			int p0 = optTop.edges[e][0];
			int p1 = optTop.edges[e][1];
			scene.addShape(new LineSegment(scPoints[p0], scPoints[p1]), new Color(100,100,250));
		}

		//Draw MST
		int[][] mstEdges = mst.getEdges();
		for(int e=0;e<mstEdges.length;e++){
			int p0 = mstEdges[e][0];
			int p1 = mstEdges[e][1];
			Shape edge = new LineSegment(scPoints[p0], scPoints[p1]);
			scene.addShape(edge, new Color(250,100,100), 0.005);

		}

		//Draw sites
		for(int i=0;i<N;i++)
			scene.addShape(new Circle(scSites[i], siteRad), Color.BLACK, 0,true);

		scene.repaint();

	}


	public static void main(String[] args) {
		ProGAL.geomNd.Point[] sites = new ProGAL.geomNd.Point[]{
				new ProGAL.geomNd.Point(new double[]{0,0}),
				new ProGAL.geomNd.Point(new double[]{0,1}),
				new ProGAL.geomNd.Point(new double[]{3,1}),
				new ProGAL.geomNd.Point(new double[]{1.5,1.5}),
				//				new ProGAL.geomNd.Point(new double[]{1.5,-0.5}),
				new ProGAL.geomNd.Point(new double[]{3,0})
		}; 
		new KinSteinerDisplayer(sites);
	}


}
