package ProGAL.geom3d.viewer;

import java.awt.Color;
import java.awt.Font;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsEnvironment;
import java.awt.event.*;
import java.util.ArrayList;
import java.util.ConcurrentModificationException;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Timer;
import java.util.TimerTask;
import java.util.Map.Entry;

import javax.media.j3d.*;
import javax.swing.*;
import javax.vecmath.Color3f;
import javax.vecmath.Matrix3f;
import javax.vecmath.Vector3d;
import javax.vecmath.Vector3f;

import ProGAL.math.Matrix;

import com.sun.j3d.utils.geometry.Primitive;
import com.sun.j3d.utils.geometry.Text2D;
import com.sun.j3d.utils.picking.PickCanvas;
import com.sun.j3d.utils.picking.PickResult;
import com.sun.j3d.utils.picking.PickTool;
import com.sun.j3d.utils.universe.SimpleUniverse;

import ProGAL.geom3d.*;
import ProGAL.geom3d.surface.ParametricParaboloid;
import ProGAL.geom3d.surface.ParametricSurface;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Lens;
import ProGAL.geom3d.volumes.OBB;
import ProGAL.geom3d.volumes.RSS;
import ProGAL.geom3d.volumes.Sphere;

/** A graphics class for viewing scenes using Java3D. 
 * All the <code>Shape</code>-subclasses specified in the <code>edu.geom3D</code> 
 * package can be added to a <code>J3DScene</code> object and are automatically 
 * painted on a <code>Canvas3D</code> object. For 
 * instance the following code creates a scene with a cylinder and a red 
 * transparent box and adds the canvas to a frame. 
 * <pre>
 * J3DScene scene = new J3DScene();
 * scene.addShape(  new Cylinder(new Vector(1,0,0), new Vector(0.5,0.5, 0.3), 0.1f) );
 * 
 * Vector boxCorner = new Vector(-1,0,0);
 * Vector[] boxBases = {new Vector(1,0,0), new Vector(0,1,0), new Vector(0,0,1)};
 * float[] boxExtents = {0.8f, 1, 2};
 * Box box = new Box( boxCorner, boxBases, boxExtents );
 * scene.addShape( box, new Color(200,0,0,100) );
 * 
 * Canvas3D canvas = scene.getCanvas();
 * 
 * JFrame frame = new JFrame();
 * frame.setSize(400,400);
 * frame.getContentPane().add( canvas );
 * frame.setVisible(true);
 * </pre>
 * Text can be added to the scene as well and will always face the camera. 
 * 
 * The <code>repaint()</code> method must be called every time the position of 
 * shapes has changed and the canvas should be updated. The pointers 
 * to added shapes are stored, so subsequent changes in the <code>box</code> 
 * object in the above code will be visible on the canvas when <code>repaint()</code> 
 * is called. The following example shows how to animate a sphere rotating around origo.
 * <pre>
 *  J3DScene scene = new J3DScene();
 *  Sphere sphere = new Sphere( new Vector(1,0,0), 0.1f); 
 *  scene.addShape(sphere);
 *  float t = 0;
 *  while(true){
 * 		t+=0.01f;
 * 		sphere.center = new Vector(Math.cos(t), Math.sin(t), 0);
 * 		scene.repaint();
 * 		try{ Thread.sleep(30); }catch(InterruptedException exc){}
 *  }
 * </pre>
 * 
 * A static method is supplied for conveniently creating a frame containing a scene-viewer. 
 * The following example shows how to quickly create a <code>J3DScene</code> object 
 * that is shown in a frame and ready for use:
 * <pre>
 * J3DScene scene = J3DScene.createJ3DSceneInFrame();
 * scene.setAxisEnabled(true);
 * scene.addShape(  new Cylinder(new Vector(1,0,0), new Vector0,1,0), 0.1f) );
 * </pre>
 * @author R. Fonseca
 */
public class J3DScene {
	public JFrame frame;
	Canvas3D canvas;
	private BranchGroup sceneRoot, scene;
	//	private CamBehavior camBehavior;
	private RebuildBehavior rebuildBehavior;
	private PickCanvas pickCanvas;
	private Timer repaintTimer;

	private final BoundingSphere bounds = new BoundingSphere(new javax.vecmath.Point3d(0,0,0), 5000);
	private Background background;
	private final LinearFog fog = new LinearFog(); 

	private final Map<Shape,BranchGroup> shapeTransforms = new HashMap<Shape,BranchGroup>();
	final Map<Shape,Color> primitives = new HashMap<Shape,Color>();
	private final Map<Node, Shape> pickMap = new HashMap<Node,Shape>();
	private final List<ClickListener> clickListeners = new LinkedList<ClickListener>();
	private final List<Shape> axisElements = new ArrayList<Shape>();
	private Camera camera;
	//	OrbitBehavior orbitBehavior;
	//	private final Point sceneCenter = new Point(0,0,0);


	/** Set color of background. */
	public void setBackgroundColor(Color c){
		background.setColor(c.getRed()/255f, c.getGreen()/255f, c.getBlue()/255f);
		fog.setColor(c.getRed()/255f, c.getGreen()/255f, c.getBlue()/255f);
	}

	/** Removes one volume from the scene. */
	public void removeShape(Shape v){
		primitives.remove(v);
		BranchGroup bg = shapeTransforms.remove(v);
		if(bg!=null){
			bg.detach();
			scene.removeChild(bg);
		}

		for(Entry<Node,Shape> entry: new LinkedList<Entry<Node,Shape>>(pickMap.entrySet())){
			if(entry.getValue()==v){ 
				pickMap.remove(entry.getKey());
			}
		}
		//		if(camera.getControlPanel()!=null && camera.getControlPanel().isVisible()) 
		//			camera.collectShapes();
	}

	/** Remove all volumes from the scene. */
	public void removeAllShapes(){
		while(!primitives.isEmpty())
			removeShape(primitives.entrySet().iterator().next().getKey());
		primitives.clear();
		shapeTransforms.clear();
		pickMap.clear();
	}


	/** Add a volume object. The standard color gray will be used */
	public void addShape(Shape v){	addShape(v,Color.gray);	}

	/** Add a volume object with a specified color */
	public void addShape(Shape v, Color c){	
		addShape(v,c,12);
	}

	/** Add a volume object with a specified color and detail-level */
	public void addShape(Shape v, Color c, int divisions){	
		primitives.put(v, c);
		Node p = genPrimitive(v, c, divisions);
		if(p!=null){
			scene.addChild(p);

			//			if(camera!=null && camera.getControlPanel()!=null && camera.getControlPanel().isVisible()) 
			//				camera.collectShapes();
		}
	}

	/** Add a text-object at the specified position. */
	public TextShape addText(String t, Point pos){ 
		TextShape text = new TextShape(t,pos);
		addShape(text, Color.GRAY); 
		return text;
	}
	public void addText(String t, Point pos, double height){ 
		addShape(new TextShape(t,pos,height), Color.GRAY); 
	}
	public TextShape addText(String t, Point pos, double height, Color c){
		TextShape text = new TextShape(t,pos,height);
		addShape(text, c);
		return text;
	}


	public void addSurface(ParametricSurface surface){
		addSurface(surface, Color.GRAY, -10, 10, 10, -10, 10, 10);
	}
	public void addSurface(ParametricSurface surface, Color col){
		addSurface(surface, col, -10, 10, 10, -10, 10, 10);
	}
	public void addSurface(ParametricSurface surface, Color col, double uMin, double uMax, int uDivs, double vMin, double vMax, int vDivs){
		primitives.put(surface, col);
		Node p = genSurface(surface,uMin, uMax, uDivs, vMin, vMax, vDivs, col);
		if(p!=null)	scene.addChild(p);
	}

	/** Sets the location that the camera looks at to the center of all the shapes added 
	 * to the scene.  */
	public void centerCamera(){
		Vector newCenter = new Vector(0,0,0);

		if(!primitives.isEmpty()){
			for(Entry<Shape, Color> entry: primitives.entrySet()){
//				System.out.println(entry.getKey()+" .. "+entry.getKey().getCenter().toVector());
				newCenter.addThis(entry.getKey().getCenter().toVector());
			}
			newCenter.multiplyThis(1f/primitives.entrySet().size());
		}
		//		centerCamera(newCenter.toPoint());

		//		Transform3D transform = new Transform3D();
		//		transform.setTranslation(new Vector3f(-(float)newCenter.x(), -(float)newCenter.y(), -(float)newCenter.z()));
		//		TransformGroup tg = ((TransformGroup)((TransformGroup)sceneRoot.getChild(0)).getChild(0));
		//		tg.setTransform(transform);
		//		sceneCenter = newCenter.toPoint();


		camera.setLookingAt(newCenter.toPoint());
	}

	public void centerCamera(Point newCenter){
		Point lookingAt = camera.getLookingAt();
		Line l = new Line(lookingAt.clone(), camera.getLookingAt().vectorTo(newCenter));
		for(double t=0;t<=1;t+=(Math.sin(t*Math.PI)/10+0.01)){
			lookingAt.set(l.getPoint(t));
			camera.updateView();
			//			float x = (float)( sceneCenter.x()*(1-t) + newCenter.x()*t );
			//			float y = (float)( sceneCenter.y()*(1-t) + newCenter.y()*t );
			//			float z = (float)( sceneCenter.z()*(1-t) + newCenter.z()*t );
			//			Transform3D transform = new Transform3D();
			//			transform.setTranslation(new Vector3f(-x,-y,-z));
			//			TransformGroup tg = ((TransformGroup)((TransformGroup)sceneRoot.getChild(0)).getChild(0));
			//			tg.setTransform(transform);
			try{Thread.sleep(50);}catch(InterruptedException exc){}
		}
		//		Transform3D transform = new Transform3D();
		//		transform.setTranslation(new Vector3f(-(float)newCenter.x(), -(float)newCenter.y(), -(float)newCenter.z()));
		//		TransformGroup tg = ((TransformGroup)((TransformGroup)sceneRoot.getChild(0)).getChild(0));
		//		tg.setTransform(transform);
		//		sceneCenter.set(newCenter);
	}

	/** Zooms such that the maximal distance between two objects is within the view */
	public void autoZoom(){
		if(primitives.isEmpty()) return;
		//View axis
		Line l = new Line(camera.getEye(), camera.getEye().vectorTo(camera.getLookingAt()).normalizeThis());
		double tanAlpha = Math.tan(0.8*camera.getViewAngle()/2);
		try{
			double minT = Double.POSITIVE_INFINITY;
			for(Entry<Shape, Color> entry: primitives.entrySet()){
				Point p = entry.getKey().getCenter();
				double tProj = l.orthogonalProjectionParameter(p);
				double h = p.distance(l.getPoint(tProj));
				double t = -h/tanAlpha+tProj;
				if(t<minT) minT = t;
			}
			camera.setLocation(l.getPoint(minT));
		}catch(ConcurrentModificationException exc){
			try{ Thread.sleep(300); }catch(InterruptedException exc2){}
			autoZoom();
			return;
		}

		//				double maxDist = 0;
		//				try{
		//					for(Entry<Shape, Color> entry: primitives.entrySet()){
		//						for(Entry<Shape, Color> entry2: primitives.entrySet()){
		//							double d = entry.getKey().getCenter().distance(entry2.getKey().getCenter());
		//							if(d>maxDist) maxDist=d;
		//						}
		//					}
		//				}catch(ConcurrentModificationException exc){
		//					try{ Thread.sleep(300); }catch(InterruptedException exc2){}
		//					autoZoom();
		//					return;
		//				}
		//				if(maxDist>0){
		//					this.camBehavior.setScale(4/(maxDist+10));
		//					this.repaint();
		//				}
	}

	private boolean parallelProjection = false;

	/** Enables and disables parallel projection (as opposed to perspective projection). */
	public void setParallelProjection(boolean enabled) {
		if(enabled && !parallelProjection){
			canvas.getView().setProjectionPolicy(View.PARALLEL_PROJECTION);
		}
		if(!enabled && parallelProjection){
			canvas.getView().setProjectionPolicy(View.PERSPECTIVE_PROJECTION);
		}
		parallelProjection = enabled;
	}


	public void setAntialiasing(boolean enabled){
		canvas.getView().setSceneAntialiasingEnable(enabled);
	}

	private class RotThread extends Thread {
		boolean stop = false;
		public void run() {
			stop = false;
			while(!stop){
				//				camBehavior.rotate(0.01f);
				double angle = 0.01;
				Point eye = camera.getEye();
				Point center = camera.getLookingAt();
				Vector up = camera.getUp();
				Vector x = eye.vectorTo(center);
				Vector y = x.cross(up);
				eye = eye.addThis(x.multiplyThis(1-Math.cos(angle))).addThis(y.multiplyThis(Math.sin(angle)));
				camera.updateView();
				try {Thread.sleep(40);} catch (InterruptedException e) {	}
			}
		}
	}

	private RotThread rotThread;

	/** Toggles rotation */
	public void toggleRotation(){
		//Thread t = new RotThread();
		//t.start();
		if(rotThread!=null && rotThread.isAlive()){
			rotThread.stop = true;
		}else{
			rotThread = new RotThread();
			rotThread.start();
		}
	}

	/**
	 * Add a click-listener that gets called every time an object or the background is clicked
	 * @param cl
	 */
	public void addClickListener(ClickListener cl){
		clickListeners.add(cl);
	}
	public List<ClickListener> getClickListeners(){
		return clickListeners;
	}


	private void updateTransforms(Shape v){
		if(v instanceof ProGAL.geom3d.volumes.Sphere) updateSphereTransforms((ProGAL.geom3d.volumes.Sphere)v);
		//		if(v instanceof ProGAL.geom3d.volumes.Cylinder) updateCylinderTransforms((ProGAL.geom3d.volumes.Cylinder)v);
		//		if(v instanceof ProGAL.geom3d.volumes.Box) updateBoxTransforms((ProGAL.geom3d.volumes.Box)v);
		//		if(v instanceof ProGAL.geom3d.volumes.Cone) updateConeTransforms((ProGAL.geom3d.volumes.Cone)v);
		if(v instanceof ProGAL.geom3d.volumes.RSS) 	updateRSSTransforms((ProGAL.geom3d.volumes.RSS)v);
		if(v instanceof ProGAL.geom3d.volumes.LSS) 	updateLSSTransforms((ProGAL.geom3d.volumes.LSS)v);
		if(v instanceof ProGAL.geom3d.volumes.Tetrahedron) updateTetrahedronTransforms((ProGAL.geom3d.volumes.Tetrahedron)v);
		if(v instanceof ProGAL.geom3d.Triangle) updateTriangleTransforms((ProGAL.geom3d.Triangle)v);
		if(v instanceof TextShape) updateTextTransforms((TextShape)v);
		if(v instanceof ParametricSurface) updateSurface((ParametricSurface)v);
	}
	private void updateSphereTransforms(ProGAL.geom3d.volumes.Sphere s){
		TransformGroup tg = (TransformGroup)shapeTransforms.get(s).getChild(0);

		Transform3D trans = new Transform3D();
		trans.setTranslation(toJ3DVec(s.getCenter()));
		trans.setScale(s.getRadius());

		tg.setTransform(trans);
	}
//	private void updateTorusTransforms(ProGAL.geom3d.volumes.Torus t){
//		Transform3D trans = new Transform3D();
//		Vector v1 = new Vector(0,0,1.0001);
//		Vector v2 = t.getNormal();
//
//		if(v2.length()>0.000001 && v1.angle(v2)>0.00001 && v1.angle(v2)<Math.PI-0.00001){ 
//			Vector v = v1.cross(v2);
//			v.normalizeThis();
//			Matrix m4 = Matrix.createRotationMatrix(v1.angle(v2), v);
//			trans.set(to4x4CoordArray(m4));
//		}
//		trans.setScale(new Vector3d(c.getRadius(), v2.length(), c.getRadius()));
//		trans.setTranslation(toJ3DVec(c.getSegment().getMidPoint().toVector()));
//
//		((TransformGroup)shapeTransforms.get(c).getChild(0)).setTransform(trans);
//	}
	private void updateLensTransforms(ProGAL.geom3d.volumes.Lens lens){
		Transform3D trans = new Transform3D();
		Vector v1 = new Vector(0,1,0);
		Vector v2 = lens.getSphereCenter(0).vectorTo(lens.getSphereCenter(1));
		double angle = v1.angle(v2);
		if(v1.length()>0 && v2.length()>0 && angle>0.00001 && angle<Math.PI-0.00001){
			Matrix m = Matrix.createRotationMatrix(angle, v1.cross(v2).scaleToLength(1));
			trans.set(to4x4CoordArray(m));
		}
		trans.setScale(new Vector3d(lens.getRadius(), 1, lens.getRadius()));
		trans.setTranslation(toJ3DVec(lens.getCenter()));


		TransformGroup tg = ((TransformGroup)shapeTransforms.get(lens).getChild(0));
		tg.setTransform(trans);

		BranchGroup bg = (BranchGroup)tg.getChild(0);

		trans = new Transform3D();
		trans.setScale(new Vector3d(1,lens.getSphereRadius(0)-lens.getFocalDistance(0),1));
		//		System.out.println(lens.getSphereRadius(0)+" "+ lens.getFocalDistance(0));
		//		trans.setTranslation(new Vector3d(0,0,0));
		((TransformGroup)bg.getChild(0)).setTransform(trans);

		trans = new Transform3D();
		trans.rotX(Math.PI);
		//		trans.setTranslation(new Vector3d(0,0,0));
		trans.setScale(new Vector3d(1,lens.getSphereRadius(1)-lens.getFocalDistance(1),1));
		((TransformGroup)bg.getChild(1)).setTransform(trans);
	}
	private void updateCylinderTransforms(ProGAL.geom3d.volumes.Cylinder c){

		Transform3D trans = new Transform3D();
		Vector v1 = new Vector(0,1.0001,0);
		Vector v2 = c.getSegment().getAToB();

		if(v2.length()>0.000001 && v1.angle(v2)>0.00001 && v1.angle(v2)<Math.PI-0.00001){ 
			Vector v = v1.cross(v2);
			v.normalizeThis();
			Matrix m4 = Matrix.createRotationMatrix(v1.angle(v2), v);
			trans.set(to4x4CoordArray(m4));
		}
		trans.setScale(new Vector3d(c.getRadius(), v2.length(), c.getRadius()));
		trans.setTranslation(toJ3DVec(c.getSegment().getMidPoint().toVector()));

		((TransformGroup)shapeTransforms.get(c).getChild(0)).setTransform(trans);
	}
	private void updateConeTransforms(ProGAL.geom3d.volumes.Cone c){

		Transform3D trans = new Transform3D();
		Vector v1 = new Vector(0,-1,0);
		Vector v2 = c.getAxis();//c.p1.vectorTo(c.p2);
		if(v2.length()>0.000001 && v1.angle(v2)>0.00001)
		{ 
			//Matrix m = Matrix.createRotationMatrix(v1.angle(v2), v1.cross(v2).normIn());
			//trans.set(m.getCoordArray());
			Vector v = v1.cross(v2).normalizeThis();
			Matrix m4 = Matrix.createRotationMatrix(v1.angle(v2), v);
//			trans.set(to4x4CoordArray(m4));
			trans.setRotation(new Matrix3f(to3x3CoordArray(m4)));
		}
		trans.setScale(new javax.vecmath.Vector3d(c.getBaseRadius(), c.getAxisLength(), c.getBaseRadius()));
		trans.setTranslation(toJ3DVec(c.getCenter().toVector()));

		((TransformGroup)shapeTransforms.get(c).getChild(0)).setTransform(trans);
	}
	private void updateTextTransforms(TextShape t){
		Transform3D transform = new Transform3D();
		transform.setTranslation(toJ3DVec(t.pos));
		transform.setScale(4*t.height);

		((TransformGroup)shapeTransforms.get(t).getChild(0)).setTransform(transform);
	}
	private void updateBoxTransforms(ProGAL.geom3d.volumes.OBB b){
		Transform3D transform = new Transform3D();
		if(b.getXDir().cross(b.getYDir()).dot(b.getZDir())<0){
			Matrix m = Matrix.createColumnMatrix(
					b.getXDir().multiply(b.extents[0]), 
					b.getZDir().multiply(b.extents[1]), 
					b.getYDir().multiply(b.extents[2])  );
			transform.set(to4x4CoordArray(m));
			//transform.setScale(new Vector3d(b.getXDir().getLength(), b.getZDir().getLength(), b.getYDir().getLength()));
			transform.setTranslation(toJ3DVec(b.getAnchor()));
		}else{
			Matrix m = Matrix.createColumnMatrix(
					b.getXDir().multiply(b.extents[0]), 
					b.getYDir().multiply(b.extents[1]), 
					b.getZDir().multiply(b.extents[2])  );
			transform.set(to4x4CoordArray(m));
			//transform.setScale(new Vector3d(b.getXDir().getLength(), b.getYDir().getLength(), b.getZDir().getLength()));
			transform.setTranslation(toJ3DVec(b.getAnchor()));
		}

		((TransformGroup)shapeTransforms.get(b).getChild(0)).setTransform(transform);
	}
	private void updateTriangleTransforms(ProGAL.geom3d.Triangle t){
		Transform3D transform = new Transform3D();
		Matrix m = Matrix.createColumnMatrix(t.getP1().vectorTo(t.getP2()), t.getP1().vectorTo(t.getP3()), t.getNormal());

		transform.set(to4x4CoordArray(m));
		transform.setTranslation(toJ3DVec(t.getP1()));

		((TransformGroup)shapeTransforms.get(t).getChild(0)).setTransform(transform);
	}
	private void updateTetrahedronTransforms(ProGAL.geom3d.volumes.Tetrahedron t){
		Transform3D transform = new Transform3D();
		Matrix m = Matrix.createColumnMatrix(
				t.getPoint(0).vectorTo(t.getPoint(1)), 
				t.getPoint(0).vectorTo(t.getPoint(2)), 
				t.getPoint(0).vectorTo(t.getPoint(3)));
		transform.set(to4x4CoordArray(m));
		transform.setTranslation(toJ3DVec(t.getPoint(0)));

		((TransformGroup)shapeTransforms.get(t).getChild(0)).setTransform(transform);
	}
	private void updateRSSTransforms(RSS r){
		Transform3D trans = new Transform3D();
		Vector v2 = r.rectangle.bases[0];
		Vector v3 = r.rectangle.bases[1];
		double width = v2.length()*2;
		double height = v3.length()*2;
		double radius = r.radius;

		Matrix m = Matrix.createColumnMatrix(v2.normalize(),v3.normalize(),v2.cross(v3).normalizeThis());
		trans.set(to4x4CoordArray(m));
		trans.setScale(new javax.vecmath.Vector3d(width, height, radius));
		trans.setTranslation(toJ3DVec(r.getCenter()));

		((TransformGroup)shapeTransforms.get(r).getChild(0)).setTransform(trans);
		TransformGroup tg = ((TransformGroup)shapeTransforms.get(r).getChild(0));

		BranchGroup bg = (BranchGroup)tg.getChild(0);


		trans = new Transform3D();
		trans.setScale(new javax.vecmath.Vector3d(radius/width,radius/height,1));
		trans.setTranslation(new javax.vecmath.Vector3d(0.5,0.5,0));
		((TransformGroup)bg.getChild(0)).setTransform(trans);

		trans = new Transform3D();
		trans.rotZ(Math.PI/2);
		trans.setTranslation(new javax.vecmath.Vector3d(-0.5,0.5,0));
		trans.setScale(new javax.vecmath.Vector3d(radius/height,radius/width,1));
		((TransformGroup)bg.getChild(1)).setTransform(trans);

		trans = new Transform3D();
		trans.rotZ(Math.PI);
		trans.setTranslation(new javax.vecmath.Vector3d(-0.5,-0.5,0));
		trans.setScale(new javax.vecmath.Vector3d(radius/width,radius/height,1));
		((TransformGroup)bg.getChild(2)).setTransform(trans);

		trans = new Transform3D();
		trans.rotZ(3*Math.PI/2);
		trans.setTranslation(new javax.vecmath.Vector3d(0.5,-0.5,0));
		trans.setScale(new javax.vecmath.Vector3d(radius/height,radius/width,1));
		((TransformGroup)bg.getChild(3)).setTransform(trans);

		//Half-cylinders
		trans = new Transform3D();
		//		trans.rotZ(3*Math.PI/2);
		trans.setTranslation(new javax.vecmath.Vector3d(0.5,0,0));
		trans.setScale(new javax.vecmath.Vector3d(radius/width,1,1));
		((TransformGroup)bg.getChild(4)).setTransform(trans);

		trans = new Transform3D();
		trans.rotZ(Math.PI/2);
		trans.setTranslation(new javax.vecmath.Vector3d(0,0.5,0));
		trans.setScale(new javax.vecmath.Vector3d(radius/height,1,1));
		((TransformGroup)bg.getChild(5)).setTransform(trans);

		trans = new Transform3D();
		trans.rotZ(Math.PI);
		trans.setTranslation(new javax.vecmath.Vector3d(-0.5,0,0));
		trans.setScale(new javax.vecmath.Vector3d(radius/width,1,1));
		((TransformGroup)bg.getChild(6)).setTransform(trans);

		trans = new Transform3D();
		trans.rotZ(3*Math.PI/2);
		trans.setTranslation(new javax.vecmath.Vector3d(0,-0.5,0));
		trans.setScale(new javax.vecmath.Vector3d(radius/height,1,1));
		((TransformGroup)bg.getChild(7)).setTransform(trans);
	}
	private void updateLSSTransforms(ProGAL.geom3d.volumes.LSS c){
		Transform3D trans = new Transform3D();
		Vector v1 = new Vector(0,1,0);
		Vector v2 = c.segment.getAToB();
		double angle = v1.angle(v2);
		if(v1.length()>0 && v2.length()>0 && angle>0.00001 && angle<Math.PI-0.00001){
			Matrix m = Matrix.createRotationMatrix(angle, v1.cross(v2).scaleToLength(1));
			trans.set(to4x4CoordArray(m));
		}
		trans.setScale(new Vector3d(c.rad, v2.length(), c.rad));
		trans.setTranslation(toJ3DVec(c.segment.getMidPoint()));


		((TransformGroup)shapeTransforms.get(c).getChild(0)).setTransform(trans);
		TransformGroup tg = ((TransformGroup)shapeTransforms.get(c).getChild(0));

		BranchGroup bg = (BranchGroup)tg.getChild(0);

		trans = new Transform3D();
		trans.setScale(new Vector3d(1,c.rad/v2.length(),1));
		trans.setTranslation(new Vector3d(0,0.5,0));
		((TransformGroup)bg.getChild(0)).setTransform(trans);

		trans = new Transform3D();
		trans.rotX(Math.PI);
		trans.setTranslation(new Vector3d(0,-0.5,0));
		trans.setScale(new Vector3d(1,c.rad/v2.length(),1));
		((TransformGroup)bg.getChild(1)).setTransform(trans);
	}

	private void updateSurface(ProGAL.geom3d.surface.ParametricSurface s){
		Surface3D surf3d = (Surface3D)((TransformGroup)shapeTransforms.get(s).getChild(0)).getChild(0);
		surf3d.update();
	}

	private Appearance genAppearance(Color color){
		Appearance app = new Appearance(); 

		Material mat = new Material();
		Color3f ambient  = new Color3f(0.2f,0.2f,0.2f);
		Color3f emissive = new Color3f(0,0,0);
		Color3f diffuse  = new Color3f(color);
		Color3f specular = new Color3f(1,1,1);
		float shininess = 64;
		mat = new Material(ambient, emissive, diffuse, specular, shininess);
		app.setMaterial(mat);

		PolygonAttributes pa = new PolygonAttributes();
		app.setPolygonAttributes(pa);

		if(color.getAlpha()<255){
			TransparencyAttributes ta = new TransparencyAttributes();
			//			ta.setTransparencyMode(TransparencyAttributes.NICEST);
			ta.setTransparencyMode(TransparencyAttributes.FASTEST);
			ta.setTransparency(1-color.getAlpha()/255f);
			app.setTransparencyAttributes(ta);
		}
		return app;
	}



	private Node genTriangle(ProGAL.geom3d.Triangle t, Color color){
		Appearance app = genAppearance(color);

		PolygonAttributes pa = new PolygonAttributes();
		pa.setCullFace(PolygonAttributes.CULL_NONE);
		pa.setBackFaceNormalFlip(true);
		app.setPolygonAttributes(pa);

		Triangle3D tri = new Triangle3D(1, app);

		pickMap.put(tri, t);

		TransformGroup tg = new TransformGroup();
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg.addChild(tri);
		//shapeTransforms.put(t, tg);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		ret.compile();
		shapeTransforms.put(t, ret);
		updateTriangleTransforms(t);
		ret.compile();
		return ret;
	}
	private Node genTetrahedron(ProGAL.geom3d.volumes.Tetrahedron t, Color color){
		Appearance app = genAppearance(color);

		PolygonAttributes pa = new PolygonAttributes();
		pa.setCullFace(PolygonAttributes.CULL_NONE);
		pa.setBackFaceNormalFlip(true);
		app.setPolygonAttributes(pa);

		List<ProGAL.geom3d.Triangle> tList = new LinkedList<ProGAL.geom3d.Triangle>();
		tList.add(new Triangle(new Point(0,0,0),new Point(1,0,0),new Point(0,1,0)));//tList.add(new Triangle(t.p1,t.p2,t.p3));
		tList.add(new Triangle(new Point(0,0,0),new Point(0,0,1),new Point(1,0,0)));//tList.add(new Triangle(t.p1,t.p4,t.p2));
		tList.add(new Triangle(new Point(0,1,0),new Point(0,0,0),new Point(0,0,1)));//tList.add(new Triangle(t.p3,t.p1,t.p4));
		tList.add(new Triangle(new Point(0,1,0),new Point(1,0,0),new Point(0,0,1)));//tList.add(new Triangle(t.p3,t.p2,t.p4));
		TriangleSet3D tset = new TriangleSet3D(tList, app);

		pickMap.put(tset, t);

		TransformGroup tg = new TransformGroup();
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg.addChild(tset);
		//shapeTransforms.put(t, tg);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		ret.compile();
		shapeTransforms.put(t, ret);
		updateTetrahedronTransforms(t);
		return ret;
	}
	private Node genRSS(RSS r, Color color, int divisions){
		Appearance app = genAppearance(color);

		BranchGroup capsGroup = new BranchGroup();
		capsGroup.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		Node shape;

		//Quarterspheres
		Transform3D trans = new Transform3D();
		trans.setTranslation(new javax.vecmath.Vector3d(0.5,0.5,0));
		TransformGroup tg  = new TransformGroup(trans);
		shape = new QuarterSphere3D(1,app,divisions);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		trans = new Transform3D();
		trans.rotZ(Math.PI/2);
		trans.setTranslation(new javax.vecmath.Vector3d(-0.5,0.5, 0));
		tg = new TransformGroup(trans);
		shape = new QuarterSphere3D(1,app,divisions);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		trans = new Transform3D();
		trans.rotZ(Math.PI);
		trans.setTranslation(new javax.vecmath.Vector3d(-0.5,-0.5, 0));
		tg = new TransformGroup(trans);
		shape = new QuarterSphere3D(1,app,divisions);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		trans = new Transform3D();
		trans.rotZ(3*Math.PI/2);
		trans.setTranslation(new javax.vecmath.Vector3d(0.5,-0.5, 0));
		tg = new TransformGroup(trans);
		shape = new QuarterSphere3D(1,app,divisions);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);


		//Half-cylinders
		trans = new Transform3D();
		//trans.rotZ(3*Math.PI/2);
		//		trans.setTranslation(new Vector3d(0.5,0, 0));
		tg = new TransformGroup(trans);
		shape = new HalfCylinder3D(1,1,app, divisions);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		trans = new Transform3D();
		tg = new TransformGroup(trans);
		shape = new HalfCylinder3D(1,1,app, divisions);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		trans = new Transform3D();
		tg = new TransformGroup(trans);
		shape = new HalfCylinder3D(1,1,app, divisions);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		trans = new Transform3D();
		tg = new TransformGroup(trans);
		shape = new HalfCylinder3D(1,1,app, divisions);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		//Bounding planes
		trans = new Transform3D();
		trans.setTranslation(new javax.vecmath.Vector3d(0,0,1));
		tg = new TransformGroup(trans);
		shape = new Rectangle3D(1,1,app);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		capsGroup.addChild(tg);

		trans = new Transform3D();
		trans.rotX(Math.PI);
		trans.setTranslation(new javax.vecmath.Vector3d(0,0,-1));
		tg = new TransformGroup(trans);
		shape = new Rectangle3D(1,1,app);
		tg.addChild(shape);
		enablePicking(shape);
		pickMap.put(shape, r);
		capsGroup.addChild(tg);



		TransformGroup tg1 = new TransformGroup();
		tg1.addChild(capsGroup);
		tg1.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg1.setCapability(TransformGroup.ALLOW_CHILDREN_READ);
		//shapeTransforms.put(c, tg1);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg1);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		shapeTransforms.put(r, ret);

		updateRSSTransforms(r);
		return ret;
	}
	private Node genLSS(ProGAL.geom3d.volumes.LSS c, Color color, int divisions){
		Appearance app = genAppearance(color);

		BranchGroup capsGroup = new BranchGroup();
		capsGroup.setCapability(BranchGroup.ALLOW_CHILDREN_READ);

		//First hemisphere
		Transform3D trans = new Transform3D();
		trans.setTranslation(new Vector3d(0,0.5,0));
		TransformGroup tg  = new TransformGroup(trans);
		Shape3D shape = new Hemisphere3D(1,Math.PI/2,app,divisions);
		enablePicking(shape);
		pickMap.put(shape, c);
		tg.addChild(shape);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		//Second hemisphere
		trans = new Transform3D();
		trans.rotX(Math.PI);
		trans.setTranslation(new Vector3d(0,-0.5, 0));
		tg = new TransformGroup(trans);
		shape = new Hemisphere3D(1,Math.PI/2,app,divisions);
		enablePicking(shape);
		pickMap.put(shape, c);
		tg.addChild(shape);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		shape = new HollowCylinder3D(1,1,app, divisions);
		enablePicking(shape);
		pickMap.put(shape, c);
		capsGroup.addChild(shape);


		TransformGroup tg1 = new TransformGroup();
		tg1.addChild(capsGroup);
		tg1.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg1.setCapability(TransformGroup.ALLOW_CHILDREN_READ);
		//shapeTransforms.put(c, tg1);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg1);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		shapeTransforms.put(c, ret);

		updateLSSTransforms(c);
		ret.compile();
		return ret;
	}
	private Node genCylinder(ProGAL.geom3d.volumes.Cylinder c, Color color, int divisions){
		Appearance app = genAppearance(color);

		com.sun.j3d.utils.geometry.Cylinder cyl = new com.sun.j3d.utils.geometry.Cylinder(1, 1,com.sun.j3d.utils.geometry.Cylinder.GENERATE_NORMALS, divisions, 1, app);
		enablePicking(cyl);
		pickMap.put(cyl, c);

		TransformGroup tg1 = new TransformGroup();
		tg1.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg1.addChild(cyl);
		//shapeTransforms.put(c, tg1);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg1);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		shapeTransforms.put(c, ret);
		updateCylinderTransforms(c);
		return ret;
	}
	//	private Node genCylinder(geom3d.Cylinder3d c, Color color){
	//		Appearance app = genAppearance(color);
	//
	//		//Cylinder cyl = new Cylinder(c.rad, c.p1.distance(c.p2), app);
	//		Cylinder cyl = new Cylinder(1, 1,Cylinder.GENERATE_NORMALS, 32, 1, app);
	//
	//
	//
	//		TransformGroup tg1 = new TransformGroup();
	//		tg1.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
	//		tg1.addChild(cyl);
	//		//shapeTransforms.put(c, tg1);
	//
	//		BranchGroup ret = new BranchGroup();
	//		ret.addChild(tg1);
	//		ret.setCapability(BranchGroup.ALLOW_DETACH);
	//		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
	//		shapeTransforms.put(c, ret);
	//		updateCylinderTransforms(c);
	//		return ret;
	//	}
	private Node genCone(ProGAL.geom3d.volumes.Cone c, Color color, int divisions){
		Appearance app = genAppearance(color);

		//Cylinder cyl = new Cylinder(c.rad, c.p1.distance(c.p2), app);
		com.sun.j3d.utils.geometry.Cone cone = new com.sun.j3d.utils.geometry.Cone(
				1,//Radius
				1,//Height
				com.sun.j3d.utils.geometry.Cone.GENERATE_NORMALS,//primflags
				divisions,//xDivisions
				1,//yDivisions
				app);
		enablePicking(cone);
		pickMap.put(cone, c);

		TransformGroup tg1 = new TransformGroup();
		tg1.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg1.addChild(cone);
		//shapeTransforms.put(c, tg1);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg1);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		ret.compile();
		shapeTransforms.put(c, ret);
		updateConeTransforms(c);
		return ret;
	}
	private Node genSphere(ProGAL.geom3d.volumes.Sphere s, Color color, int divisions){
		Appearance app = genAppearance(color);

		com.sun.j3d.utils.geometry.Sphere sphere = 
				new com.sun.j3d.utils.geometry.Sphere(
						1,
						com.sun.j3d.utils.geometry.Sphere.GENERATE_NORMALS, 
						divisions,
						app);

		enablePicking(sphere);
		pickMap.put(sphere, s);

		/*Transform3D trans = new Transform3D();
		trans.setTranslation(toJ3DVec(s.center));
		trans.setScale(s.radius);*/

		//TransformGroup tg1 = new TransformGroup(trans);
		TransformGroup tg1 = new TransformGroup();
		tg1.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg1.addChild(sphere);
		//shapeTransforms.put(s, tg1);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg1);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		shapeTransforms.put(s, ret);
		updateSphereTransforms(s);
		ret.compile();
		return ret;
	}
//	private Node genTorus(ProGAL.geom3d.volumes.Torus t, Color color, int divisions){
//		Appearance app = genAppearance(color);
//
//		BranchGroup torusGroup = new BranchGroup();
//		torusGroup.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
//
//		//First hemisphere
//		Transform3D trans = new Transform3D();
//		trans.setTranslation(new Vector3d(0,0.5,0));
//		TransformGroup tg  = new TransformGroup(trans);
//		Shape3D shape = new Torus3D(t.getMajorRadius(), t.getMinorRadius(),app,divisions);
//		enablePicking(shape);
//		pickMap.put(shape, t);
//		tg.addChild(shape);
//		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
//		torusGroup.addChild(tg);
//
//		updateTorusTransforms(t);
//		torusGroup.compile();
//		return torusGroup;
//	}
	private Node genLens(ProGAL.geom3d.volumes.Lens lens, Color color, int divisions){
		Appearance app = genAppearance(color);
		BranchGroup capsGroup = new BranchGroup();
		capsGroup.setCapability(BranchGroup.ALLOW_CHILDREN_READ);

		//First hemisphere
		Transform3D trans = new Transform3D();
		trans.setTranslation(new Vector3d(0,0,0));
		TransformGroup tg  = new TransformGroup(trans);
		float angle = (float)(Math.PI/2+Math.atan(lens.getRadius()/lens.getFocalDistance(0)));
		//		System.out.println(angle);
		Shape3D shape = new Hemisphere3D(1f,angle, app, divisions);//new SphereTop3D(1,angle,app,divisions);
		enablePicking(shape);
		pickMap.put(shape, lens);
		tg.addChild(shape);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		//Second hemisphere
		trans = new Transform3D();
		trans.rotX(Math.PI);
		trans.setTranslation(new Vector3d(0,0, 0));
		tg = new TransformGroup(trans);
		angle = (float)(Math.PI/2+Math.atan(lens.getRadius()/lens.getFocalDistance(1)));
		//		shape = new SphereTop3D(1,angle,app,divisions);
		shape = new Hemisphere3D(1,angle, app,divisions);
		enablePicking(shape);
		pickMap.put(shape, lens);
		tg.addChild(shape);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		capsGroup.addChild(tg);

		TransformGroup tg1 = new TransformGroup();
		tg1.addChild(capsGroup);
		tg1.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg1.setCapability(TransformGroup.ALLOW_CHILDREN_READ);
		//shapeTransforms.put(c, tg1);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg1);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		shapeTransforms.put(lens, ret);

		updateLensTransforms(lens);
		ret.compile();
		return ret;
	}
	private Node genBox(ProGAL.geom3d.volumes.OBB b, Color color) {
		Appearance app = genAppearance(color);

		com.sun.j3d.utils.geometry.Box box = new com.sun.j3d.utils.geometry.Box(1,1,1, app);
		enablePicking(box);
		pickMap.put(box, b);

		TransformGroup tg = new TransformGroup();
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg.addChild(box);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		shapeTransforms.put(b, ret);
		updateBoxTransforms(b);
		return ret;
	}
	private Node genSurface(ProGAL.geom3d.surface.ParametricSurface s, double uMin, double uMax, int uDiv, double vMin, double vMax, int vDiv, Color color){
		Appearance app = genAppearance(color);

		PolygonAttributes pa = new PolygonAttributes();
		pa.setCullFace(PolygonAttributes.CULL_NONE);
		pa.setBackFaceNormalFlip(true);
		app.setPolygonAttributes(pa);

		Surface3D surf3d = new Surface3D(s, uMin, uMax, uDiv, vMin, vMax, vDiv, app);

		pickMap.put(surf3d, s);

		TransformGroup tg = new TransformGroup();
		//		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tg.setCapability(TransformGroup.ALLOW_CHILDREN_READ);
		tg.addChild(surf3d);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		ret.compile();
		shapeTransforms.put(s, ret);
		updateSurface(s);
		return ret;
	}
	//	private Node genPrimitive(Shape v, Color c){
	//		//		if(		v instanceof ProGAL.geom3d.volumes.Cylinder)		return genCylinder((geom3d.Cylinder)v, c);
	//		if(v instanceof ProGAL.geom3d.volumes.Sphere)				return genSphere((ProGAL.geom3d.volumes.Sphere)v, c, 32);
	//		else if(v instanceof ProGAL.geom3d.volumes.LSS)				return genLSS((ProGAL.geom3d.volumes.LSS)v, c, 32);
	//		else if(v instanceof ProGAL.geom3d.volumes.RSS)				return genRSS((ProGAL.geom3d.volumes.RSS)v, c, 32);
	//		else if(v instanceof ProGAL.geom3d.volumes.Cylinder)		return genCylinder((ProGAL.geom3d.volumes.Cylinder)v, c, 32);
	//		else if(v instanceof ProGAL.geom3d.volumes.Cone)			return genCone((ProGAL.geom3d.volumes.Cone)v, c);
	//		else if(v instanceof ProGAL.geom3d.volumes.OBB)				return genBox((ProGAL.geom3d.volumes.OBB)v, c);
	//		else if(v instanceof ProGAL.geom3d.Triangle)				return genTriangle((ProGAL.geom3d.Triangle)v, c);
	//		else if(v instanceof ProGAL.geom3d.volumes.Tetrahedron)		return genTetrahedron((ProGAL.geom3d.volumes.Tetrahedron)v, c);
	//		else if(v instanceof TextShape)		return genText((TextShape)v, c);
	//		else{ System.err.println("Warning: unknown primitive: "+v.getClass().getName()); return null; }
	//	}
	private Node genPrimitive(Shape v, Color c, int divisions){
		//if(v instanceof ProGAL.geom3d.volumes.Cylinder)	return genCylinder((ProGAL.geom3d.volumes.Cylinder)v, c, divisions);
		if(v instanceof ProGAL.geom3d.volumes.Sphere)		return genSphere((ProGAL.geom3d.volumes.Sphere)v, c, divisions);
		else if(v instanceof ProGAL.geom3d.volumes.Lens)		return genLens((ProGAL.geom3d.volumes.Lens)v, c, divisions);
		else if(v instanceof ProGAL.geom3d.volumes.LSS)		return genLSS((ProGAL.geom3d.volumes.LSS)v, c, divisions);
		else if(v instanceof ProGAL.geom3d.volumes.RSS)		return genRSS((ProGAL.geom3d.volumes.RSS)v, c, divisions);
		else if(v instanceof ProGAL.geom3d.volumes.Cylinder)		return genCylinder((ProGAL.geom3d.volumes.Cylinder)v, c, divisions);
		else if(v instanceof ProGAL.geom3d.volumes.Cone)			return genCone((ProGAL.geom3d.volumes.Cone)v, c, divisions);
		else if(v instanceof ProGAL.geom3d.volumes.OBB)				return genBox((ProGAL.geom3d.volumes.OBB)v, c);
		else if(v instanceof ProGAL.geom3d.Triangle)				return genTriangle((ProGAL.geom3d.Triangle)v, c);
		else if(v instanceof ProGAL.geom3d.volumes.Tetrahedron)		return genTetrahedron((ProGAL.geom3d.volumes.Tetrahedron)v, c);
		else if(v instanceof TextShape)		return genText((TextShape)v, c);
		else{ System.err.println("Warning: unknown primitive: "+v.getClass().getName()); return null; }
	}
	private Node genText(TextShape t, Color c){
		//Color c = new Color(255-c.getRed(), 255-c.getGreen(), 255-c.getBlue());
		Shape3D text3D = new Text2D(t.text, new Color3f(c),"Arial",48, Font.BOLD );

		enablePicking(text3D);
		pickMap.put(text3D, t);

		TransformGroup subTg = new TransformGroup();
		subTg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		Billboard billboard = new Billboard(subTg);
		billboard.setSchedulingBounds( bounds );
		subTg.addChild( billboard );

		TransformGroup tg = new TransformGroup();
		tg.addChild(subTg);
		subTg.addChild(text3D);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		//shapeTransforms.put(t, tg);

		BranchGroup ret = new BranchGroup();
		ret.addChild(tg);
		ret.setCapability(BranchGroup.ALLOW_DETACH);
		ret.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		shapeTransforms.put(t, ret);
		updateTextTransforms(t);
		return ret;
	}
	private Node genBackground() {
		background = new Background(new Color3f(Color.white)); 
		background.setCapability(Background.ALLOW_COLOR_WRITE);
		background.setCapability(Background.ALLOW_COLOR_READ);
		background.setApplicationBounds(bounds);
		return background;
	}
	private void genLights(BranchGroup bg){

		Color3f lightColor = new Color3f(1f, 1f, 1f);

		Vector3f lightDirection = new Vector3f(0,-5,-5);
		DirectionalLight light  = new DirectionalLight(lightColor, lightDirection);
		light.setInfluencingBounds(bounds);
		bg.addChild(light);

		lightDirection = new Vector3f(0,-5,5);
		light  = new DirectionalLight(lightColor, lightDirection);
		light.setInfluencingBounds(bounds);
		bg.addChild(light);

		lightDirection = new Vector3f(0,5,0);
		light  = new DirectionalLight(lightColor, lightDirection);
		light.setInfluencingBounds(bounds);
		bg.addChild(light);

	}

	private boolean axisEnabled = false;
	/** Enables or disables xyz-axis from the origo */
	public void setAxisEnabled(boolean axisEnabled){
		if(axisEnabled && axisElements.isEmpty()){
			float rad = 0.02f;
			axisElements.add(new ProGAL.geom3d.volumes.Cylinder(new Point(0,0,0),new Point(1-2*rad,0,0), rad));
			axisElements.add(new ProGAL.geom3d.volumes.Cylinder(new Point(0,0,0),new Point(0,1-2*rad,0), rad));
			axisElements.add(new ProGAL.geom3d.volumes.Cylinder(new Point(0,0,0),new Point(0,0,1-2*rad), rad));
			axisElements.add(new ProGAL.geom3d.volumes.Cone(new Point(1,0,0), new Point(1-2*rad,0,0), 2*rad));
			axisElements.add(new ProGAL.geom3d.volumes.Cone(new Point(0,0,1), new Point(0,0,1-2*rad), 2*rad));
			axisElements.add(new ProGAL.geom3d.volumes.Cone(new Point(0,1,0), new Point(0,1-2*rad,0), 2*rad));
//			axisElements.add(new ProGAL.geom3d.volumes.Cone(new Point(1-2*rad,0,0),new Point(1,0,0), 2*rad));
//			axisElements.add(new ProGAL.geom3d.volumes.Cone(new Point(0,0,1-2*rad),new Point(0,0,1), 2*rad));
//			axisElements.add(new ProGAL.geom3d.volumes.Cone(new Point(0,1-2*rad,0),new Point(0,1,0), 2*rad));

			axisElements.add(new TextShape("x", new Point(1,0,0), 0.3));
			axisElements.add(new TextShape("y", new Point(0,1,0), 0.3));
			axisElements.add(new TextShape("z", new Point(0,0,1), 0.3));
		}
		if(axisEnabled && !this.axisEnabled){
			for(Shape s: axisElements) addShape(s, Color.GRAY);
		}
		if(!axisEnabled && this.axisEnabled){
			for(Shape s: axisElements) removeShape(s);
		}
		this.axisEnabled = axisEnabled;
	}


	private void initialBuild(){
		sceneRoot = new BranchGroup();
		sceneRoot.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		TransformGroup tgroup = new TransformGroup();
		sceneRoot.addChild(tgroup);
		tgroup.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		tgroup.setCapability(TransformGroup.ALLOW_CHILDREN_EXTEND);
		tgroup.setCapability(TransformGroup.ALLOW_CHILDREN_WRITE);
		tgroup.setCapability(TransformGroup.ALLOW_CHILDREN_READ);

		//		camBehavior = new CamBehavior(tgroup);
		//		camBehavior.setSchedulingBounds(bounds);
		//		sceneRoot.addChild(camBehavior);

		//rebuildBehavior = new RebuildBehavior(tgroup);
		rebuildBehavior = new RebuildBehavior();
		rebuildBehavior.setSchedulingBounds(bounds);
		sceneRoot.addChild(rebuildBehavior);


		//BranchGroup scene = buildScene();
		//BranchGroup scene = new BranchGroup();

		Transform3D transform = new Transform3D();
		transform.setTranslation(toJ3DVec(new Vector(0,0,0)));
		TransformGroup tg = new TransformGroup(transform);
		tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);

		scene = new BranchGroup();
		//scene.setCapability(BranchGroup.ALLOW_DETACH);
		scene.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
		scene.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);

		for(Entry<Shape, Color> entry: primitives.entrySet())
			scene.addChild(genPrimitive(entry.getKey(), entry.getValue(), 32));
		//for(TextPrimitive tp: texts) movedScene.addChild(genTextPrimitive(tp));
		genLights(scene);
		scene.addChild(genBackground());
		//if(paintAxis) scene.addChild(genAxis());
		tg.addChild(scene);
		//scene.addChild(tg);
		tgroup.addChild(tg);

		fog.setColor(new Color3f(Color.WHITE));
		fog.setFrontDistance(9);
		fog.setBackDistance(10);
		//	    fog.setCapability(Fog.ALLOW_COLOR_WRITE);
		fog.setCapability(LinearFog.ALLOW_DISTANCE_WRITE);
		fog.setCapability(LinearFog.ALLOW_DISTANCE_READ);
		fog.setCapability(LinearFog.ALLOW_COLOR_WRITE);
		fog.setInfluencingBounds(bounds);
		scene.addChild(fog);

		//tgroup.addChild(scene);
		scene.compile();
		sceneRoot.compile();


	}

	public J3DScene(){
		this.initialBuild();
	}

	public Camera getCamera(){
		return camera;
	}


	//	private class HudCanvas3D extends Canvas3D implements MouseListener{
	//		private static final long serialVersionUID = 1L;
	//
	//		public HudCanvas3D(GraphicsConfiguration arg0) {
	//			super(arg0);
	//			this.addMouseListener(this);
	//		}
	//
	//		public void postRender(){
	//			super.postRender();
	//			J3DGraphics2D g = super.getGraphics2D();
	//			g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
	//			int w = 30;
	//			int h = 30;
	//			int x = getWidth()-w-2;
	//			int y = getHeight()-h-2;
	//
	//			//Fill circle
	//			g.setColor(new Color(200,200,200));
	//			g.fillOval(x, y, w, h);
	//			x+=2;y+=2;w-=4;h-=4;
	//			g.setColor(new Color(100,100,200));
	//			g.fillOval(x, y, w, h);
	//
	//			//Highlights
	//			for(int i=1;i<10;i++){
	//				g.setColor(new Color(200-10*i,200-10*i,250-5*i));
	//				g.fillOval((int)(x+5), y+i, (int)(w-10), h/2-i);
	//			}
	//			for(int i=1;i<6;i++){
	//				g.setColor(new Color(200-20*i,200-20*i,250-10*i));
	//				g.fillOval((int)(x+5), y-i+h/2, (int)(w-10), h/2-i);
	//			}
	//
	//			//Camera
	//			g.setColor(Color.WHITE);
	//			g.fillRoundRect((int)(x+w*0.16), (int)(y+h*0.47), (int)(w*0.5), (int)(h*0.3), (int)(w*0.15), (int)(w*0.15));
	//			g.fillOval((int)(x+w*0.13), (int)(y+h*0.27), (int)(h*0.25), (int)(h*0.25));
	//			g.fillOval((int)(x+w*0.4), (int)(y+h*0.22), (int)(h*0.3), (int)(h*0.3));
	//			g.fillPolygon(
	//					new int[]{(int)(x+w*0.69), (int)(x+w*0.69), (int)(x+w*0.84), (int)(x+w*0.84)}, 
	//					new int[]{(int)(y+h*0.68), (int)(y+h*0.50), (int)(y+h*0.45), (int)(y+h*0.71)}, 4 );
	//
	//			g.flush(false);
	//		}
	//
	//		public void mouseClicked(MouseEvent arg0) {	}
	//		public void mouseEntered(MouseEvent arg0) {}
	//		public void mouseExited(MouseEvent arg0) {	}
	//		public void mousePressed(MouseEvent arg0) {	}
	//		public void mouseReleased(MouseEvent arg0) {
	//			java.awt.Point p = arg0.getPoint();
	//
	//			int w = 30;
	//			int h = 30;
	//			int x = getWidth()-w-2;
	//			int y = getHeight()-h-2;
	//			java.awt.Point buttonCenter = new java.awt.Point(x+w/2, y+h/2);
	//			if(p.distance(buttonCenter)<15){
	//				JFrame ctrlPanel = camera.getControlPanel();
	//				if(ctrlPanel==null) camera.createControlPanel();
	//				ctrlPanel = camera.getControlPanel();
	//				ctrlPanel.setVisible(true);
	//				java.awt.Point loc = arg0.getLocationOnScreen();
	//				loc.x-=ctrlPanel.getWidth();
	//				loc.y-=ctrlPanel.getHeight();
	//				ctrlPanel.setLocation(loc);
	//				camera.collectShapes();
	//			}
	//		}
	//
	//	}

	/** Get the canvas that displays the scene. If this method is called 
	 * several times the same <code>Canvas3D</code> object will be returned 
	 * every time.*/
	public Canvas3D getCanvas(){
		if(canvas==null){
			//initialBuild();

			//			GraphicsConfiguration config = SimpleUniverse.getPreferredConfiguration();
			GraphicsConfigTemplate3D template = new GraphicsConfigTemplate3D();
			template.setSceneAntialiasing(GraphicsConfigTemplate3D.PREFERRED);
			template.setRedSize(6);
			template.setGreenSize(6);
			template.setBlueSize(6);//Fixes the weird ugly rastering on mac 
			GraphicsConfiguration config =
					GraphicsEnvironment.getLocalGraphicsEnvironment().
					getDefaultScreenDevice().getBestConfiguration(template);

			//			canvas = new HudCanvas3D(config);
			canvas = new Canvas3D(config);

			SimpleUniverse universe = new SimpleUniverse(canvas);
			//universe.getViewer().getView().setProjectionPolicy(View.PARALLEL_PROJECTION);
			universe.addBranchGraph(sceneRoot);
			universe.getViewingPlatform().setNominalViewingTransform();
			universe.getViewer().getView().setLocalEyeLightingEnable(true);

			camera = new Camera(this, universe.getViewingPlatform(), fog);
			//						universe.getViewer().getView().setSceneAntialiasingEnable(true);

			//			CamListener cl = new CamListener();
			//			canvas.addMouseListener(cl);
			//			canvas.addMouseMotionListener(cl);
			//			canvas.addKeyListener(cl);
			//			canvas.addMouseWheelListener(cl);

			//			orbitBehavior = new OrbitBehavior(canvas,
			//					OrbitBehavior.PROPORTIONAL_ZOOM | OrbitBehavior.REVERSE_ROTATE
			//					| OrbitBehavior.REVERSE_TRANSLATE );
			//			orbitBehavior.setSchedulingBounds(bounds);    
			//			universe.getViewingPlatform().setViewPlatformBehavior(orbitBehavior);

			pickCanvas = new PickCanvas(canvas, sceneRoot);
			pickCanvas.setMode(PickCanvas.GEOMETRY);
			addClickListener(new ClickListener(){
				public void shapeClicked(Shape shape, MouseEvent e) {
					if(e.getClickCount()==2 && shape!=null){
						centerCamera(shape.getCenter());
						//						camera.setLookingAt(shape.getCenter());
					}
				}});

			canvas.addMouseListener(new PickListener());

			canvas.getView().setTransparencySortingPolicy(View.TRANSPARENCY_SORT_GEOMETRY);

		}

		return canvas;
	}

	/** Repaint the canvas. If the scene has been changed in any way the 
	 * scene displayer will update the view when <code>repaint()</code> is called 
	 * and no sooner. If the scene is repeatedly changed, and repaint repeatedly 
	 * called the viewer will show an animation. */
	public void repaint(){
		rebuildBehavior.rebuild();

		//		if(camera.getControlPanel()!=null && camera.getControlPanel().isVisible()) 
		//			camera.collectShapes();
	}

	/** Repaint the canvas repeatedly every <code>millisecondDelay</code> milliseconds. */
	public void repaintRepeatedly(long millisecondDelay){
		if(repaintTimer!=null){
			repaintTimer.cancel();
		}else{
			repaintTimer = new Timer();
		}
		class RepaintTask extends TimerTask{
			public void run() {
				repaint();
			}
		}
		repaintTimer.schedule(new RepaintTask(), 1, millisecondDelay);
	}

	private static Vector3f toJ3DVec(Vector v){ return new Vector3f((float)v.x(), (float)v.y(), (float)v.z() ); }
	private static Vector3f toJ3DVec(Point v){ return new Vector3f((float)v.x(), (float)v.y(), (float)v.z() ); }
	private static float[] to4x4CoordArray(Matrix m){
		float[] ret = new float[16];
		for(int r=0;r<m.getM();r++) 
			for(int c=0;c<m.getN();c++)
				ret[r*4+c]=(float)m.get(r, c);

		ret[15] = 1; 
		return ret;
	}
	private static float[] to3x3CoordArray(Matrix m){
		float[] ret = new float[9];
		for(int r=0;r<m.getM();r++) 
			for(int c=0;c<m.getN();c++)
				ret[r*3+c]=(float)m.get(r, c);

		return ret;
	}

	private class PickListener extends MouseAdapter{

		@Override
		public void mouseClicked(MouseEvent e) {
			pickCanvas.setShapeLocation(e);
			PickResult result = pickCanvas.pickClosest();
			if(result==null){ 
				for(ClickListener cl: clickListeners) cl.shapeClicked(null, e);
				return;
			}

			Primitive p = (Primitive)result.getNode(PickResult.PRIMITIVE);
			Shape3D s = (Shape3D)result.getNode(PickResult.SHAPE3D);

			Shape clickedShape = p!=null?pickMap.get(p):s!=null?pickMap.get(s):null;
			for(ClickListener cl: clickListeners)
				cl.shapeClicked(clickedShape, e);
		}
	}

	private void enablePicking(Node node) {
		node.setPickable(true);
		node.setCapability(Node.ENABLE_PICK_REPORTING);
		try {
			Group group = (Group) node;
			for (Enumeration<?> e = group.getAllChildren(); e.hasMoreElements();) {
				enablePicking((Node)e.nextElement());
			}
		}catch(ClassCastException e) {
			// if not a group node, there are no children so ignore exception
		}catch(RestrictedAccessException exc){}
		try {
			Shape3D shape = (Shape3D) node;
			PickTool.setCapabilities(node, PickTool.INTERSECT_FULL);
			for (Enumeration<?> e = shape.getAllGeometries(); e.hasMoreElements();) {
				Geometry g = (Geometry)e.nextElement();
				g.setCapability(Geometry.ALLOW_INTERSECT);
			}
		}
		catch(ClassCastException e) {
			// not a Shape3D node ignore exception
		}
	} 

	//		private class CamListener extends MouseAdapter implements MouseMotionListener, MouseWheelListener, KeyListener{
	//			private boolean shiftPressed = false;
	//			private java.awt.Point lastPoint = null;
	//			private long lastTime = System.currentTimeMillis();
	//			public void mousePressed(MouseEvent e) {	lastPoint = e.getPoint();		}
	//			public void mouseReleased(MouseEvent e) {	lastPoint = null; }
	//			public void mouseClicked(MouseEvent e){
	//				rebuildBehavior.rebuild();
	//			}
	//	
	//			public void mouseDragged(MouseEvent e) {
	//				if(lastPoint==null) {
	//					lastPoint = e.getPoint();
	//					lastTime = System.currentTimeMillis();
	//					return;
	//				}
	//				java.awt.Point point = e.getPoint();
	//				float dX = point.x-lastPoint.x;
	//				float dY = point.y-lastPoint.y;
	//				float damper = Math.max(10, (float)(System.currentTimeMillis()-lastTime))*10f;
	//	
	//				if(shiftPressed){
	//					Vector delta = new Vector(dX, -dY, 0).multiplyThis(1/damper);
	//					camBehavior.translate(delta);
	//				}else{
	//					
	//					camBehavior.rotate(dX*(float)Math.PI/damper);
	//				}
	//				lastPoint = point;
	//				lastTime = System.currentTimeMillis();
	//	
	//			}
	//	
	//			public void mouseWheelMoved(MouseWheelEvent e){
	//				float damper = Math.max(10, (float)(System.currentTimeMillis()-lastTime))*10f;
	//				camBehavior.scale(e.getWheelRotation()/damper);
	//				
	////				orbitBehavior.setZoomFactor(orbitBehavior.getZoomFactor()+e.getWheelRotation());
	////				System.out.println(orbitBehavior.getZoomFactor());
	//				lastTime = System.currentTimeMillis();
	//			}
	//	
	//			public void mouseMoved(MouseEvent e) {}
	//			public void keyPressed(KeyEvent e) {
	//				if( e.getKeyCode()==KeyEvent.VK_SHIFT )	shiftPressed = true;
	//	
	//				if(e.getKeyCode()==KeyEvent.VK_DOWN && shiftPressed){
	//					float damper = Math.max(10, (float)(System.currentTimeMillis()-lastTime));
	//					camBehavior.scale(10f/damper);
	//					lastTime = System.currentTimeMillis();
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_UP && shiftPressed){
	//					float damper = Math.max(10, (float)(System.currentTimeMillis()-lastTime));
	//					camBehavior.scale(-10f/damper);
	//					lastTime = System.currentTimeMillis();
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_LEFT && shiftPressed){
	//					camBehavior.rotate(0.1f);
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_RIGHT && shiftPressed){
	//					camBehavior.rotate(-0.1f);
	//				}
	//	
	//				if(e.getKeyCode()==KeyEvent.VK_UP && !shiftPressed){
	//					camBehavior.translate(new Vector(0,-0.1,0));
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_DOWN && !shiftPressed){
	//					camBehavior.translate(new Vector(0,0.1,0));
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_LEFT && !shiftPressed){
	//					camBehavior.translate(new Vector(0.1,0,0));
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_RIGHT && !shiftPressed){
	//					camBehavior.translate(new Vector(-0.1,0,0));
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_S){
	//					J3DImageFileWriter.writeJPEGFile("J3DScene.jpg", canvas);
	//					System.out.println("Stored view to J3DScene.jpg");
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_E){
	//					J3DImageFileWriter.writeEPSFile("J3DScene.eps", canvas);
	//					System.out.println("Stored view to J3DScene.eps");
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_C){
	//					J3DScene.this.centerCamera();
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_Z){
	//					J3DScene.this.autoZoom();
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_R){
	//					J3DScene.this.toggleRotation();
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_P){
	//					J3DScene.this.setParallelProjection(!parallelProjection);
	//				}
	//				if(e.getKeyCode()==KeyEvent.VK_A){
	//					J3DScene.this.setAxisEnabled(!axisEnabled);
	//				}
	//			}
	//			public void keyReleased(KeyEvent e) {
	//				if( e.getKeyCode()==KeyEvent.VK_SHIFT ) shiftPressed = false;
	//			}
	//			public void keyTyped(KeyEvent e) {}
	//	
	//		}

	//	private static class CamBehavior extends Behavior {
	//
	//		private TransformGroup transformGroup;
	//		private Transform3D trans = new Transform3D();
	//		private WakeupCriterion criterion;
	//		private double yAngle = 0.0f;
	//		private Vector3f translation = new Vector3f(0,0,0);
	//		private double scale = 1f;
	//		
	////		private Point3d eye, center;
	////		private Vector3d up;
	//
	//
	//
	//		private final int ROTATE = 1;
	//
	//		// create a new RotateBehavior
	//		CamBehavior(TransformGroup tg) {	
	//			transformGroup = tg;
	////			eye = new Point3d(0,0,1);
	////			center = new Point3d(0,0,0);
	////			up = new Vector3d(0,1,0);
	//		}
	//
	//		// initialize behavior to wakeup on a behavior post with id = ROTATE
	//		public void initialize() {
	//			criterion = new WakeupOnBehaviorPost(this, ROTATE);
	//			wakeupOn(criterion);
	//		}
	//
	//		// processStimulus to rotate the cube
	//		@SuppressWarnings("rawtypes")
	//		public void processStimulus(Enumeration criteria) {
	//			trans.rotY(yAngle);
	//			trans.setTranslation(translation);
	//			trans.setScale(scale);
	////			trans.lookAt(eye, center, up);
	//			transformGroup.setTransform(trans);
	//			wakeupOn(criterion);
	//			//System.out.println("Scale "+scale);
	//		}
	//
	//		// when the mouse is clicked, postId for the behavior
	//		void rotate(float dY) {
	//			yAngle+=dY;
	//			postId(ROTATE);
	//		}
	//		void translate(Vector delta){
	//			translation.add(new Vector3f((float)delta.x(), (float)delta.y(), (float)delta.z()));
	//			postId(ROTATE);
	//		}
	//		void scale(double s){
	//			scale-=s;
	//			if(scale<=0.001) scale=0.001f;
	//			postId(ROTATE);
	//		}
	//		void setScale(double s){
	//			scale=s;
	//			if(scale<=0.001) scale=0.001f;
	//			postId(ROTATE);
	//		}
	//	}

	private class RebuildBehavior extends Behavior {
		private boolean rebuilding = false;
		//private TransformGroup tgroup;
		private WakeupCriterion criterion;

		private final int REBUILD = 5;

		// initialize behavior to wakeup on a behavior post with id = ROTATE
		public void initialize() {
			criterion = new WakeupOnBehaviorPost(this, REBUILD);
			wakeupOn(criterion);
		}

		@SuppressWarnings("rawtypes")
		public void processStimulus(Enumeration criteria) {
			try{
				for(Entry<Shape, BranchGroup> entry: shapeTransforms.entrySet()){
					updateTransforms(entry.getKey());
				}

			}catch(ConcurrentModificationException exc){}
			wakeupOn(criterion);
			rebuilding = false;
		}

		// when the mouse is clicked, postId for the behavior
		synchronized void rebuild() {
			if(rebuilding) return;
			rebuilding = true;
			postId(REBUILD);
		}
	}


	/** 
	 * Create a frame containing a canvas, display it and return the  J3DScene object shown in the frame. 
	 * The frame can be retrieved using the <code>J3DScene.frame</code> field.  
	 */
	public static J3DScene createJ3DSceneInFrame() {
		JFrame f = new JFrame("J3DScene-viewer");
		f.setSize(1000,800);
		JPopupMenu.setDefaultLightWeightPopupEnabled(false);

		J3DScene j3ds = new J3DScene();
		j3ds.frame = f;
		Canvas3D canvas = j3ds.getCanvas();

		JMenuBar menubar = new JMenuBar();
		f.setJMenuBar(menubar);
		JMenu view = new JMenu("View");
		menubar.add(view);

		JMenuItem item;
		view.add(item = new JMenuItem("Export JPG"));
		class ExportActionListener implements ActionListener{
			J3DScene j3ds;
			ExportActionListener(J3DScene j3ds){ this.j3ds = j3ds; }
			public void actionPerformed(ActionEvent e) {
				JFileChooser chooser = new JFileChooser();
				ExampleFileFilter filter = new ExampleFileFilter(new String[]{"jpg","jpeg"}, "JPG images");
				chooser.setFileFilter(filter);
				int returnVal = chooser.showSaveDialog(j3ds.canvas);
				if(returnVal == JFileChooser.APPROVE_OPTION) {
					J3DImageFileWriter.writeJPEGFile(chooser.getSelectedFile().getAbsolutePath(), j3ds.canvas);
				}	
			}	
		}
		item.addActionListener(new ExportActionListener(j3ds));
		view.add(item = new JMenuItem("Export EPS"));
		class ExportEPSActionListener implements ActionListener{
			J3DScene j3ds;
			ExportEPSActionListener(J3DScene j3ds){ this.j3ds = j3ds; }
			public void actionPerformed(ActionEvent e) {
				JFileChooser chooser = new JFileChooser();
				ExampleFileFilter filter = new ExampleFileFilter(new String[]{"eps"}, "Encapsulated Postscript images");
				chooser.setFileFilter(filter);
				int returnVal = chooser.showSaveDialog(j3ds.canvas);
				if(returnVal == JFileChooser.APPROVE_OPTION) {
					J3DImageFileWriter.writeEPSFile(chooser.getSelectedFile().getAbsolutePath(), j3ds.canvas);
				}	
			}
		}
		item.addActionListener(new ExportEPSActionListener(j3ds));
		view.addSeparator();

		view.add(item = new JMenuItem("Auto-zoom (z)"));
		class AutozoomActionListener implements ActionListener{
			J3DScene j3ds;
			AutozoomActionListener(J3DScene j3ds){ this.j3ds = j3ds; }
			public void actionPerformed(ActionEvent e) {
				j3ds.autoZoom();
			}
		}
		item.addActionListener(new AutozoomActionListener(j3ds));

		view.add(item = new JMenuItem("Center view (c)"));
		class CenterActionListener implements ActionListener{
			J3DScene j3ds;
			CenterActionListener(J3DScene j3ds){ this.j3ds = j3ds; }
			public void actionPerformed(ActionEvent e) {
				j3ds.centerCamera();
			}
		}
		item.addActionListener(new CenterActionListener(j3ds));

		view.add(item = new JMenuItem("Toggle rotation (r)"));
		class RotateActionListener implements ActionListener{
			J3DScene j3ds;
			RotateActionListener(J3DScene j3ds){ this.j3ds = j3ds; }
			public void actionPerformed(ActionEvent e) {
				j3ds.toggleRotation();
			}
		}
		item.addActionListener(new RotateActionListener(j3ds));


		view.add(item = new JMenuItem("Toggle parallel projection (p)"));
		class ParallelActionListener implements ActionListener{
			J3DScene j3ds;
			ParallelActionListener(J3DScene j3ds){ this.j3ds = j3ds; }
			public void actionPerformed(ActionEvent e) {
				j3ds.setParallelProjection(!j3ds.parallelProjection);
			}
		}
		item.addActionListener(new ParallelActionListener(j3ds));


		view.add(item = new JMenuItem("Toggle anti-aliasing"));
		class AntialiasActionListener implements ActionListener{
			J3DScene j3ds;
			AntialiasActionListener(J3DScene j3ds){ this.j3ds = j3ds; }
			public void actionPerformed(ActionEvent e) {
				j3ds.setAntialiasing(true);
			}
		}
		item.addActionListener(new AntialiasActionListener(j3ds));


		JMenu help = new JMenu("Help");
		menubar.add(help);

		help.add(item = new JMenuItem("Scene navigation"));
		class HelpActionListener implements ActionListener{
			J3DScene j3ds;
			HelpActionListener(J3DScene j3ds){ this.j3ds = j3ds; }
			public void actionPerformed(ActionEvent e) {
				String message = "Left-click-dragging on the screen will rotate around the 'scene center'.\n\n"
						+ "The 'scene center' can be changed by double-clicking an object, by holding \n"
						+ "shift while left-click-dragging, or by pressing the 'c' button.\n\n"
						+ "The mouse wheel or the 'z' button is used to zoom\n\n";

				JOptionPane.showMessageDialog(j3ds.frame, message, "Navigation", JOptionPane.INFORMATION_MESSAGE);
			}	
		}
		item.addActionListener(new HelpActionListener(j3ds));



		f.getContentPane().add(canvas);
		f.setVisible(true);
		f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		return j3ds;
	}




	public static void main(String[] args) {
		J3DScene j3ds = createJ3DSceneInFrame();
		OBB box = new OBB(
				new Point(-0.5,0.5,-0.5),
				new Vector(-1,-1,0).scaleToLength(0.5),
				new Vector(-1,1,0).scaleToLength(0.2),
				new Vector(0,0,-1).scaleToLength(0.1) );
		j3ds.addShape(box, Color.GREEN.darker());

		j3ds.setAxisEnabled(true);

		j3ds.addShape(new LSS(new Point(1,1,0), new Point(1,0,0), 0.1), new Color(20,200,20, 100));
		ProGAL.geom3d.volumes.Cylinder cyl = new ProGAL.geom3d.volumes.Cylinder(new Point(0.4,0,0.1), new Point(0.4,0.5,0), 0.1);
		j3ds.addShape(cyl, Color.RED.darker().darker());
		ProGAL.geom3d.volumes.Sphere s = new ProGAL.geom3d.volumes.Sphere(new Point(-1,-0.2,0), 0.3f);
		j3ds.addShape(s, Color.MAGENTA);
		ProGAL.geom3d.volumes.Tetrahedron tetr = new ProGAL.geom3d.volumes.Tetrahedron(new Point(0,0,1), new Point(-0.5,0,1), new Point(-0.5,0.5,1), new Point(-0.25,0.25, 1.25));
		j3ds.addShape(tetr, new Color(50,50,255));
		j3ds.addShape(new ProGAL.geom3d.Triangle(new Point(0.2f, -0.2f, 0.1f), new Point(0.8, -0.8, 0.1), new Point(1,-0.3, 0.1)), new Color(150,50,255));
		j3ds.addShape(new RSS(new Point(0.5,0,0),new Vector[]{new Vector(0.23,0.23,0),new Vector(-0.15,0.15,0)},0.1), new Color(200,200,50));

		Sphere s0 = new Sphere(new Point(0,-1,0), 1);
		Sphere s1 = new Sphere(new Point(10,-1,0), 9.2);
		Lens l = new Lens(s0,s1);
		j3ds.addShape(l, Color.ORANGE.darker());

		ParametricParaboloid surf = new ParametricParaboloid(0.1, 0.1, 0.01, new Vector(0,0,-2)); 
		j3ds.addSurface(surf, new Color(255,165,0));


		//		j3ds.getCamera().setLocation(new Point(1,0,2));
		//		j3ds.getCamera().setLocation(new Point(1,0,2));
		//		j3ds.getCamera().setFrontFog(1);

		//Small animation example
		double t=0;
		while(true){
			t+=0.05;
			Matrix m = Matrix.createRotationMatrix(t, new Vector(0,0,1));
			surf.setRotation(m);
			j3ds.repaint();
			try {
				Thread.sleep(50);
			} catch (InterruptedException e) {
			}
		}
	}

}
