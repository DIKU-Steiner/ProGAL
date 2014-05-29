package ProGAL.geom3d.viewer;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.util.ConcurrentModificationException;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;

import javax.media.j3d.LinearFog;
import javax.media.j3d.Transform3D;
import javax.media.j3d.View;
import javax.swing.ButtonGroup;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JLayeredPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import ProGAL.geom3d.Line;
import ProGAL.geom3d.Point;
import ProGAL.geom3d.Shape;
import ProGAL.geom3d.Triangle;
import ProGAL.geom3d.Vector;
import ProGAL.geom3d.volumes.LSS;
import ProGAL.geom3d.volumes.Sphere;
import ProGAL.geom3d.volumes.Tetrahedron;

import com.sun.j3d.utils.universe.ViewingPlatform;

/**
 * Represents a camera located somewhere in the scene. This object completely decides the viewpoint 
 * shown in the canvas. The camera is specified by the following properties
 * <ul>
 * <li>The position-point indicates the position of the eye-point.</li>
 * <li>The lookingAt-point indicates what the eye is looking at.</li>
 * <li>The up-vector indicates the up-direction of the screen.</li>
 * <li>The boolean perspective-parameter indicates if perspective or parallel projection is used</li>
 * <ul><li>If perspective-projection is enabled, the view-angle indicates the zoom level</li></ul>
 * <ul><li>If perspective-projection is disabled, the view-size indicates the zoom level</li></ul>
 * <li>The front- and back-clip variables indicate at which distances the scene is clipped.</li>
 * <li>The front- and back-fog variables indicate at which distances the scene starts and ends fog.</li>
 * </ul>
 */
public class Camera implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	private final J3DScene scene;
	private final ViewingPlatform platform;
	private final View view;
	private final LinearFog fog;
	private final Transform3D transform;
	protected boolean rotateAroundObject = true;

	private Point eye, lookingAt;
	private Vector up;

	Camera(J3DScene j3dScene, ViewingPlatform vp, LinearFog fog){
		this.scene = j3dScene;
		this.view = scene.canvas.getView();
		this.fog = fog;
		this.platform = vp;
		transform = new Transform3D();
		platform.getViewPlatformTransform().getTransform(transform);
		double[] entries = new double[16];
		transform.get(entries);

		eye = new Point(entries[3], entries[7], entries[11]);
		lookingAt = new Point(0,0,0);
		up = new Vector(entries[1], entries[5], entries[9]);

		this.scene.canvas.addMouseListener(this);
		this.scene.canvas.addMouseMotionListener(this);
		this.scene.canvas.addMouseWheelListener(this);
		this.scene.canvas.addKeyListener(this);

		setBackClip(10000);
		setBackFog(1000);
		setFrontFog(50);
	}

	public Point getEye(){ return eye; }
	public Point getLookingAt(){ return lookingAt; }
	public Vector getUp(){ return up; }
	public double getViewAngle(){ return view.getFieldOfView(); }
	public double getFrontClip(){ return view.getFrontClipDistance()*6.7; }
	public double getBackClip(){ return view.getBackClipDistance()*6.7; }
	public double getFrontFog(){ return fog.getFrontDistance(); }
	public double getBackFog(){ return fog.getBackDistance(); }
	public boolean getParallel(){ return (view.getProjectionPolicy()==View.PARALLEL_PROJECTION); }

	public void setLocation(Point p){ 
		this.eye = p;
		updateView();
	}
	public void setLookingAt(Point p){ 
		this.lookingAt = p;
		updateView();
	}
	public void setUp(Vector v){
		this.up = v;
		updateView();
	}
	public void setViewAngle(double a){		view.setFieldOfView(a);		}
	public void setFrontClip(double s){		view.setFrontClipDistance(s/6.7);	}
	public void setBackClip(double s){ 		view.setBackClipDistance(s/6.7);	}
	public void setFrontFog(double s){		fog.setFrontDistance(s);	}
	public void setBackFog(double s){ 		fog.setBackDistance(s);	}
	//	public void setParallel(boolean p){		
	//		if(p)	view.setProjectionPolicy(View.PARALLEL_PROJECTION);
	//		else 	view.setProjectionPolicy(View.PERSPECTIVE_PROJECTION);
	//	}

	public void updateView(){
		platform.getViewPlatformTransform().getTransform(transform);
		//Prepare camera rotation (determined from location.VectorTo(lookingAt) and up).
		Vector c = eye.vectorTo(lookingAt).multiplyThis(-1/eye.distance(lookingAt));
		Vector b = up.subtract( c.multiply(up.dot(c)) ).normalizeThis();
		Vector a = b.cross(c);

		transform.set(new double[]{
				a.x(),b.x(),c.x(),eye.x(), 
				a.y(),b.y(),c.y(),eye.y(),
				a.z(),b.z(),c.z(),eye.z(),
				0,0,0,1 });
		platform.getViewPlatformTransform().setTransform(transform);
		if(controlPanel!=null && controlPanel.isVisible()) {
			collectShapes();
			controlPanel.repaint();
		}
	}
	public void mouseWheelMoved(MouseWheelEvent e) {
		int rot = e.getWheelRotation();
		double dist = eye.distance(lookingAt);
//		Line l = new Line(eye, eye.vectorTo(lookingAt));
//		if(shiftPressed){
//			double t;
//			if(rot<0){
//				t = 1-Math.pow(0.99, -rot);
//			}else{
//				t = 1-Math.pow(1.0/0.99, rot);
//			}
//			location = l.getPoint(t);
//			lookingAt = l.getPoint(1+t);
//		}else{
		
			if(rot<0){
//				System.out.println(eye.distanceSquared(lookingAt));
				double lineVal = 1-Math.pow(0.99, -rot);
				if(dist<0.001) return;
//				eye = l.getPoint(lineVal);
				eye = eye.add(eye.vectorTo(lookingAt).multiplyThis(lineVal));
			}else{
				double lineVal = 1-Math.pow(1.0/0.95, rot);
				if(dist>1000) return;
//				eye = l.getPoint(lineVal);
				eye = eye.add(eye.vectorTo(lookingAt).multiplyThis(lineVal));
			}
//		}
		updateView();
	}
	public void mouseClicked(MouseEvent e) {}
	public void mouseEntered(MouseEvent e) {}
	public void mouseExited(MouseEvent e) {}

	private java.awt.Point lastPoint;
	public void mousePressed(MouseEvent e) {	lastPoint = e.getLocationOnScreen();	}
	public void mouseReleased(MouseEvent e) {	lastPoint = null;}
	public void mouseDragged(MouseEvent e) {
		boolean rotateAroundObject = this.rotateAroundObject^shiftPressed;
		java.awt.Point p = e.getLocationOnScreen();
		int xdiff = p.x-lastPoint.x;
		int ydiff = p.y-lastPoint.y;
		lastPoint = p;
		double xangle = -xdiff*(rotateAroundObject?0.01:0.003);
		double yangle = -ydiff*(rotateAroundObject?0.01:0.003);
		Vector dir = eye.vectorTo(lookingAt).normalizeThis();
		Vector up = this.up;
		Vector right = dir.cross(up);
		up.rotateIn(dir, xangle);
		right.rotateIn(dir, yangle);
		right.rotateIn(up, yangle);
		if(rotateAroundObject){
			Line l = new Line(lookingAt, dir.normalizeThis());
			eye = l.getPoint(-eye.distance(lookingAt));
		}else{
			Line l = new Line(eye, dir.normalizeThis());
			lookingAt = l.getPoint(eye.distance(lookingAt));
		}
		updateView();
	}
	public void mouseMoved(MouseEvent e) {}

	private boolean shiftPressed = false;
	public void keyPressed(KeyEvent e) {
		if(e.getKeyCode()==KeyEvent.VK_SHIFT) shiftPressed = true;
	}
	public void keyReleased(KeyEvent e) {
		if(e.getKeyCode()==KeyEvent.VK_SHIFT) {shiftPressed = false; return; }
		switch(e.getKeyCode()){
		case KeyEvent.VK_Z: scene.autoZoom();break;
		case KeyEvent.VK_C: scene.centerCamera();break;
		case KeyEvent.VK_R: scene.toggleRotation();break;
		case KeyEvent.VK_A: scene.setAntialiasing(true);break;
		}
	}
	public void keyTyped(KeyEvent e) {}

	private JFrame controlFrame;
	private ControlPanel controlPanel;
	public JFrame getControlPanel(){
		return controlFrame;
	}
	public void createControlPanel(){
		controlPanel = new ControlPanel();
		controlFrame = new ControlFrame();
	}

	private class ControlFrame extends JFrame{
		private static final long serialVersionUID = 1L;

		ControlFrame(){
			setSize(new Dimension(600,250));
			setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
			setAlwaysOnTop(true);
			getContentPane().add(controlPanel);
			addKeyListener(Camera.this);
			controlPanel.addKeyListener(Camera.this);
		}

	}

	private class ControlPanel extends JPanel{
		private static final long serialVersionUID = 1L;
		private JLayeredPane pane = new JLayeredPane();
		private final int panelWidth = 200;
		List<ProGAL.geom2d.Shape> shapes = new LinkedList<ProGAL.geom2d.Shape>(); 
		List<java.awt.Color> shapeColors = new LinkedList<java.awt.Color>(); 


		ControlPanel(){
			super(new GridLayout(1,1));
			add(pane);

			ButtonGroup group = new ButtonGroup();
			JRadioButton rb = new JRadioButton("Rotate around object",rotateAroundObject);
			pane.add(rb);
			rb.setLocation(10,10);
			rb.setSize(panelWidth-10,20);
			rb.setSelected(rotateAroundObject);
			group.add(rb);
			rb.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					rotateAroundObject = true;
				}
			});
			rb = new JRadioButton("Rotate around camera",!rotateAroundObject);
			pane.add(rb);
			rb.setLocation(10,30);
			rb.setSize(panelWidth-10,20);
			rb.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent arg0) {
					rotateAroundObject = false;
				}
			});
			group.add(rb);

			int spinnerWidth = 80;
			int spinnerHeight = 25;

			//Angle spinner
			JSpinner spinner = new JSpinner(new SpinnerNumberModel((int)(getViewAngle()*180/Math.PI), 0, 360, 1));
			pane.add(spinner);
			spinner.setSize(spinnerWidth,spinnerHeight);
			spinner.setLocation(10, 55+0*spinnerHeight);
			spinner.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent arg0) {
					int a = (Integer)((JSpinner)arg0.getSource()).getValue();
					setViewAngle(a*Math.PI/180);
					repaint();
				}
			});
			JLabel l = new JLabel("View angle");
			pane.add(l);
			l.setLocation(15+spinnerWidth, 55+0*spinnerHeight);
			l.setSize(panelWidth-spinnerWidth, 20);

			spinner = new JSpinner(new SpinnerNumberModel(getFrontClip(), 0, 1000, 0.1));
			pane.add(spinner);
			spinner.setSize(spinnerWidth,spinnerHeight);
			spinner.setLocation(10, 60+1*spinnerHeight);
			spinner.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent arg0) {
					double a = (Double)((JSpinner)arg0.getSource()).getValue();
					setFrontClip(a);
					repaint();
				}
			});
			l = new JLabel("Front clipping");
			pane.add(l);
			l.setLocation(15+spinnerWidth, 60+1*spinnerHeight);
			l.setSize(panelWidth-spinnerWidth, 20);

			spinner = new JSpinner(new SpinnerNumberModel(getBackClip(), 0, Math.max(getBackClip(), 10000), 0.1));
			pane.add(spinner);
			spinner.setSize(spinnerWidth,spinnerHeight);
			spinner.setLocation(10, 60+2*spinnerHeight);
			spinner.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent arg0) {
					double a = (Double)((JSpinner)arg0.getSource()).getValue();
					setBackClip(a);
					repaint();
				}
			});
			l = new JLabel("Back clipping");
			pane.add(l);
			l.setLocation(15+spinnerWidth, 60+2*spinnerHeight);
			l.setSize(panelWidth-spinnerWidth, 20);


			spinner = new JSpinner(new SpinnerNumberModel(getFrontFog(), 0, 1000, 0.01));
			pane.add(spinner);
			spinner.setSize(spinnerWidth,spinnerHeight);
			spinner.setLocation(10, 65+3*spinnerHeight);
			spinner.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent arg0) {
					double a = (Double)((JSpinner)arg0.getSource()).getValue();
					setFrontFog(a);
					repaint();
				}
			});
			l = new JLabel("Front fog");
			pane.add(l);
			l.setLocation(15+spinnerWidth, 65+3*spinnerHeight);
			l.setSize(panelWidth-spinnerWidth, 20);

			spinner = new JSpinner(new SpinnerNumberModel(getBackFog(), 0, Math.max(getBackFog(), 10000), 0.01));
			pane.add(spinner);
			spinner.setSize(spinnerWidth,spinnerHeight);
			spinner.setLocation(10, 65+4*spinnerHeight);
			spinner.addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent arg0) {
					double a = (Double)((JSpinner)arg0.getSource()).getValue();
					setBackFog(a);
					repaint();
				}
			});
			l = new JLabel("Back fog");
			pane.add(l);
			l.setLocation(15+spinnerWidth, 65+4*spinnerHeight);
			l.setSize(panelWidth-spinnerWidth, 20);



		}


		public void paintComponent(Graphics g){
			super.paintComponent(g);
			Graphics2D g2d = (Graphics2D)g;
			int w = getWidth();
			int h = getHeight();
			g2d.setColor(Color.BLACK);
			g2d.drawLine(panelWidth, 5, panelWidth, h-5);
			g2d.drawLine(panelWidth, 25, w-5, 25);

			g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			//Paint camera field
			int x = panelWidth+5, y = 30;
			w = w-panelWidth-10;
			h = h-35;
			java.awt.Point camPos = new java.awt.Point(x,y+h/2);
			java.awt.Point locPos = new java.awt.Point(x+w/2,y+h/2);


//			java.awt.Shape oldClip = g2d.getClip();
//			Area area = new Area(new Rectangle(new java.awt.Point(x,y), new Dimension(w,h)));
//			g2d.setClip(area);
			double realDist = eye.distance(lookingAt);
			double planeDist = camPos.distance(locPos);
			for(int i=0;i<shapes.size();i++){//TODO: Improve
				g2d.setColor(shapeColors.get(i));
				ProGAL.geom2d.Shape shape = shapes.get(i);
				if(shape instanceof ProGAL.geom2d.Circle){
					ProGAL.geom2d.Circle circ = (ProGAL.geom2d.Circle)shape;
					int diam = (int)(circ.getRadius()*2*planeDist/realDist);
					int xt = x+(int)(circ.center().x()*planeDist/realDist)-diam/2;
					int yt = y+h/2-(int)(circ.center().y()*planeDist/realDist)-diam/2;
					g2d.fillOval(xt, yt, diam, diam);
				}
				if(shape instanceof ProGAL.geom2d.LSC){
					ProGAL.geom2d.LSC lsc = (ProGAL.geom2d.LSC)shape;
					float diam = (float)(lsc.getRadius()*2*planeDist/realDist);
					int x1t = x+(int)(lsc.getSegment().getA().x()*planeDist/realDist);
					int y1t = y+h/2-(int)(lsc.getSegment().getA().y()*planeDist/realDist);
					int x2t = x+(int)(lsc.getSegment().getB().x()*planeDist/realDist);
					int y2t = y+h/2-(int)(lsc.getSegment().getB().y()*planeDist/realDist);
					Stroke old = g2d.getStroke();
					g2d.setStroke(new BasicStroke(diam, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND));
					g2d.drawLine(x1t, y1t, x2t, y2t);
					g2d.setStroke(old);
				}
				if(shape instanceof ProGAL.geom2d.Triangle){
					ProGAL.geom2d.Triangle tri = (ProGAL.geom2d.Triangle)shape;
					int x1t = x+(int)(tri.getCorner(0).x()*planeDist/realDist);
					int y1t = y+h/2-(int)(tri.getCorner(0).y()*planeDist/realDist);
					int x2t = x+(int)(tri.getCorner(1).x()*planeDist/realDist);
					int y2t = y+h/2-(int)(tri.getCorner(1).y()*planeDist/realDist);
					int x3t = x+(int)(tri.getCorner(2).x()*planeDist/realDist);
					int y3t = y+h/2-(int)(tri.getCorner(2).y()*planeDist/realDist);
					g2d.fillPolygon(new int[]{x1t,x2t,x3t}, new int[]{y1t,y2t,y3t}, 3);
				}
				
			}
//			g2d.setClip(oldClip);



			g2d.setColor(Color.BLACK);
			g2d.fillOval(camPos.x-2, camPos.y-2, 5, 5);
			g2d.fillOval(locPos.x-2, locPos.y-2, 5, 5);

			//Draw view cone
			int d = (int)(h/(2.0*Math.tan(getViewAngle()/2)));
			g2d.drawLine(camPos.x, camPos.y, camPos.x+d, y);
			g2d.drawLine(camPos.x, camPos.y, camPos.x+d, y+h);


			g2d.setColor(new Color(50,50,50));
			//Draw front clip
			int lineX = (int)( getFrontClip()*planeDist/realDist );
			g2d.drawLine(x+lineX, y, x+lineX, y+h);
			
			//Draw back clip
			lineX = (int)( getBackClip()*planeDist/realDist );
			g2d.drawLine(x+lineX, y, x+lineX, y+h);
			
			//Draw front fog
			
			//Draw back fog


		}
	}

	void collectShapes(){
		if(controlFrame==null || !controlFrame.isVisible()) return;

		Color background = controlPanel.getBackground();
		List<ProGAL.geom2d.Shape> newShapes = new LinkedList<ProGAL.geom2d.Shape>();
		List<java.awt.Color> newColors = new LinkedList<java.awt.Color>();
		Vector x = eye.vectorTo(lookingAt).normalizeThis();
		Vector y = up.subtract( x.multiply(up.dot(x)) ).normalizeThis();
		try{
			for(Entry<Shape,Color> entry: scene.primitives.entrySet()){
				int objectsAdded = 0;
				
				Shape shape = entry.getKey();
				if(shape instanceof Sphere){ 
					Sphere s = (Sphere)shape;
					double xT = x.dot(eye.vectorTo(s.getCenter()));
					double yT = y.dot(eye.vectorTo(s.getCenter()));
					newShapes.add(new ProGAL.geom2d.Circle(new ProGAL.geom2d.Point(xT,yT), s.getRadius()));
					objectsAdded = 1;
				}
				//Tets
				if(shape instanceof Tetrahedron){ 
					Tetrahedron s = (Tetrahedron)shape;
					double x1T = x.dot(eye.vectorTo(s.getPoint(0)));
					double y1T = y.dot(eye.vectorTo(s.getPoint(0)));
					double x2T = x.dot(eye.vectorTo(s.getPoint(1)));
					double y2T = y.dot(eye.vectorTo(s.getPoint(1)));
					double x3T = x.dot(eye.vectorTo(s.getPoint(2)));
					double y3T = y.dot(eye.vectorTo(s.getPoint(2)));
					double x4T = x.dot(eye.vectorTo(s.getPoint(3)));
					double y4T = y.dot(eye.vectorTo(s.getPoint(3)));
					newShapes.add(new ProGAL.geom2d.Triangle(new ProGAL.geom2d.Point(x1T,y1T),new ProGAL.geom2d.Point(x2T,y2T),new ProGAL.geom2d.Point(x3T,y3T) ));
					newShapes.add(new ProGAL.geom2d.Triangle(new ProGAL.geom2d.Point(x1T,y1T),new ProGAL.geom2d.Point(x2T,y2T),new ProGAL.geom2d.Point(x4T,y4T) ));
					newShapes.add(new ProGAL.geom2d.Triangle(new ProGAL.geom2d.Point(x1T,y1T),new ProGAL.geom2d.Point(x3T,y3T),new ProGAL.geom2d.Point(x4T,y4T) ));
					newShapes.add(new ProGAL.geom2d.Triangle(new ProGAL.geom2d.Point(x2T,y2T),new ProGAL.geom2d.Point(x3T,y3T),new ProGAL.geom2d.Point(x4T,y4T) ));
					objectsAdded = 4;
				}
				//Tris
				if(shape instanceof Triangle){ 
					Triangle s = (Triangle)shape;
					double x1T = x.dot(eye.vectorTo(s.getPoint(0)));
					double y1T = y.dot(eye.vectorTo(s.getPoint(0)));
					double x2T = x.dot(eye.vectorTo(s.getPoint(1)));
					double y2T = y.dot(eye.vectorTo(s.getPoint(1)));
					double x3T = x.dot(eye.vectorTo(s.getPoint(2)));
					double y3T = y.dot(eye.vectorTo(s.getPoint(2)));
					newShapes.add(new ProGAL.geom2d.Triangle(new ProGAL.geom2d.Point(x1T,y1T),new ProGAL.geom2d.Point(x2T,y2T),new ProGAL.geom2d.Point(x3T,y3T) ));
					objectsAdded = 1;
				}
				//LSSs
				if(shape instanceof LSS){ 
					LSS s = (LSS)shape;
					double x1T = x.dot(eye.vectorTo(s.segment.getA()));
					double y1T = y.dot(eye.vectorTo(s.segment.getA()));
					double x2T = x.dot(eye.vectorTo(s.segment.getB()));
					double y2T = y.dot(eye.vectorTo(s.segment.getB()));
					newShapes.add(new ProGAL.geom2d.LSC(new ProGAL.geom2d.Point(x1T,y1T),new ProGAL.geom2d.Point(x2T,y2T), s.rad));
					objectsAdded = 1;
				}

				for(int i=0;i<objectsAdded;i++){
					Color col = entry.getValue();
					double blend = 0.7;
					newColors.add(new Color(
							(int)(col.getRed()*(1-blend)+background.getRed()*blend), 
							(int)(col.getGreen()*(1-blend)+background.getGreen()*blend), 
							(int)(col.getBlue()*(1-blend)+background.getBlue()*blend) 
							));
				}

			}
			controlPanel.shapes = newShapes;
			controlPanel.shapeColors = newColors;
		}catch(ConcurrentModificationException exc){}
	}

}