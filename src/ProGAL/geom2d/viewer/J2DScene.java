package ProGAL.geom2d.viewer;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JPanel;

import ProGAL.geom2d.*;

/**
 * A graphics class for viewing 2D scenes. 
 * Many of the <code>Shape</code>-subclasses specified in the <code>ProGAL.geom2d</code> 
 * package can be added to a <code>J2DScene</code> object and are automatically 
 * painted. For instance the following code creates a scene with 3 circles, two lines 
 * indicating unit axis, and a text label. 
 * <pre>
 * J2DScene scene = new J2DScene();
 * scene.addShape(new LineSegment(new Point(0,0), new Point(1,0)), Color.BLACK);
 * scene.addShape(new LineSegment(new Point(0,0), new Point(0,1)), Color.BLACK);
 * scene.addShape(new Circle(new Point(0,0), 1), new Color(250,0,0,100), 0.05, false);
 * scene.addShape(new Circle(new Point(1,0), 0.2), Color.GREEN);
 * scene.addShape(new Circle(new Point(0,1), 0.2), Color.BLUE);
 * scene.addShape(new TextShape("(1,0)", new Point(1,0), 0.2), Color.BLACK);
 * 
 * JPanel canvas = scene.getCanvas();
 * 
 * JFrame frame = new JFrame();
 * frame.setSize(400,400);
 * frame.getContentPane().add( canvas );
 * frame.setVisible(true);
 * </pre>
 * 
 * The <code>repaint()</code> method must be called every time the position of 
 * shapes has changed and the canvas should be updated. The pointers 
 * to added shapes are stored, so subsequent changes will be visible on the canvas when 
 * <code>repaint()</code> is called. The following example shows how to animate a circle 
 * rotating around origo.
 * <pre>
 * J2DScene scene = J2DScene.createJ2DSceneInFrame();
 * scene.addShape(new LineSegment(new Point(0,0), new Point(1,0)), Color.BLACK);
 * scene.addShape(new LineSegment(new Point(0,0), new Point(0,1)), Color.BLACK);
 * Circle c = new Circle(new Point(0,1), 0.2);
 * scene.addShape(c, Color.BLUE);
 * 
 * double t = 0;
 * while(true){
 * 	c.center().setCoord(0, Math.cos(t));
 * 	c.center().setCoord(1, Math.sin(t));
 * 	scene.repaint();
 * 	t+=0.01;
 * 	try {Thread.sleep(50);} catch (InterruptedException e) {}
 * }
 * </pre>
 * As shown in this example, a static method is supplied for conveniently creating a frame 
 * containing a scene-viewer.
 * 
 * Shapes are painted in the order they are inserted. The <code>addShape</code> method has 
 * support for border width and filling out a shape. 
 * 
 * @author R. Fonseca
 */
public class J2DScene {
	public JFrame frame;
	private final PaintPanel canvasPanel;
	private final List<ClickListener> clickListeners = new LinkedList<ClickListener>();
	private final List<ShapeOptions> shapes = new ArrayList<ShapeOptions>();
	private Point camCenter;
	double scale;

	/**
	 * Construct a representation of a 2D scene. 
	 */
	public J2DScene(){
		canvasPanel = new PaintPanel();
		camCenter = new Point(0,0);
		scale = 100;

		canvasPanel.addKeyListener(new KeyListener() {
			public void keyTyped(KeyEvent e) {}
			public void keyReleased(KeyEvent e) {}
			public void keyPressed(KeyEvent e) {
				if(e.getKeyCode()==KeyEvent.VK_C)	centerCamera();
				if(e.getKeyCode()==KeyEvent.VK_Z)	autoZoom();
				if(e.getKeyCode()==KeyEvent.VK_S)	savePng("J2DSnapshot.png");
			}
		});
		canvasPanel.addMouseListener(canvasPanel);
		canvasPanel.addMouseMotionListener(canvasPanel);
		canvasPanel.addMouseWheelListener(canvasPanel);
		canvasPanel.setFocusable(true);
	}

	/** Return the panel that this scene is painted on. */
	public JPanel getCanvas(){
		return canvasPanel;
	}

	/** Center the view on the objects in the scene */
	public void centerCamera(){
		double minX=Double.POSITIVE_INFINITY;
		double maxX=Double.NEGATIVE_INFINITY;
		double minY=Double.POSITIVE_INFINITY;
		double maxY=Double.NEGATIVE_INFINITY;
		for(ShapeOptions so: shapes){
			Point p = so.shape.getCenter();
			if(p.x()<minX) minX = p.x();
			if(p.y()<minY) minY = p.y();
			if(p.x()>maxX) maxX = p.x();
			if(p.y()>maxY) maxY = p.y();
		}
		camCenter = new Point((minX+maxX)/2, (minY+maxY)/2);
		repaint();
	}

	/** Zoom the view to enclose all the objects in the scene. */
	public void autoZoom(){
		double minX=Double.POSITIVE_INFINITY;
		double maxX=Double.NEGATIVE_INFINITY;
		double minY=Double.POSITIVE_INFINITY;
		double maxY=Double.NEGATIVE_INFINITY;
		for(ShapeOptions so: shapes){
			Point p = so.shape.getCenter();
			if(p.x()<minX) minX = p.x();
			if(p.y()<minY) minY = p.y();
			if(p.x()>maxX) maxX = p.x();
			if(p.y()>maxY) maxY = p.y();
		}
		int w = canvasPanel.getWidth();
		int h = canvasPanel.getHeight();
		if(w==0) w = 1000;
		if(h==0) h = 700;
		scale = Math.min(w/(maxX-minX), h/(maxY-minY))*0.9;
		repaint();
	}

	/** Add a shape to this scene. Currently, 
	 * <ul><li><code>ProGAL.geom2d.Circle</code></li>
	 * <li><code>ProGAL.geom2d.LineSegment</code> and </li>
	 * <li><code>ProGAL.geom2d.viewer.TextShape</code></li>
	 * <li><code>ProGAL.geom2d.LSC</code></li></ul>
	 * are supported. 
	 */
	public void addShape(Shape s){
		addShape(s,Color.GRAY,0.01,false);
	}

	/** Add a shape to this scene with the specified color. Currently, 
	 * <ul><li><code>ProGAL.geom2d.Circle</code></li>
	 * <li><code>ProGAL.geom2d.LineSegment</code> and </li>
	 * <li><code>ProGAL.geom2d.viewer.TextShape</code></li>
	 * <li><code>ProGAL.geom2d.LSC</code></li></ul>
	 * are supported. 
	 */
	public void addShape(Shape s, Color c){
		addShape(s,c,0.01,false);
	}
	public void addShape(Shape s, Color c, double border){
		addShape(s,c,border,false);
	}
	public void addShape(Shape s, Color c, double border, boolean fill){
		shapes.add(new ShapeOptions(s,c,border,fill));
		//		repaint();
	}


	public void savePng(String fName)
	{
		
		BufferedImage bImg = new BufferedImage(canvasPanel.getWidth(), canvasPanel.getHeight(), BufferedImage.TYPE_INT_RGB);
		Graphics2D cg = bImg.createGraphics();
		canvasPanel.paintAll(cg);
		try {
			if (ImageIO.write(bImg, "png", new File(fName)))
			{
				System.out.println("Successfully wrote "+fName);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	/** Remove the specified shape from the scene */
	public void removeShape(Shape s) {
		ShapeOptions opt = null;
		for(ShapeOptions so: shapes) if(so.shape==s){ opt = so; break; }
		shapes.remove(opt);
		//		repaint();
	}

	/** Remove all shapes from the scene */
	public void removeAllShapes() {
		//		while (!shapes.isEmpty()) shapes.remove(0);
		shapes.clear();
		//		repaint();
	}

	/** Add a click-listener that gets called every time an object or the background is clicked	 */
	public void addClickListener(ClickListener cl){
		clickListeners.add(cl);
	}


	/** Repaint the scene */
	public void repaint(){
		//		canvasPanel.paintImmediately(0, 0, canvasPanel.getWidth(), canvasPanel.getHeight());
		canvasPanel.repaint();
	}

	private final static ShapePainter[] shapePainters = {
		new CirclePainter(),
		new LineSegmentPainter(),
		new TextPainter(),
		new LSCPainter(),
		new TrianglePainter(),
		new PolygonPainter(),
		new LinePainter()
	};
	private static ShapePainter getShapePainter(Shape s){
		if(s instanceof Circle) return shapePainters[0];
		if(s instanceof LineSegment) return shapePainters[1];
		if(s instanceof TextShape) return shapePainters[2];
		if(s instanceof LSC) return shapePainters[3];
		if(s instanceof Triangle) return shapePainters[4];
		if(s instanceof Polygon) return shapePainters[5];
		if(s instanceof Line) return shapePainters[6];
		return null;
	}



	private class PaintPanel extends JPanel implements MouseListener, MouseMotionListener, MouseWheelListener{
		private static final long serialVersionUID = 1L;

		public void paint(Graphics g){
			Graphics2D g2d = (Graphics2D)g;
			g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, 	RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
			g2d.setRenderingHint(RenderingHints.KEY_RENDERING, 			RenderingHints.VALUE_RENDER_QUALITY);
			g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 		RenderingHints.VALUE_ANTIALIAS_ON);
			g2d.setColor(Color.WHITE);
			g.fillRect(0, 0, getWidth(), getHeight());
			for(ShapeOptions so: new LinkedList<ShapeOptions>(shapes)){
				if(so==null) continue;
				Stroke oldStroke = g2d.getStroke();
				g2d.setStroke(  new BasicStroke(
						(float)(scale*so.borderWidth),
						BasicStroke.CAP_ROUND, 
						BasicStroke.JOIN_ROUND)  );
//				g2d.setColor(so.color);
				g2d.setPaint(so.color);
				ShapePainter sp = getShapePainter(so.shape);
				if(sp==null) System.err.println("J2DScene: ShapePainter not implemented for "+so.shape.getClass().getSimpleName());
				else sp.paintShape(so, g2d);
				g2d.setStroke(oldStroke);
			}
		}

		java.awt.Point lastPoint = null;
		public void mouseDragged(MouseEvent e) {
			java.awt.Point p = e.getLocationOnScreen();
			//			java.awt.Point dP = new java.awt.Point(p.x-lastPoint.x, p.y-lastPoint.y);
			Point p0 = transformPoint(lastPoint);
			Point p1 = transformPoint(p);
			camCenter.addThis(new Vector(p0.x()-p1.x(), p0.y()-p1.y()));
			repaint();
			lastPoint = p;
		}

		Point transformPoint(java.awt.Point p){
			int w = getWidth();
			int h = getHeight();
			double pX = (p.x-w/2)/scale + camCenter.x();
			double pY = -(p.y-h/2)/scale + camCenter.y();
			return new Point(pX,pY);
		}

		public void mouseMoved(MouseEvent e) {	}
		public void mouseClicked(MouseEvent e) {
			Point p = null;
			Shape shapeClicked = null;

			for(ShapeOptions so: new LinkedList<ShapeOptions>(shapes)){
				if(so.fill){
					if(p==null)	p = transformPoint(e.getPoint());
					if(so.shape.contains(p))
						shapeClicked = so.shape;
				}
			}

			for(ClickListener cl: clickListeners){
				cl.shapeClicked(shapeClicked, e);
			}
		}
		public void mouseEntered(MouseEvent e) {	}
		public void mouseExited(MouseEvent e) {		}
		public void mousePressed(MouseEvent e) { lastPoint = e.getLocationOnScreen(); }
		public void mouseReleased(MouseEvent e) {	}
		public void mouseWheelMoved(MouseWheelEvent e) {
			double factor =Math.pow(0.97, e.getWheelRotation()); 
			scale*=factor;
			repaint();
		}
	}


	public java.awt.Point transformPoint(Point p){
		int w = canvasPanel.getWidth();
		int h = canvasPanel.getHeight();
		int gX = (int)(scale*(p.x()-camCenter.x())+w/2);
		int gY = (int)(-scale*(p.y()-camCenter.y())+h/2);
		return new java.awt.Point(gX,gY);
	}

	public Point transformPoint(java.awt.Point p){
		int w = canvasPanel.getWidth();
		int h = canvasPanel.getHeight();
		int gX = p.x;
		int gY = p.y;
		double x = ((gX-w/2)/scale)+camCenter.x();
		double y = ((gY-h/2)/-scale)+camCenter.y(); 
		return new Point(x,y);
	}

	class ShapeOptions{
		Shape shape;
		Color color;
		double borderWidth;
		boolean fill;

		ShapeOptions(Shape s, Color c, double bw, boolean f){
			this.shape = s;
			this.color = c;
			this.borderWidth = bw;
			this.fill = f;
		}

		java.awt.Point transformPoint(Point p){
			return J2DScene.this.transformPoint(p);
		}

		Point transformPoint(java.awt.Point p){
			return J2DScene.this.transformPoint(p);
		}


		double getScale(){ return scale; }
	}

	/** 
	 * Create a frame containing a canvas, display it and return the J2DScene object shown in the frame. 
	 * The frame can be retrieved using the <code>J2DScene.frame</code> field.  
	 */
	public static J2DScene createJ2DSceneInFrame(){
		J2DScene scene = new J2DScene();
		JFrame frame = new JFrame("J2DScene");
		frame.getContentPane().add(scene.canvasPanel);
		frame.setSize(1000,1000);
		frame.setVisible(true);
		scene.frame = frame;
		return scene;
	}


	public static void main(String[] args){
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		
		scene.addShape(new LineSegment(new Point(0,0), new Point(1,0)), Color.BLACK);
		scene.addShape(new LineSegment(new Point(0,0), new Point(0,1)), Color.BLACK);
		scene.addShape( new TextShape("(1,0)", new Point(1,0), 0.2) );
		scene.addShape( new TextShape("(0,1)", new Point(0,1), 0.2) );

		scene.addShape(new Circle(new Point(0.1,0.1), 0.01), Color.BLACK, 0, true);
		scene.addShape(new Circle(new Point(0.2,0.2), 0.01), Color.BLACK, 0, true);

		scene.addShape(new LSC(new Point(0,-0.5), new Point(1,-0.4), 0.2), Color.RED);
		scene.addShape(new LSC(new Point(0,-0.9), new Point(1,-0.9), 0.1), Color.RED, 0, true);

		Point p = new Point(-2,0);
		scene.addShape(new Triangle(new Point(-1,0), new Point(-1,1), p), Color.GREEN, 0, true);

		scene.addShape(new Polygon(
				new Point[]{
						new Point(-1,2), 
						new Point(-1,2.4), 
						new Point(-1.4,2.1),
						new Point(-1.6,2.8),
				}
				), Color.BLUE,0,true);

		scene.addShape(new Line(new Point(3,0), new Vector(0.01,-0.03)));

		scene.centerCamera();

		while(true){
			p.addThis(new Vector(0.01,0));
			scene.repaint();
			try {
				Thread.sleep(30);
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}


	}

}
