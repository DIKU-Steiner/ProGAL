package ProGAL.geom2d.viewer;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import ProGAL.geom2d.Circle;
import ProGAL.geom2d.LineSegment;
import ProGAL.geom2d.Point;
import ProGAL.geom2d.Shape;

public class J2DScene {
	public JFrame frame;
	private PaintPanel canvasPanel;
	private List<ShapeOptions> shapes = new ArrayList<ShapeOptions>();
	private Point camCenter;
	double scale;
	
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
			}
		});
		canvasPanel.setFocusable(true);
	}

	public JPanel getCanvas(){
		return canvasPanel;
	}
	
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
	
	public void addShape(Shape s, Color c){
		addShape(s,c,0.01,false);
	}
	public void addShape(Shape s, Color c, double border){
		addShape(s,c,border,false);
	}
	public void addShape(Shape s, Color c, double border, boolean fill){
		shapes.add(new ShapeOptions(s,c,border,fill));
		repaint();
	}
	public void removeShape(Shape s) {
		shapes.remove(s);
		repaint();
	}

	
	public void repaint(){
		canvasPanel.repaint();
	}
	
	private final static ShapePainter[] shapePainters = {
		new CirclePainter(),
		new LineSegmentPainter()
	};
	private static ShapePainter getShapePainter(Shape s){
		if(s instanceof Circle) return shapePainters[0];
		if(s instanceof LineSegment) return shapePainters[1];
			return null;
	}

	
	
	private class PaintPanel extends JPanel{
		private static final long serialVersionUID = 1L;

		public void paint(Graphics g){
			Graphics2D g2d = (Graphics2D)g;
			g2d.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, 	RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
			g2d.setRenderingHint(RenderingHints.KEY_RENDERING, 			RenderingHints.VALUE_RENDER_QUALITY);
			g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 		RenderingHints.VALUE_ANTIALIAS_ON);
			g2d.setColor(Color.WHITE);
			g.fillRect(0, 0, getWidth(), getHeight());
			for(ShapeOptions so: new LinkedList<ShapeOptions>(shapes)){
				Stroke oldStroke = g2d.getStroke();
				g2d.setStroke(  new BasicStroke(
						(float)(scale*so.borderWidth),
						BasicStroke.CAP_ROUND, 
						BasicStroke.JOIN_MITER)  );
				g2d.setColor(so.color);
				getShapePainter(so.shape).paintShape(so, g2d);
				g2d.setStroke(oldStroke);
			}
		}
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
			int w = canvasPanel.getWidth();
			int h = canvasPanel.getHeight();
			int gX = (int)(scale*(p.x()-camCenter.x())+w/2);
			int gY = (int)(-scale*(p.y()-camCenter.y())+h/2);
			return new java.awt.Point(gX,gY);
		}
		
		double getScale(){ return scale; }
	}
	
	
	public static J2DScene createJ2DSceneInFrame(){
		J2DScene scene = new J2DScene();
		JFrame frame = new JFrame("J2DScene");
		frame.getContentPane().add(scene.canvasPanel);
		frame.setSize(1000,1000);
		frame.setVisible(true);
		return scene;
	}
	
	public static void main(String[] args){
		J2DScene scene = J2DScene.createJ2DSceneInFrame();
		scene.addShape(new LineSegment(new Point(0,0), new Point(1,0)), Color.BLACK);
		scene.addShape(new LineSegment(new Point(0,0), new Point(0,1)), Color.BLACK);
		scene.addShape(new Circle(new Point(0,0), 1), new Color(250,0,0,100), 0.05, false);
		scene.addShape(new Circle(new Point(1,0), 0.2), Color.GREEN);
		scene.addShape(new Circle(new Point(0,1), 0.2), Color.BLUE);
	}

}
