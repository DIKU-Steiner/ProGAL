package ProGAL.dataStructures.viewer;

import java.awt.Color;

public interface InteractiveBinaryTree {
	public InteractiveBinaryTree left();
	public InteractiveBinaryTree right();
	public Color leftLegColor();
	public Color rightLegColor();
	public Color nodeColor();
	public String label();
	public void click();
	
}
