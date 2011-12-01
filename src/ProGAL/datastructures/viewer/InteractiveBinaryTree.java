package ProGAL.datastructures.viewer;

import java.awt.Color;

import ProGAL.datastructures.BinaryTree;

public interface InteractiveBinaryTree extends BinaryTree {
	public Color leftLegColor();
	public Color rightLegColor();
	public Color nodeColor();
	public String label();
	public void click();
	
}
