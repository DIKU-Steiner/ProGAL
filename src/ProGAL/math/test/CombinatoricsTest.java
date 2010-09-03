package ProGAL.math.test;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;

import org.junit.Test;

import ProGAL.math.Combinatorics;

public class CombinatoricsTest {

	@Test
	public void testGetAllPermutations() {
		Throwable e=null;
		//Should give empty-permutation if max<=0
		try{ Combinatorics.getAllPermutations(-5); }catch(Exception ex){ e=ex; }
		assertTrue(e instanceof RuntimeException);

		assertTrue(listContains(Combinatorics.getAllPermutations(0),new int[]{}));
		assertEquals(1,Combinatorics.getAllPermutations(0).size());
		

		List<int[]> perms = Combinatorics.getAllPermutations(3);
		assertTrue(listContains(perms,new int[]{0,1,2}));
		assertTrue(listContains(perms,new int[]{0,2,1}));
		assertTrue(listContains(perms,new int[]{1,0,2}));
		assertTrue(listContains(perms,new int[]{1,2,0}));
		assertTrue(listContains(perms,new int[]{2,0,1}));
		assertTrue(listContains(perms,new int[]{2,1,0}));
		
		//And just a touch of stress-testing as well
		perms = Combinatorics.getAllPermutations(4);
		assertTrue(listContains(perms,new int[]{2,1,0,3}));
		assertTrue(listContains(perms,new int[]{3,2,1,0}));
		assertTrue(listContains(perms,new int[]{0,1,2,3}));
	}
	
	private boolean listContains(List<int[]> list, int[] elem){
		for(int[] el: list)
			if(Arrays.equals(el, elem)) return true;
		return false;
 	}

	@Test
	public void testBinom() {
		//Cover meaningless cases (n<0 & m>n)
		assertEquals(0,Combinatorics.binom(10, 11));
		assertEquals(0,Combinatorics.binom(9, 11));
		assertEquals(0,Combinatorics.binom(-5, -2));
		assertEquals(0,Combinatorics.binom(5, -2));
		assertEquals(Combinatorics.binom(5, 2),Combinatorics.binom(-5, 2));
		assertEquals(Combinatorics.binom(5, 0),Combinatorics.binom(-5, 0));
		
		//Cover n=0
		assertEquals(1,Combinatorics.binom(0, 0));
		assertEquals(0,Combinatorics.binom(0, 1));
		assertEquals(0,Combinatorics.binom(0, 3));
		assertEquals(0,Combinatorics.binom(0, -3));

		//Cover n>0 & m<=n
		assertEquals(1,Combinatorics.binom(6, 6));
		assertEquals(6,Combinatorics.binom(6, 5));
		assertEquals(20,Combinatorics.binom(6, 3));
		assertEquals(1,Combinatorics.binom(6, 0));
	}

}
