package BilevelBOBILIB_DISCRETE;

/**
 * Helper object used to sort variable indices by a score.
 */
public class Item implements Comparable<Item>{

	int pos;
	int val;
	
	public Item(int p, int v) {
		pos = p;
		val = v;
	}
	
	/**
	 * Orders items in descending score.
	 */
	@Override
	public int compareTo(Item a) {
		return -Integer.compare(val, a.val);
	}
	
}
