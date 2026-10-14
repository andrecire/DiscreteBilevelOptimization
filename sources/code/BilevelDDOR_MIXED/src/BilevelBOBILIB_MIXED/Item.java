package BilevelBOBILIB_MIXED;

/**
 * Helper object used to sort variable indices by a score.
 */
public class Item implements Comparable<Item>{

	int pos;
	double val;
	
	public Item(int p, double v) {
		pos = p;
		val = v;
	}
	
	/**
	 * Orders items in descending score.
	 */
	@Override
	public int compareTo(Item a) {
		return -Double.compare(val, a.val);
	}
	
}
