package org.stjude.compbio.rnapeg;

import java.util.*;

public class HashCounter {
  // FIX ME: generics version

  HashMap<String,Integer> counts;

  public HashCounter() {
    counts = new HashMap<String,Integer>();
  }

  public int add (String s) {
    Integer count = counts.get(s);
    if (count == null) count = Integer.valueOf(0);
    count++;
    counts.put(s, count);
    return count;
  }

  public void add (char c) {
    add(Character.toString(c));
  }

  public void add_ignore_case (char c) {
    add(Character.toString(Character.toUpperCase(c)));
  }

  public HashMap<String,Integer> get_counts() {
    return counts;
  }

}