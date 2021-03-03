package org.stjude.compbio.rnapeg;

public class Timer {
  private String name;
  private long start_time;

  public Timer (String name) {
    this.name = name;
    start_time = System.currentTimeMillis();
  }

  public void finish () {
    System.err.println(name + " took " + 
		       (System.currentTimeMillis() - start_time) + " ms");
  }
}
