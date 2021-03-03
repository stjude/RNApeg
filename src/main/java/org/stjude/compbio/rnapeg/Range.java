package org.stjude.compbio.rnapeg;

public class Range {
  public int start, end;

  public Range() {
    start = end = -1;
  }
  
  public Range (int start, int end) {
    this.start = start;
    this.end = end;
  }

  public boolean isValid() {
    return start != -1 && end != -1 && end >= start;
  }

  public String toString() {
    return start + "-" + end;
  }

  public int size() {
    return (end - start) + 1;
  }

}
