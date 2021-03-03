package org.stjude.compbio.rnapeg;
// simple minimum/maximum value tracker

import java.awt.Rectangle;

public class RangeDigest implements Comparable {
  public int start, end;
  // public for speed
  private boolean is_first_start, is_first_end;

  public RangeDigest() {
    setup();
  }

  public RangeDigest(int start, int end) {
    setup();
    add_start(start);
    add_end(end);
  }

  public RangeDigest(TranscriptFeature tf) {
    setup();
    add_start(tf.start);
    add_end(tf.end);
  }

  // begin Comparable interface
  public int compareTo(Object o) {
    RangeDigest rd_other = (RangeDigest) o;
    int result;
    if (rd_other.start == start) {
      result = 0;
    } else if (start < rd_other.start) {
      result = -1;
    } else {
      result = 1;
    }
    return result;
  }
  // end Comparable interface


  private void setup() {
    is_first_start = is_first_end = true;
  }

  public void add (Rectangle r) {
    add_start(r.x);
    add_end(r.x + (r.width - 1));
  }

  public void add (RangeDigest rd) {
    add_start(rd.get_minimum());
    add_end(rd.get_maximum());
  }

  public void add_start (int x) {
    if (is_first_start || x < start) {
      start = x;
      is_first_start = false;
    }
  }

  public void add_end (int x) {
    if (is_first_end || x > end) {
      end = x;
      is_first_end = false;
    }
  }

  public void add_start_and_length (int x, int len) {
    add_start(x);
    add_end(x + (len - 1));
  }

  public void add_start_and_end (int start, int end) {
    add_start(start);
    add_end(end);
  }
  
  public int get_minimum () {
    return start;
  }

  public int get_maximum () {
    return end;
  }

  public int get_size() {
    return (end - start) + 1;
  }

  public boolean intersects (int other_start, int other_end) {
    return !(other_end < start || other_start > end);
  }

  public boolean intersects (RangeDigest other) {
    return intersects(other.start, other.end);
  }

  public boolean intersects (Rectangle r) {
    return intersects(r.x, r.x + (r.width - 1));
  }

  public boolean intersects (int i) {
    return intersects(i, i);
  }

  public boolean is_adjacent (RangeDigest other) {
    return other.end == start - 1 ||
      other.start == end + 1;
  }

  public String toString() {
    return start + "-" + end;
  }

  public boolean equals (int start, int end) {
    return this.start == start && this.end == end;
  }

  public int get_overlap (RangeDigest other) {
    int overlap = 0;
    if (intersects(other)) {
      Rectangle r_this = generate_rectangle(this);
      Rectangle r_other = generate_rectangle(other);
      Rectangle r_overlap = r_this.intersection(r_other);
      overlap = r_overlap.width;
    }
    return overlap;
  }

  private Rectangle generate_rectangle(RangeDigest rd) {
    return new Rectangle(rd.start,
			 1,
			 (rd.end - rd.start) + 1,
			 1);
  }


}
