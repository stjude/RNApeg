package org.stjude.compbio.rnapeg;

public class SplicedReadInfo implements Comparable {
  String reference_name;
  Integer segment_1_end, segment_2_start;
  // reference base numbers (1-based)
  int counter, counter_plus, counter_minus;
  // count of supporting reads, and reads with support on + and - strands
  int counter_perfect, counter_clean, counter_flanking_qc;

  public SplicedReadInfo(String reference_name, int segment_1_end, int segment_2_start) {
    this.reference_name = reference_name;
    this.segment_1_end = segment_1_end;
    this.segment_2_start = segment_2_start;
    counter = counter_plus = counter_minus = 0;
    counter_perfect = counter_clean = 0;
    counter_flanking_qc = 0;
  }

  public String get_name() {
    return reference_name + ":" + segment_1_end + "-" + segment_2_start;
  }
  public void increment_counter() {
    counter++;
  }

  public void increment_plus() {
    counter_plus++;
  }

  public void increment_minus() {
    counter_minus++;
  }

  public void increment_flanking_qc() {
    counter_flanking_qc++;
  }

  public int get_count() {
    return counter;
  }

  // begin Comparable interface
  public int compareTo(Object other) {
    SplicedReadInfo o = (SplicedReadInfo) other;
    int diff = reference_name.compareTo(o.reference_name);
    if (diff == 0) {
      diff = segment_1_end.compareTo(o.segment_1_end);
      if (diff == 0) diff = segment_2_start.compareTo(o.segment_2_start);
    }
    return diff;
  }
  // end Comparable interface

  public boolean overlaps (SplicedReadInfo other) {
    return !(other.segment_2_start < segment_1_end ||
	     other.segment_1_end > segment_2_start);
  }

  public int get_span_size() {
    return segment_2_start - segment_1_end - 1;
  }

}