package org.stjude.compbio.rnapeg;

public class TwoBitBlock {
  int start, size, end;

  public TwoBitBlock (int start, int size) {
    this.start = start;
    this.size = size;
    this.end = start + size;
  }

  public void set_n (byte[] full_sequence) {
    int i, end;
    for (i = start, end = start + size; i < end; i++) {
      full_sequence[i] = 'N';
    }
  }

  public void mask (byte[] full_sequence) {
    int i, end;
    for (i = start, end = start + size; i < end; i++) {
      full_sequence[i] = (byte) Character.toLowerCase((char) full_sequence[i]);
    }
  }

  public void mask (byte[] region_buf, int start_base, int region_size) {
    int si = start_base - 1;
    int ei = (si + region_size);
    int i;

    if (!(end < si || start > ei)) {
      // block overlaps unless:
      //  - block ends before region
      //  - block starts after region
      
      int effect_start = start - si;
      if (effect_start < 0) effect_start = 0;
      int effect_end = end - si;
      //      if (effect_end >= ei) effect_end = ei;
      if (effect_end >= region_buf.length) effect_end = region_buf.length;
//       System.err.println("si="+ si + " ei=" + ei);  // debug
//       System.err.println("ee = " + effect_end);  // debug
//       System.err.println("mask block from " + start + " => " + end);  // debug
//       System.err.println("effect = " + effect_start + " " + effect_end);  // debug

      for (i=effect_start; i < effect_end; i++) {
	region_buf[i] = (byte) Character.toLowerCase((char) region_buf[i]);
      }

    }
    
  }
  
  public void set_n (byte[] region_buf, int start_base, int region_size) {
    int si = start_base - 1;
    int ei = (si + region_size);
    int i;

    if (!(end < si || start > ei)) {
      int effect_start = start - si;
      if (effect_start < 0) effect_start = 0;
      int effect_end = end - si;
      if (effect_end >= region_buf.length) effect_end = region_buf.length;
      for (i=effect_start; i < effect_end; i++) {
	region_buf[i] = 'N';
      }

    }
    
  }


}