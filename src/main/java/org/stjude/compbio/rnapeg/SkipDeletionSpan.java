package org.stjude.compbio.rnapeg;
// for a given CIGAR skip, find ultimate upstream and downstream
// boundaries, considering deletions and additional skips
// MNE 9/2026

import htsjdk.samtools.*;
import java.io.*;
import java.util.*;

public class SkipDeletionSpan {
  private ArrayList<CigarElement> ces;
  private int ces_max_lookahead_i;
  // max index where we are able to look ahead to the next element
  private int ces_size;
  private SAMRecord sr;
  private int seq1_end, seq2_start;
  private boolean usable;
  
   
  public SkipDeletionSpan (ArrayList<CigarElement> ces, SAMRecord sr) {
    this.ces = ces;
    this.sr = sr;
    ces_size = ces.size();
    ces_max_lookahead_i = ces_size - 2;
  }

  public boolean find_boundaries (int ref_base, int cei) {
    // TO DO: return boolean?

    usable = true;
    seq1_end = boundary_check(ref_base, cei, -1);
    seq2_start = boundary_check(ref_base, cei, 1);

    return usable;
  }

  private int boundary_check (int ref_base, int cei, int direction) {
    // start index is of a SKIPPED_REGION operator
    int pos = ref_base;

    if (direction == -1) {
      pos = ref_base - 1;
      // base # before skip starts
    } else if (direction == 1) {
      CigarElement ce = ces.get(cei);
      pos = ref_base + ce.getLength();
      // base # after skip ends
    } else {
      System.err.println("ERROR, illegal direction " + direction);  // debug
      System.exit(1);
    }

    int extra_skips = 0;
    while (true) {
      cei += direction;
      // move to next neighboring operator in specified direction
      
      if (direction == -1) {
	if (cei < 0) break;
      } else if (direction == 1) {
	if (cei >= ces_size) break;
      } else {
	System.err.println("ERROR, illegal direction " + direction);  // debug
	System.exit(1);
      }

      // OK
      CigarElement ce_adj = ces.get(cei);
      CigarOperator co = ce_adj.getOperator();
      
      if (co.equals(CigarOperator.SKIPPED_REGION)) {
	extra_skips++;
	pos += ce_adj.getLength() * direction;
      } else if (co.equals(CigarOperator.DELETION)) {
	pos += ce_adj.getLength() * direction;
      } else if (co.equals(CigarOperator.MATCH_OR_MISMATCH)) {
	break;
      }
      // need to handle other operators?
    }

    if (extra_skips > 0) {
      //      System.err.println("read with multiple skips: " + sr.getReadName());  // debug

      if (direction == -1) {
	// in the case of multiple skips, only allow use of the first one
	usable = false;
      }
    }

    return pos;
  }

  public int get_start() {
    return seq1_end;
  }

  public int get_end() {
    return seq2_start;
  }

}


  
