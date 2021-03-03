package org.stjude.compbio.rnapeg;
// information about sequence flanking a splice event.
// each splice has 2 copies, 1 for upstream sequence and 1 for downstream.

public class SplicedReadFlankingInfo {
  int count_aligned_bases = 0;
  // could be either bases before the start/end of the sequence, 
  // or bases until the next splice.
  // Might be a good idea to split this into 2 categories:
  // for bases before a read start/end we probably want to be
  // fairly strict, because short runs might indicate a 
  // lower-confidence mapping.
  // Cases where a single read contains multiple junctions
  // are a tougher call because some isoforms contain very small
  // exons and we don't want to reject them outright.

  int count_aligned_mismatched_bases = 0;
  // mismatched bases in aligned read blocks only

  int count_soft_clip_bases = 0;
  int count_inserted_bases = 0;
  int count_deleted_bases = 0;

  int last_read_i = 0;
  // the rightmost index of aligned sequence:
  //
  //   SplicedReadInfo_left         intron           SplicedReadInfo_right
  // ACGTACGTACGTCGTACGTACGTAC-------------------->AGCGTACGTACGTACGTACGTACGTA
  //        last_read_i----->|                     |<----index 0
  //
  // these positions can be used to determine whether the read 
  // has mismatches at the splice site, a QC metric which we might 
  // be able to use to help resolve ambiguous novel junctions later.

  public double get_junk_ratio() {
    if (false) {
      System.err.println("DEBUG, limiting");  // debug
      return (double) count_aligned_mismatched_bases / count_aligned_bases;
    } else {
      return (double) (
		       count_aligned_mismatched_bases +
		       count_soft_clip_bases +
		       count_inserted_bases +
		       count_deleted_bases
		       ) / count_aligned_bases;
    }
  }

  public double get_junk_ratio_soft_clip() {
    return (double) count_soft_clip_bases / count_aligned_bases;
  }

  public double get_junk_ratio_indel() {
    return (double) (count_inserted_bases + count_deleted_bases) / count_aligned_bases;
  }

  public double get_junk_ratio_mismatches() {
    return (double) count_aligned_mismatched_bases / count_aligned_bases;
  }

}
