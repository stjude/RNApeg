package org.stjude.compbio.rnapeg;

import java.io.*;

public class FASTAFileTools {

  private long fa_ptr_id;
  // file pointer for sequence ID line

  private FAIIndexRecord index;

  //  private long fa_ptr_sequence;
  // file pointer for start of sequence

  //  private long fa_cooked_line_nt;
  // number of nucleotides per line
  //  private long fa_raw_line_bytes;
  // number of file bytes per line
  //  private int line_overhead;

  //  private long sequence_length;
  // length (in nt) of sequence (assuming 1 sequence in file)

  public FAIIndexRecord detect_record_properties (RandomAccessFile raf) throws IOException {
    // given a RandomAccessFile pointed to start of a FASTA record (">" line),
    // set class variables for
    //
    // - file pointer of ID line
    // - file pointer of first sequence line
    // - number of nucleotides per line (cooked)
    // - number of file bytes per line (includes line break)
    //
    index = new FAIIndexRecord();

    fa_ptr_id = raf.getFilePointer();
    String id_line = raf.readLine();
    if (id_line.indexOf(">") == 0) {
      index.file_sequence_offset = raf.getFilePointer();
      String first_sequence_line = raf.readLine();
      long ptr_second = raf.getFilePointer();
      index.bytes_per_line = (int) (ptr_second - index.file_sequence_offset);
      //      System.err.println("first="+first_sequence_line);  // debug
      index.nt_per_line = first_sequence_line.length();
      //      System.err.println("nt/line="+index.nt_per_line);  // debug
      //      System.err.println("bytes/line="+index.bytes_per_line);  // debug

      long flen = raf.length();
      //      System.err.println("flen="+flen);  // debug

      long body_size = flen - index.file_sequence_offset;
      //      System.err.println("body="+body_size);  // debug

      long lines = body_size / index.bytes_per_line;
      int line_overhead = index.get_line_overhead();

      long leftover = body_size % index.bytes_per_line;
      if (leftover > 0) leftover -= line_overhead;
      // if partial line, make sure to account for line-terminator overhead

      index.sequence_length = (int) ((lines * index.nt_per_line) + leftover);
      // only valid if 1 sequence in this file!
    } else {
      System.err.println("ERROR: expected FASTA ID line, got " + id_line);  // debug
    }
    return index;
  }

  public FAIIndexRecord get_index_info() {
    return index;
  }

  public byte[] get_region (RandomAccessFile raf, int start_base, int length) throws IOException {
    byte[] results = null;

    if (start_base > index.sequence_length) {
      System.err.println("ERROR: start base past end of sequence");  // debug
    } else if (start_base < 1) {
      System.err.println("ERROR: start base must be at least 1");  // debug
    } else {
      //	System.err.println("FIX ME: TRIM if asking for too much");  // debug

      int overflow = (start_base + length - 1) - index.sequence_length;
      //	System.err.println("overflow="+overflow);  // debug
      if (overflow > 0) {
	System.err.println("WARNING: requested region past end of sequence");  // debug
	length -= overflow;
      }

      int needed = length;

      int base_index = start_base - 1;
      // convert from base # to sequence index

      int start_line_num = base_index / index.nt_per_line;
      int start_line_index = base_index % index.nt_per_line;

      long start_ptr = index.file_sequence_offset +
	(start_line_num * index.bytes_per_line) +
	start_line_index;

      raf.seek(start_ptr);

      String line;
      results = new byte[length];
      int buf_ptr = 0;
      int i;

      while (true) {
	line = raf.readLine();
	if (line == null) break;
	char[] chars = line.toCharArray();
	// slow, but a PITA to deal with line endings if in middle of line
	for (i=0; i < chars.length && needed > 0; needed--) {
	  results[buf_ptr++] = (byte) chars[i++];
	}
	if (needed <= 0) break;
      }
    }

    return results;
  }

  public byte[] get_all (RandomAccessFile raf) throws IOException {
    int needed = index.sequence_length;
    byte[] results = new byte[needed];
    int buf_ptr = 0;

    raf.seek(index.file_sequence_offset);
    // move file pointer to start of sequence

    int read;
      
    while (true) {
      if (needed <= 0) break;

      try {
	//	  System.err.println("needed="+needed + " lb="+ line_bytes);  // debug
	read = needed < index.bytes_per_line ? needed : index.bytes_per_line;
	raf.readFully(results, buf_ptr, read);
	buf_ptr += index.nt_per_line;
	// only move pointer up by the number of nt in the line
	// (i.e. not the line break at the end)
	needed -= index.nt_per_line;
	//	  System.err.println("buf_ptr="+buf_ptr);  // debug
	//	  System.err.println("needed="+needed);  // debug
      } catch (EOFException ex) {
	//	  System.err.println("EOF");  // debug
	break;
      }
    }

    return results;
  }

  public void set_index (FAIIndexRecord index) {
    this.index = index;
  }



}