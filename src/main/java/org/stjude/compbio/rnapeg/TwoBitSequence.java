package org.stjude.compbio.rnapeg;
// .2bit files: http://genome.ucsc.edu/FAQ/FAQformat#format7
// implementation of full sequence record
// MNE 2/2010

import java.io.*;
import java.util.*;

public class TwoBitSequence {

  private ArrayList<TwoBitBlock> n_blocks, mask_blocks;
  private int sequence_length;
  private EndianReader er;
  private long packed_dna_offset;

  private static final int BASES_PACKED_PER_BYTE = 4;
  private static final boolean VERBOSE = false;
  private boolean parse_error = false;

  public TwoBitSequence (EndianReader er, int offset, boolean read_fully) {
    this.er = er;
    setup(offset, read_fully);
  }

  public void set_randomaccessfile (RandomAccessFile raf) {
    // hacktacular
    er.set_randomaccessfile(raf);
  }

  private void setup (int offset, boolean read_fully) {
    //
    //  decode sequence header section
    //
    try {
      er.seek(offset);

      sequence_length = er.readInt();

      if (!read_fully) {
	//	System.err.println("TwoBitSequence: reading length only!");  // debug
	return;
      }

      //      System.err.println("N blocks:");  // debug
      n_blocks = read_block_set();
      if (VERBOSE) System.err.println("N block count=" + n_blocks.size());  // debug

      //      System.err.println("mask blocks:");  // debug
      mask_blocks = read_block_set();
      if (VERBOSE) System.err.println("mask block count=" + mask_blocks.size());  // debug

      int reserved = er.readInt();
      if (reserved != 0) {
	System.err.println("WARNING: nonzero reserved value");  // debug
      }

      packed_dna_offset = er.getFilePointer();
      //      BigEndianReader ber = new BigEndianReader(er.get_randomaccessfile());
      //      ber.seek(offset);
      //      System.err.println("slen2=" + ber.readInt());  // debug
      //      System.err.println("bcount2=" + ber.readInt());  // debug

    } catch (Exception e) {
      System.err.println("ERROR: " + e);  // debug
      e.printStackTrace();
    }
  }

  public byte[] get_region (int start_base, int region_size) throws IOException {
    //
    //  decode a portion of the sequence.
    //
    Timer t;
    // t = new Funk.Timer("sequence read");
    byte[] result = read_sequence(start_base, region_size);
    // raw sequence
    //    t.finish();

    //    t = new Funk.Timer("block apply");

    //
    //  Apply N regions:
    //
    if (n_blocks != null) {
      for (TwoBitBlock nb : n_blocks) {
	nb.set_n(result, start_base, region_size);
      }
    } else {
      System.err.println("no N blocks in .2bit file");  // debug
    }

    //
    //  apply masking:
    //
    if (mask_blocks != null) {
      for (TwoBitBlock mb : mask_blocks) {
	mb.mask(result, start_base, region_size);
      }
    } else {
      System.err.println("no mask blocks in .2bit file");  // debug
    }
    //    t.finish();

    return result;
  }

  private byte[] read_sequence (int start_base, int bases_needed) throws IOException {
    //
    //  parse packed DNA; no masking or N regions processed.
    //
    byte[] sequence_buf = new byte[bases_needed];
    //    Arrays.fill(sequence_buf, (byte) 'N');
    RandomAccessFile raf = er.get_randomaccessfile();

    //    new Exception().printStackTrace();

    //    System.err.println("valid handle? " + raf.getFD().valid());  // debug

    long seek_start = packed_dna_offset + ((start_base - 1) / BASES_PACKED_PER_BYTE);
    raf.seek(seek_start);

    int start_block_base = (start_base - 1) % BASES_PACKED_PER_BYTE;

    //    System.err.println("start_base: " + start_base + " seek: " + seek_start + " SBB: " + start_block_base);  // debug

    int block;
    //
    // not sure if bases are packed in an architecture-dependent fashion or not;
    // assume ordering as shown in spec example.
    //
    int bases_wanted_in_block, shift, value;
    int seq_ptr = 0;
    byte base;

    BufferedInputStream bis = new BufferedInputStream(new FileInputStream(raf.getFD()));
    // buffer for (much) better read performance

    boolean first_block = true;

    //
    //  decode bases:
    //
    parse_error = false;

    while (bases_needed > 0) {
      block = bis.read();
      if (block == -1) {
	System.err.println("ERROR: still need " + bases_needed + " bases, but EOF!");  // debug
	parse_error = true;
	break;
      }
      if (first_block) {
	bases_wanted_in_block = BASES_PACKED_PER_BYTE - start_block_base;
	shift = (bases_wanted_in_block - 1) * 2;
	//	System.err.println("start wanted " + bases_wanted_in_block + " shift=" + shift);  // debug
      } else {
	bases_wanted_in_block = BASES_PACKED_PER_BYTE;
	shift = 6;
	// hack: adapt from constant?
	// => no, this is slower for subsequent loops
      }

      if (bases_wanted_in_block > bases_needed) bases_wanted_in_block = bases_needed;

      for (; bases_wanted_in_block > 0;
	   bases_wanted_in_block--, shift -= 2, bases_needed--) {
	//	  value = (block >> shift) & 0x03;
	//	  switch (value) {
	switch ((block >> shift) & 0x03) {
	case 0: base = 'T'; break;
	case 1: base = 'C'; break;
	case 2: base = 'A'; break;
	case 3: base = 'G'; break;
	default: base = '!'; 
	}
	sequence_buf[seq_ptr++] = base;
	//	  System.err.println("shift=" + shift + " value=" + ((block >> shift) & 0x03) + " base="+ base);  // debug

	//	  if (seq_ptr % 10000000 == 0) System.err.println("ptr="+seq_ptr);  // debug
      }
      first_block = false;
    }
    //      System.err.println("done");  // debug
    
    return parse_error ? null : sequence_buf;
  }

  public byte[] get_full_sequence () throws IOException {
    //
    //  decode and return entire sequence.
    //  Slow and memory intensive.
    //
    byte[] full_sequence = read_sequence(1, sequence_length);

    //
    //  Apply N regions:
    //
    
    if (n_blocks != null) {
      for (TwoBitBlock nb : n_blocks) {
	nb.set_n(full_sequence);
      }
    } else {
      System.err.println("no N blocks in .2bit file");  // debug
    }

    //
    //  apply masking:
    //
    if (mask_blocks != null) {
      for (TwoBitBlock mb : mask_blocks) {
	mb.mask(full_sequence);
      }
    } else {
      System.err.println("no mask blocks in .2bit file");  // debug
    }

    //      System.err.println(new String(full_sequence, 34950, 2000));  // debug
    return full_sequence;
  }

  private ArrayList<TwoBitBlock> read_block_set() throws IOException {
    int block_count = er.readInt();
    int[] starts = read_int_set(block_count);
    int[] sizes = read_int_set(block_count);

    ArrayList<TwoBitBlock> set = new ArrayList<TwoBitBlock>();
    for (int i = 0; i < block_count; i++) {
      TwoBitBlock block = new TwoBitBlock(starts[i], sizes[i]);
      //      System.err.println("  start="+block.start + " size="+block.size);  // debug
      set.add(block);
    }
    return set;
  }

  private int[] read_int_set (int count) throws IOException {
    int[] set = new int[count];
    for (int i = 0; i < count; i++) {
      set[i] = er.readInt();
    }
    return set;
  }

  public int get_length() {
    return sequence_length;
  }
			 

}
