package org.stjude.compbio.rnapeg;

import java.io.*;
import java.util.*;

public interface ReferenceSequence {
  public static final int NULL_LENGTH = -1;

  public byte[] get_region (String sequence_name, int start_base, int length) throws IOException;
  // fetch a region of a reference sequence
  // start_base is base NUMBER (i.e. starts with 1), NOT 0-based index

  public byte[] get_all (String sequence_name) throws IOException;
  // fetch entire reference sequence

  public int get_length (String sequence_name) throws IOException;
  // sequence length

  public boolean supports_sequence_list();
  public ArrayList<String> get_sequence_names();

}
