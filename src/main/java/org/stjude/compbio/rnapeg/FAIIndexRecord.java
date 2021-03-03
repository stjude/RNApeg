package org.stjude.compbio.rnapeg;

public class FAIIndexRecord {
  String sequence_id;
  int sequence_length;
  long file_sequence_offset;
  int nt_per_line;
  // count of nucleotides per line (cooked)
  int bytes_per_line;
  // count of file bytes per line (raw, includes line breaks)

  public FAIIndexRecord () {
    // manual
  }

  public FAIIndexRecord (String s) {
    // .fai file line
    String[] fields = s.split("\t");
    if (fields.length == 5) {
      sequence_id = new String(fields[0]);
      sequence_length = Integer.parseInt(fields[1]);
      file_sequence_offset = Long.parseLong(fields[2]);
      nt_per_line = Integer.parseInt(fields[3]);
      bytes_per_line = Integer.parseInt(fields[4]);
    } else {
      System.err.println("ERROR: need 5-field FAI record");  // debug
    }
  }

  public int get_line_overhead () {
    return (int) (bytes_per_line - nt_per_line);
  }

}