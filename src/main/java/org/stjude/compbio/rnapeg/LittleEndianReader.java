package org.stjude.compbio.rnapeg;

import java.io.*;
import java.util.*;

public class LittleEndianReader extends EndianReader {

  private byte[] int_buf = new byte[4];

  public LittleEndianReader (RandomAccessFile raf) {
    super(raf);
  }

  public int readInt() throws IOException {
    int result = 0;
    if (raf.read(int_buf) == 4) {
      //      result = int_buf[0] + (int_buf[1] << 8) + (int_buf[2] << 16) + (int_buf[3] << 24);
      result = (int_buf[0] & 0xff) +
	((int_buf[1] & 0xff) << 8) +
	((int_buf[2] & 0xff) << 16) +
	((int_buf[3] & 0xff) << 24);
      // masking required because java bytes are signed!

      //      System.err.println("buf=" + int_buf[0] + " " + int_buf[1] + " " + int_buf[2] + " " + int_buf[3] + "; result="+result);  // debug


    } else {
      throw new IOException("insufficient bytes to read 32-bit int");
    }
    return result;
  }


}