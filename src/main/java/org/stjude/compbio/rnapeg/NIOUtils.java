package org.stjude.compbio.rnapeg;
// TO DO:
// ADD SANITY CHECKING FOR SOURCE/TARGET LENGTHS!

import java.nio.channels.*;
import java.io.*;

public class NIOUtils {

  public static void copyFile(File sourceFile, File destFile) throws IOException {
    // http://stackoverflow.com/questions/106770/standard-concise-way-to-copy-a-file-in-java
    // http://www.javalobby.org/java/forums/t17036.html
    if(!destFile.exists()) {
      destFile.createNewFile();
    }

    FileChannel source = null;
    FileChannel destination = null;

    try {
      source = new FileInputStream(sourceFile).getChannel();
      destination = new FileOutputStream(destFile).getChannel();
      destination.transferFrom(source, 0, source.size());
      // FAIL for files larger than signed int limit??
      // observed truncated files of size:
      //     2,147,479,552
      // signed int limit is:
      //     2,147,483,647
      //
      // Java API appears to specify longs, not sure what the problem is
    }
    finally {
      if(source != null) {
	source.close();
      }
      if(destination != null) {
	destination.close();
      }
    }
  }

}
