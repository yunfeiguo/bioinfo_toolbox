package experiment;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Created by guoy28 on 12/6/16.
 */
public class BinaryIO {
  public static void main(String[] args) {
    try (DataOutputStream lenOut = new DataOutputStream(new BufferedOutputStream(new FileOutputStream("/tmp/1.tmp")))) {
    lenOut.writeInt(123);
      lenOut.close();
    } catch (IOException e) {
      e.printStackTrace();
    }

  }
}
