package main.java.xyz.experiment;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

/**
 * Created by guoy28 on 12/21/16.
 */
public class TestAnnotation {

  public static void main (String[] args) throws Exception {
    int tests = 0;
    int passed = 0;
    Class testClass = Class.forName("main.java.xyz.experiment.IhavesomeAnnotation");
    for (Method m : testClass.getDeclaredMethods()) {
      if (m.isAnnotationPresent(Test.class)) {
        tests++;
        try {
          m.invoke(null);
          passed++;
        } catch (InvocationTargetException e) {
          Throwable t = e.getCause();
          System.out.println("valid exception: " + t);
        } catch (Exception e) {
          System.out.println("wrong: " + e);
        }
      }
    }
    System.out.println("tests: " + tests + "passed " + passed);
  }
}
